package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"
)

type config struct {
	hPath      string
	sPath      string
	maxThreads int
}

func setup(args []string) *config {
	config := &config{}
	flag.StringVar(&config.hPath, "haps", "", "The pattern of haps files")
	flag.StringVar(&config.sPath, "sample", "", "The pattern of samples files")
	flag.IntVar(&config.maxThreads, "threads", 1, "Num threads")
	// allows args to be mocked https://github.com/nwjlyons/email/blob/master/inputs.go
	// can only run 1 such test, else, redefined flags error
	a := os.Args[1:]
	if args != nil {
		a = args
	}

	flag.CommandLine.Parse(a)

	return config
}

func init() {
	log.SetFlags(0)
}

func main() {
	config := setup(nil)

	hapsFiles, err := filepath.Glob(config.hPath)

	if err != nil {
		log.Fatal(err)
	}

	sampleFiles, err := filepath.Glob(config.sPath)

	if err != nil {
		log.Fatal(err)
	}

	if len(hapsFiles) != len(sampleFiles) {
		log.Fatal("Haps/Samples lengths don't match")
	}

	// They both start with chr... ; so get them to pair up properly
	// in case glob returns different orders
	sort.Strings(hapsFiles)
	sort.Strings(sampleFiles)

	for idx, hPath := range hapsFiles {
		wmFh, err := os.Create(strconv.Itoa(idx) + ".map")

		if err != nil {
			log.Fatal(err)
		}

		wpFh, err := os.Create(strconv.Itoa(idx) + ".ped")

		if err != nil {
			log.Fatal(err)
		}

		wMap := bufio.NewWriter(wmFh)
		wPed := bufio.NewWriter(wpFh)

		sPath := sampleFiles[idx]

		samples, err := getSamples(sPath)

		if err != nil {
			log.Fatal(err)
		}

		if len(samples) == 0 {
			log.Fatalf("No samples found for file %s", sPath)
		}

		wFunc, completer := createWriteFunc(wMap, wPed, idx, samples)

		fmt.Println("Reading", hPath, sPath)
		read(config, hPath, sPath, wFunc)

		wMap.Flush()
		wPed.Flush()

		completer()
	}
}

func createWriteFunc(mWriter *bufio.Writer, pWriter *bufio.Writer, baseChr int, samples []string) (accumulator func(chrHaps map[string][]uint8, subChrNum int), completer func()) {
	log.Println("CALLED WITH THIS MANY SAMPLES", len(samples), baseChr)

	chr := baseChr

	var id int
	// var hIdx int
	// var sPed []string
	ped := make([]bytes.Buffer, len(samples), len(samples))
	// Keep track of fake chrs, over sum{1 .. n}(numSubChr)

	for _, famLine := range samples {
		var f bytes.Buffer

		f.WriteString(famLine)

		ped = append(ped, f)
	}

	completer = func() {
		for _, lBuffer := range ped {
			lBuffer.WriteString("\n")
			pWriter.WriteString(lBuffer.String())
		}
	}

	accumulator = func(chrHaps map[string][]uint8, subChrNum int) {
		for _, hSamples := range chrHaps {
			// the id of the haplotype; since we're not going to write
			// the full haplotype to the map file
			id++

			if len(hSamples) != len(samples) {
				log.Fatal("WTFF ", len(hSamples), len(samples), baseChr)
			}

			for idx, genotype := range hSamples {
				f := ped[idx]

				if genotype == 2 {
					f.WriteString(" 1 1")
				} else if genotype == 1 {
					f.WriteString(" 1 0")
				} else if genotype == 1 {
					f.WriteString(" 0 0")
				}
			}

			// Write the haplotype to the map file
			// place into func
			var l bytes.Buffer
			l.WriteString(strconv.Itoa(chr))
			l.WriteString("\t")
			l.WriteString(strconv.Itoa(id))
			l.WriteString("\t")
			l.WriteString(string(id))
			// Can we use numerical id as position? arbitrary right?
			l.WriteString("\t")
			l.WriteString("0")
			// Can we use numerical id as position? arbitrary right?
			l.WriteString("\t")
			l.WriteString(string(id))
			l.WriteString("\n")

			mWriter.WriteString(l.String())

			chr++
		}
	}

	return accumulator, completer
}

func read(config *config, hPath string, sPath string, resultFunc func(r map[string][]uint8, n int)) {
	reader, err := getReader(hPath)

	if err != nil {
		log.Fatal(err)
	}

	// store the genotypes tensor shape [ [ [], [] ] ]
	// sample: [geno1 array, geno2 array]
	// we use buffered channels
	// https://stackoverflow.com/questions/45366954/golang-channel-wierd-deadlock-related-with-make?rq=1
	genotypeQueue := make(chan [][][]string, 100)

	// The haplotypes, with sample counts, for a given sub chr
	results := make(chan map[string][]uint8, 100)
	// var wg sync.WaitGroup

	var wg sync.WaitGroup

	wg.Add(config.maxThreads)
	for i := 0; i < config.maxThreads; i++ {
		go processGenotypes(genotypeQueue, results, &wg)
	}

	// this could have come before processGenotypes
	go readFile(reader, genotypeQueue)

	// https://play.golang.org/p/nU2Rq6PJdO
	go func() {
		// wait for each thread to signal completion, then close the channel
		wg.Wait()
		close(results)
	}()

	var chr int
	for chrHaps := range results {
		// fmt.Println(chrHaps)
		resultFunc(chrHaps, chr)
		chr++
	}
}

func getSamples(sPath string) ([]string, error) {
	var err error
	reader, err := getReader(sPath)

	var samples []string
	if err != nil {
		return samples, err
	}

	for {
		row, err := reader.ReadString('\n') // 0x0A separator = newline

		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatal(err)
		} else if row == "" {
			// We may have not closed the pipe, but not have any more information to send
			// Wait for EOF
			continue
		}

		samples = append(samples, row[:len(row)-1])
	}

	return samples, err
}

func getReader(fPath string) (*bufio.Reader, error) {
	var err error
	var reader *bufio.Reader

	inFh, err := os.Open(fPath)

	if err != nil {
		return reader, err
	}

	// It's a relatively difficult thing to guess compression from header
	// probably some packages help
	if fPath[len(fPath)-3:] == ".gz" {
		gz, err := gzip.NewReader(inFh)

		if err != nil {
			return reader, err
		}

		reader = bufio.NewReader(gz)
	} else {
		reader = bufio.NewReader(inFh)
	}

	return reader, err
}

func readFile(reader *bufio.Reader, genotypeQueue chan [][][]string) {
	// How big the haplotype blocks are
	blockSize := int(1e6)
	// Tensor shape: [ [ [], [] ] ]
	// l1: sample
	// l2,l3: the two haplotypes, in R = 2 x n_sites_in_haplotype
	var genotypes [][][]string

	// for each subchr, for each row, we append an array containing that row's first 5 fields
	var prefixes [][][]string

	// Read the lines into the work queue.
	var lastStart int
	var subChrCutoff int
	subChr := -1
	for {
		row, err := reader.ReadString('\n') // 0x0A separator = newline

		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatal(err)
		} else if row == "" {
			// We may have not closed the pipe, but not have any more information to send
			// Wait for EOF
			continue
		}

		record := strings.Fields(row[:len(row)-1])
		pos, err := strconv.Atoi(record[2])

		if err != nil {
			log.Fatal(err)
		}

		if pos > subChrCutoff {
			fmt.Println("About to make new shit", subChr, lastStart, pos)

			lastStart = pos
			subChrCutoff = lastStart + blockSize
			subChr++

			// process the genotypes from this sub chromosome into haplotypes
			if len(genotypes) > 0 {
				genotypeQueue <- genotypes
			}

			// clear genotypes, start over
			genotypes = [][][]string{}
			prefixes = append(prefixes, [][]string{})

			fmt.Println("Made new shit", subChr, lastStart, pos)
		}

		prefixes[subChr] = append(prefixes[subChr], record[:5])

		sampleIdx := -1
		for i := 5; i < len(record); i += 2 {
			sampleIdx++

			if len(genotypes) <= sampleIdx {
				genotypes = append(genotypes, [][]string{[]string{}, []string{}})
			}

			genotypes[sampleIdx][0] = append(genotypes[sampleIdx][0], record[i])
			genotypes[sampleIdx][1] = append(genotypes[sampleIdx][1], record[i+1])
		}
	}

	if len(genotypes) > 0 {
		genotypeQueue <- genotypes
	}

	close(genotypeQueue)
}

func processGenotypes(gQueue chan [][][]string, results chan<- map[string][]uint8, wg *sync.WaitGroup) {
	defer wg.Done()
	haps := make(map[string][]uint8)

	for genotypes := range gQueue {
		for sIdx, sample := range genotypes {
			// fmt.Println("GENOTYPE LENGTH", len(genotypes))
			nSamples := len(genotypes)

			for _, hap := range makeHaps(sample[0]) {
				if haps[hap] == nil {
					haps[hap] = make([]uint8, nSamples, nSamples)
				}

				haps[hap][sIdx]++
			}

			for _, hap := range makeHaps(sample[1]) {
				if haps[hap] == nil {
					haps[hap] = make([]uint8, nSamples, nSamples)
				}

				haps[hap][sIdx]++
			}
		}

		results <- haps
	}
}

func makeHaps(genotype []string) []string {
	var haps []string

	full := strings.Join(genotype, "")

	// fmt.Println("Real length is", len(full))
	haps = append(haps, full)

	// From left to right replace letters with '-'
	// and in inner loop, from right to left do the same
	lastIdx := len(full) - 1
	breakIdx := len(full) - 2
	breakReverseIdx := len(full) - 3
	for idx := range full {
		if idx == 0 {
			full = "-" + full[idx+1:]

			// fmt.Println(full)
			haps = append(haps, full)
			continue
		}

		if idx == breakIdx {
			break
		}

		// for idx 0: -111 == "-" + full[0 + 1:]
		// for idx 1: "-" + full[:len(genotype)-1]
		// == full[:0+1] == "-" , + "-" + full[2:]
		full = full[:idx] + "-" + full[idx+1:]
		haps = append(haps, full)

		revH := full
		if idx == breakReverseIdx {
			break
		}

		fmt.Println(full)

		// initial()
		for y := lastIdx; y > idx+2; y-- {
			if y == lastIdx {
				revH = revH[:y] + "-"
				haps = append(haps, revH)

				fmt.Println(revH)
				continue
			}

			revH = revH[:y] + "-" + revH[y+1:]
			haps = append(haps, revH)

			fmt.Println(revH)
		}
	}

	return haps
}
