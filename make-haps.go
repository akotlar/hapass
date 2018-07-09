package main

// This is Dave Cutler's idea, as a Go package. Really fucking cool compression based on the mutation/recombination events
// over the history of each sample, which will happen to be a function of 4Ne

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/signal"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"syscall"
	// "unsafe"
)

type config struct {
	hPath      string
	sPath      string
	maxThreads int
	cpuProfile string
}

type node struct {
	sync.RWMutex
	left   *node
	right  *node
	parent *node
	val    uint8
}

type allele struct {
	chr string

	// The final combination; will walk up the tree using parentNode to actually get the full word, to avoid storing it
	// SO FUCKING COOL!
	// In the case of a single position "haplotype" (start == stop), this is the allele of 2N chromosomes (N diploid samples)
	// In the case of a multiple position haplotype (start < stop), this is whether sample N had this haplotype (1, else 0))
	finalNode *node
}

// track updates to tree
var updates uint64

func setup(args []string) *config {
	config := &config{}
	flag.StringVar(&config.hPath, "haps", "", "The pattern of haps files")
	flag.StringVar(&config.sPath, "sample", "", "The pattern of samples files")
	flag.IntVar(&config.maxThreads, "threads", runtime.GOMAXPROCS(0), "Num threads")
	flag.StringVar(&config.cpuProfile, "cprof", "", "Cpu profile path")

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

	if config.cpuProfile != "" {
		f, err := os.Create(config.cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	if config.maxThreads > runtime.GOMAXPROCS(0) {
		runtime.GOMAXPROCS(config.maxThreads)
	}

	// https://stackoverflow.com/questions/11268943/is-it-possible-to-capture-a-ctrlc-signal-and-run-a-cleanup-function-in-a-defe
	c := make(chan os.Signal)
	signal.Notify(c, os.Interrupt, syscall.SIGTERM)

	cleanup := func() {
		log.Println("\n\nStopped early!\n\n")
		pprof.StopCPUProfile()
	}

	go func() {
		<-c
		cleanup()
		os.Exit(1)
	}()

	run(config)
}

func run(config *config) {
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

	// This is our root node
	// It allows us to re-use our tree
	// Can be re-used across all chromosomes,
	// even with different numbers of samples per chr
	// since evolution of chromosomes is similar
	root := new(node)
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
		read(config, hPath, sPath, root, wFunc)

		wMap.Flush()
		wPed.Flush()

		completer()
	}
}

func createWriteFunc(mWriter *bufio.Writer, pWriter *bufio.Writer, baseChr int, samples []string) (accumulator func(words []*allele), completer func()) {
	log.Println("CALLED WITH THIS MANY SAMPLES", len(samples), baseChr)

	// chr := baseChr

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

	var id int
	accumulator = func(alts []*allele) {
		for _, alt := range alts {
			// log.Println("Called accumulator", alt)

			id++

			id := strconv.Itoa(id)
			// Write the haplotype to the map file
			// place into func
			var l bytes.Buffer
			l.WriteString(alt.chr)
			l.WriteString("\t")
			l.WriteString(id)
			l.WriteString("\t")
			// centimorgans... dummy
			l.WriteString("0")
			l.WriteString("\t")
			// position is basically a position; the positions aren't real
			// they're the order of the haplotypes we generate
			l.WriteString(id)

			wordNode := alt.finalNode
			for {
				l.WriteString("\t")

				if wordNode.val == 0 {
					l.WriteString("0")
				} else {
					l.WriteString("1")
				}

				if wordNode.parent == nil {
					break
				}

				wordNode = wordNode.parent
			}

			mWriter.WriteString(l.String())
		}
	}

	return accumulator, completer
}

func read(config *config, hPath string, sPath string, btree *node, resultFunc func(alts []*allele)) {
	reader, err := getReader(hPath)

	if err != nil {
		log.Fatal(err)
	}

	// May be better to encapsulate genotypeQueue into readFile
	// store the genotypes tensor shape [ [ [], [] ] ]
	// sample: [geno1 array, geno2 array]
	// we use buffered channels
	// https://stackoverflow.com/questions/45366954/golang-channel-wierd-deadlock-related-with-make?rq=1
	genotypeQueue := make(chan [][]string, 100)

	// The haplotypes, with sample counts, for a given sub chr
	results := make(chan []*allele, 100)

	var wg sync.WaitGroup

	wg.Add(config.maxThreads)
	for i := 0; i < config.maxThreads; i++ {
		go processGenotypes(genotypeQueue, btree, results, &wg)
	}

	// this could have come before processGenotypes
	go readFile(reader, genotypeQueue)

	// https://play.golang.org/p/nU2Rq6PJdO
	go func() {
		// wait for each thread to signal completion, then close the channel
		wg.Wait()
		close(results)
	}()

	for chrHaps := range results {
		// fmt.Println(chrHaps)
		resultFunc(chrHaps)
	}
}

// read the fam file
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

// Read a file and make haplotype blocks
func readFile(reader *bufio.Reader, genotypeQueue chan [][]string) {
	// How big the haplotype blocks are
	blockSize := int(1e6)

	// rows
	var records [][]string

	// Read the lines into the work queue.
	var lastStart int
	var subChrCutoff int

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

		if pos >= subChrCutoff {
			fmt.Println("About to make new", record[0], lastStart, pos, subChrCutoff, len(records))

			lastStart = pos
			subChrCutoff = lastStart + blockSize

			// process the genotypes from this sub chromosome into haplotypes
			if len(records) > 0 {
				log.Println("Sending this many sites: ", len(records))
				genotypeQueue <- records
			}

			records = [][]string{}
		}

		records = append(records, record)
	}

	if len(records) > 0 {
		fmt.Println("Sending remaining sites: ", lastStart, subChrCutoff, len(records))
		genotypeQueue <- records
	}

	close(genotypeQueue)
}

// Process the haplotype blocks into haplotypes that are defined
// by slices of each sample's column vector over some range within that block
// then build a binary tree (shared across all threads) that diffs
// against every other sample's vector slice, tracking diffs to avoid overwork
func processGenotypes(genotypeQueue chan [][]string, root *node, results chan<- []*allele, wg *sync.WaitGroup) {
	defer wg.Done()
	// A new item in the queue is in fact the entire cluster, block, of <= blockSize (default 1e6)
	// Genotypes contains a list of rows, of 2N + 5 length, exact length as .tped file
	// The first 4 columns are chr, id, pos, allele sequence,
	// int because more costly to convert int16 or uint16 to string; cannot use strconv.Itoa or strconf.FormatInt
	for genotypes := range genotypeQueue {
		// lSiteIdx := len(genotypes) - 1
		nSites := len(genotypes)

		log.Println("Received this many sites: ", nSites)
		// log.Println(genotypes[0])
		nSamples := len(genotypes[0]) - 5

		// these are the 2N samples, i.e all chromosomes
		// which matches nicely to the .tped format of 2 columns per sample, in .fam order
		sampleGenos := make([][]byte, nSamples, nSamples)

		// Just build a vector of sites for
		for i := range sampleGenos {
			sampleGenos[i] = make([]byte, nSites, nSites)
		}

		var chr string
		var hasShared bool
		var hasDiff bool
		var curr *node
		var alleles []*allele
		for rowIdx, row := range genotypes {
			chr = row[0]

			//Build the vector for each sample, since haplotypes consist of this vector
			//over some start .. stop
			for i := 5; i < len(row); i++ {
				// Build up the sequence
				// We do this for all sites, even after only 1 sample has every haplotype,
				// because in the backwards pass unique combinations may still result
				sampleGenos[i-5][rowIdx] = row[i][0]
			}

			// iterate overall of the samples, setting each sample's column as the "seqence"
			// then for every other sample, calculate a diff; that diff is a binary tree "word"

			// This is done in a backward pass, starting from rowIdx:
			// ---------x (rowIdx)
			// --------xx
			// -------xxx
			// ------xxxx
			// We stop when we hit max entropy (all alleles are unique)

			// Each loop of INNER gives a set of haplotypes of len(rowIdx - start + 1)
		INNER:
			for start := rowIdx; start >= 0; start-- {
				samplesToSkip := make([]bool, len(sampleGenos))

				// rootSample[start .. rowIdx] is a haplotype
				var sharedCount int
				var skipCount int
				for rootIdx, rootSample := range sampleGenos {
					if samplesToSkip[rootIdx] == true {
						skipCount++
						continue
					}

					curr = root

					hasShared = false
					hasDiff = false
					for sIdx, sample := range sampleGenos {
						// now we're here
						curr = updateTree(curr, start, rowIdx, rootSample, sample, rootIdx == sIdx)

						if sIdx == rootIdx {
							continue
						}

						if curr.val == 1 {
							if !hasShared {
								hasShared = true
							}

							// If the sample has a shared value, no need to ever
							// check it again
							samplesToSkip[sIdx] = true
							sharedCount++
							continue
						}

						// currPoint.val == 0 here
						if !hasDiff {
							hasDiff = true
						}

					}

					// log.Println(skipCount, len(alleles))
					// to avoid needing to count # of shared alleles
					if hasShared && hasDiff {
						allele := new(allele)
						allele.chr = chr
						allele.finalNode = curr

						alleles = append(alleles, allele)
						continue
					}

				}

				// reached max entropy, no need to go further for this rowIdx
				if sharedCount == 0 {
					break INNER
				}
			}
		}

		if len(alleles) > 0 {
			results <- alleles
		}

		// For this to work, uncomment some stuff above
		numUpdates := atomic.LoadUint64(&updates)
		log.Println("Done with this many tree updates: ", numUpdates, " for this many sites", nSites)
	}
}

// func reverse(s string) (result string) {
// 	for _, v := range s {
// 		result = string(v) + result
// 	}
// 	return
// }

// func makeAllele(chr string, finalNode *node) *allele {
// 	allele := new(allele)
// 	allele.chr = chr
// 	allele.finalNode = finalNode

// 	return allele
// }

func updateTree(btree *node, start, stop int, ref, alt []byte, isRef bool) *node {
	// sample index identical to the reference sample index, so obviously identical
	if !isRef {
		for i := start; i <= stop; i++ {
			if ref[i] == alt[i] {
				continue
			}

			// we found a difference, alt != ref which means a "left" split
			// if we already have a left split, returned that node
			if left := getNode(btree, 0); left != nil {
				return left
			}

			// atomic.AddUint64(&updates, 1)
			return newNode(btree, 0)
		}
	}

	// the ref and alt are identical over start .. stop
	// if we already have a left split, returned that node
	if right := getNode(btree, 1); right != nil {
		return right
	}

	// atomic.AddUint64(&updates, 1)
	return newNode(btree, 1)
}

func newNode(btree *node, val uint8) *node {
	nNode := new(node)
	nNode.val = val

	// lock for as little time as needed
	btree.Lock()
	// allow us to walk back up the tree
	// not sure if lock must be before this
	// dont think so, since we don't care what btree's
	// state is here, as long as it isn't deleted
	nNode.parent = btree

	if val == 0 {
		btree.left = nNode
	} else {
		btree.right = nNode
	}

	btree.Unlock()

	return nNode
}

// mutex unlocks in function scope;
// we can get rid of the function if it really leads to noticeable performance overhead
func getNode(btree *node, val uint8) *node {
	var n *node

	btree.RLock()
	// Manually manage locks, beceause defer takes longer
	if val == 0 {
		n = btree.left
	} else {
		n = btree.right
	}

	btree.RUnlock()

	return n
}
