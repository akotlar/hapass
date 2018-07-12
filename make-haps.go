package main

// This is Dave Cutler's idea, as a Go package. Really fucking cool compression based on the mutation/recombination events
// over the history of each sample, which will happen to be a function of 4Ne

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/binary"
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

type sequenceTree struct {
	sync.RWMutex
	// if this is a leaf/endpoint for some sequence, store the value
	val byte
	// how many times this particular sequence has recurred
	seen   uint32
	next   map[byte]*sequenceTree
	parent *sequenceTree
}

type node struct {
	sync.RWMutex
	left   *node
	right  *node
	parent *node
	val    byte
	seq    *sequenceTree
	// seq    [][]byte
}

// type allele struct {
// 	start     uint32
// 	finalNode *node
// }

type alleles struct {
	chr string
	// start uint32

	// The final combination; will walk up the tree using parentNode to actually get the full word, to avoid storing it
	// SO FUCKING COOL!
	// In the case of a single position "haplotype" (start == stop), this is the allele of 2N chromosomes (N diploid samples)
	// In the case of a multiple position haplotype (start < stop), this is whether sample N had this haplotype (1, else 0))
	// finalNode node

	// can't pass node by value; will complain about passing locks
	wordLeaves []*node
	sequences  [][]byte
	starts     []uint32
}

const diffVal = byte('0')
const sameVal = byte('1')
const tabVal = byte('\t')

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

	var mySlice = []byte{0, 1, 0, 1, 1, 0, 1, 1}
	data := binary.BigEndian.Uint64(mySlice)
	fmt.Println(data)
	// os.Exit(1)

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

	log.Println("ALT VALUE", diffVal, sameVal)
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
	diffRoot := new(node)

	// like the difference tree, but uses hashmap to track next branch
	// heavier, but still cheaper than repeating
	seqRoot := new(sequenceTree)
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
		read(config, hPath, sPath, diffRoot, seqRoot, wFunc)

		wMap.Flush()
		wPed.Flush()

		completer()
	}
}

// Future idea: pre-make the tree
// func warmTree(tree *node, nSamples int, firstNode bool) {
// 	for i := 0; i < nSamples; i++ {
// 		if tree.left == nil {
// 			tree.left = newNode(tree, diffVal)
// 		}

// 		if tree.right == nil {
// 			tree.right = newNode(tree, sameVal)
// 		}
// 	}
// }

// Generates a write function that must be called from a single thread
// for a single file
func createWriteFunc(mWriter *bufio.Writer, pWriter *bufio.Writer, baseChr int, samples []string) (accumulator func(words alleles), completer func()) {
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

	nAlleles := len(samples) * 2

	// since samples is fixed for any call to createWriteFunc,
	// we can expect (but will check)
	// that each word is of same length
	// and can therefore write this once, and over-write
	// Since we need tab characters, we double that once - 1
	word := make([]byte, nAlleles*2-1)
	for i := 1; i < len(word); i += 2 {
		word[i] = tabVal
	}

	// log.Println("word is length", word)
	// var id int
	// var idStr string
	var i int

	// var wordNode *node
	// var chr string
	// Write the haplotype to the tped file
	accumulator = func(alts alleles) {
		// Build the allele slice backwards, since we operate
		// from leaf of the binary tree
		// Every other (from index 0) entry is an allele diff (1 has, 0 no)
		// Every other (from index 1) entry is a byte('\t)
		// We track the number of values in the alt allele
		// to ensure that we overwrite every allele entry

		for idx, wordNode := range alts.wordLeaves {
			// id++

			// idStr = strconv.Itoa(id)

			mWriter.WriteString(alts.chr)
			mWriter.WriteString("\t")
			mWriter.Write(decodeSequence(alts.sequences[idx]))
			mWriter.WriteString("\t")
			// centimorgans... dummy
			mWriter.WriteString("0")
			mWriter.WriteString("\t")
			// position is basically a position; the positions aren't real
			// they're the order of the haplotypes we generate
			mWriter.WriteString(strconv.Itoa(int(alts.starts[idx])))

			i = len(word) - 1
			seen := 0
			for {
				if wordNode.parent == nil {
					break
				}

				// log.Println(i)
				word[i] = wordNode.val

				wordNode = wordNode.parent
				i -= 2
				seen++
			}

			// i should be -2 since decrement runs one extra time
			if i != -2 {
				log.Fatal("Word length must equal sample length: ", nAlleles, " got: ", i, "seen: ", seen)
			}

			mWriter.WriteString("\t")
			mWriter.Write(word)
			mWriter.WriteString("\n")
		}
	}

	return accumulator, completer
}

func read(config *config, hPath string, sPath string, btree *node, stree *sequenceTree, resultFunc func(alts alleles)) {
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
	results := make(chan alleles, 100)

	var wg sync.WaitGroup

	wg.Add(config.maxThreads)
	for i := 0; i < config.maxThreads; i++ {
		go processGenotypes(genotypeQueue, btree, stree, results, &wg)
	}

	// this could have come before processGenotypes
	go readFile(reader, genotypeQueue)

	// https://play.golang.org/p/nU2Rq6PJdO
	go func() {
		// wait for each thread to signal completion, then close the channel
		wg.Wait()
		close(results)
	}()

	// Called sequentially; critical that this is not multithreaded
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
func processGenotypes(genotypeQueue chan [][]string, root *node, stree *sequenceTree, results chan<- alleles, wg *sync.WaitGroup) {
	defer wg.Done()

	// bufferSize := 1000
	// A new item in the queue is in fact the entire cluster, block, of <= blockSize (default 1e6)
	// Genotypes contains a list of rows, of 2N + 5 length, exact length as .tped file
	// The first 4 columns are chr, id, pos, allele sequence,
	// int because more costly to convert int16 or uint16 to string; cannot use strconv.Itoa or strconf.FormatInt
	for genotypes := range genotypeQueue {
		// lSiteIdx := len(genotypes) - 1
		nSites := len(genotypes)

		log.Println("Received this many sites: ", nSites)
		// log.Println(genotypes[0])
		// 5 non-sample fields
		nSamples := len(genotypes[0]) - 5

		// these are the 2N samples, i.e all chromosomes
		// which matches nicely to the .tped format of 2 columns per sample, in .fam order
		sampleGenos := make([][]byte, nSamples)

		// Just build a vector of sites for
		for i := range sampleGenos {
			sampleGenos[i] = make([]byte, nSites)
		}

		var chr string

		var wordLeaves []*node
		var sequences [][]uint8
		var starts []uint32
		// var read int
		// var alts [][]byte
		// var alt []byte
		// var alts []uint64
		// var starts []uint32
		for rowIdx, row := range genotypes {
			chr = row[0]

			pos, err := strconv.Atoi(row[2])

			if err != nil {
				log.Fatal("couldn't read position")
			}

			// if pos < 147019380 || pos > 202133398 {
			// 	continue
			// }

			// read++

			// log.Println("Reading", pos, read)
			//Build the vector for each sample, since haplotypes consist of this vector
			//over some start .. stop
			// TODO: double check that len(row) has the expected nSamples - 5
			if len(row)-5 != nSamples {
				log.Fatal("Row is shorter than expected")
			}

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
			// var stuff []byte
		WALK:
			for start := rowIdx; start >= 0; start-- {
				samplesToSkip := make([]bool, len(sampleGenos))

				var hapCount int
				var uniqueCount int

				// haplotypeBlock[start .. rowIdx] is a haplotype
				for hIdx, haplotypeBlock := range sampleGenos {
					if samplesToSkip[hIdx] == true {
						continue
					}

					// seqLeaf, duplicateHap := getMakeSequence(stree, haplotypeBlock, start, rowIdx)

					// if duplicateHap == true {
					// 	continue
					// }

					// log.Println(seqLeaf, duplicateHap)
					// hapCount++

					// allow more than 32,768 samples
					var diffs int

					// Build from value-less root; golang pointers passed by value
					// (are not references), so ok to modify
					btree := root
					// packed := []uint8
				S_LOOP:
					for j := 0; j < len(sampleGenos); j++ {
						// if we're here, and the sample has had a haplotype
						// identical to another, it is by definition impossible
						// for it to have a haplotype that is identical to another
						// samples
						// Therefore, split left!
						if samplesToSkip[j] == true {
							diffs++

							btree = getMakeNode(btree, diffVal)

							continue
						}

						// Check the haplotype[start : rowIdx + 1] against the sample
						// If identical across entire length, split right
						// else, split left
						for i := start; i <= rowIdx; i++ {
							if haplotypeBlock[i] == sampleGenos[j][i] {
								continue
							}

							// Diff, since haplotype[i] != sampleGenoes[j][i]
							diffs++

							btree = getMakeNode(btree, diffVal)

							continue S_LOOP
						}

						// same, since haplotype[i] == sampleGenos[j][i] for all i
						samplesToSkip[j] = true
						btree = getMakeNode(btree, sameVal)
					}

					// If we're here, we've finished comparing all samples
					// to 1 haplotype
					// If that haplotype
					if diffs >= nSamples {
						log.Fatal("wtf", chr, start, rowIdx)
					}

					// skip if all same?
					// if diffs == 0 {
					// 	continue
					// }

					// Every sample unique, except reference
					if diffs == nSamples-2 {
						uniqueCount++
					}

					encoded := encodeSequence(haplotypeBlock, start, rowIdx)
					// log.Println("\n\nEncoded: ", encoded)

					// fmt.Printf("Value is %08b\n", seqPart)

					sequences = append(sequences, encoded)

					// seq := decodeSequence(encoded)

					// log.Println("Trying to reproduce of len: ", len(haplotypeBlock[start:rowIdx+1]), " of seq : ", haplotypeBlock[start:rowIdx+1])
					// log.Println("We reproduced of len: ", len(seq), " the sequence : ", seq)
					// btree.Lock()
					// if len(btree.seq) > 0 {
					// 	log.Println("HAAASSS!", btree)
					// }
					// btree.seq = append(btree.seq, haplotypeBlock[start:rowIdx+1])
					// btree.Unlock()

					wordLeaves = append(wordLeaves, btree)
					starts = append(starts, uint32(pos))
					// zstd.Compress(stuff, haplotypeBlock[start:rowIdx+1])
					// alts = append(alts, alt)
					// alts = append(alts, allele{
					// 	// seq:       getMakeSequence(btree, haplotypeBlock, start, rowIdx),
					// 	start:     uint32(start),
					// 	finalNode: btree,
					// })
				}

				// results <- alleles{
				// 	start:      uint32(start),
				// 	wordLeaves: wordLeaves,
				// 	// seqs:       alts,
				// }
				// if len(alts) > bufferSize {
				// 	results <- alleles{
				// 		chr:  chr,
				// 		alts: alts,
				// 	}

				// 	alts = make([]allele, 0, bufferSize)
				// }

				if hapCount == uniqueCount {
					log.Println("Reached max entropy", chr, start, rowIdx)
					break WALK
				}
			}
		}

		results <- alleles{
			chr:        chr,
			wordLeaves: wordLeaves,
			starts:     starts,
			sequences:  sequences,
		}

		// if len(alts) > 0 {
		// 	results <- alleles{
		// 		chr:  chr,
		// 		alts: alts,
		// 	}
		// }

		// For this to work, uncomment some stuff above
		numUpdates := atomic.LoadUint64(&updates)
		log.Println("Done with this many tree updates: ", numUpdates, " for this many sites", nSites)
	}
}

func encodeSequence(haplotypeBlock []byte, startIdx, endIdx int) []uint8 {
	var sequence []uint8

	var modulo uint8
	var y int
	// log.Println("About to : ", rowIdx, start)
	for i := endIdx; i >= startIdx; i -= 8 {
		var seqPart uint8
		y = 0
		// log.Println("Started : ", i, rowIdx, start)
		for {
			if i-y < startIdx {
				modulo = uint8(8 - y)
				sequence = append(sequence, seqPart, modulo)
				// log.Println("Broke out: ", start, rowIdx, sequence)
				// log.Printf("Value is %d (%08b) and modulo is %d\n", seqPart, seqPart, modulo)
				break
			}

			if y == 8 {
				sequence = append(sequence, seqPart)
				break
			}

			if haplotypeBlock[i-y] == '1' {
				seqPart = seqPart | 1<<uint8(y)
			} else {
				seqPart = seqPart | 0<<uint8(y)
			}

			y++
		}
	}

	return sequence
}

func decodeSequence(sequence []uint8) []byte {
	var seq []byte
	// seen := 0
	mod := sequence[len(sequence)-1]
	// first := len(sequence) - 2
	// var seqPart uint8

	// log.Printf("\nDecoding. The first %d bits of the last block %08b are omitted/masked\n", mod, sequence[len(sequence)-2])
	// log.Println("Building in reverse")
	for idx := len(sequence) - 2; idx >= 0; idx-- {
		// seen++
		var i uint8
		var mask uint8
		// seqPart = sequence[idx]
		// log.Printf("Testing %08b\n", seqPart)
		// https://stackoverflow.com/questions/4465488/how-to-go-through-each-bit-of-a-byte
		for mask = 0x80; mask != 0; mask >>= 1 {
			if idx == len(sequence)-2 && i < mod {
				i++
				continue
			}
			if sequence[idx]&mask > 0 {
				// bit is 1
				seq = append(seq, 49)
				// log.Printf("next bit is 1 of: %08b", seqPart)
			} else {
				// bit is 0
				seq = append(seq, 48)
				// log.Printf("next bit is 0 of: %08b", seqPart)
			}
		}
	}

	// log.Printf("Saw %d chunks\n\n", seen)

	// log.Println("Thinkng we reproduced of len: ", len(seq), "of seq: ", seq)

	return seq
}

func getMakeSequence(stree *sequenceTree, haplotypeBlock []byte, startIdx, endIdx int) (*sequenceTree, bool) {
	// finalNode := stree
	var diffs int
	for i := startIdx; i <= endIdx; i++ {
		if stree.next != nil && stree.next[haplotypeBlock[i]] != nil {
			stree = stree.next[haplotypeBlock[i]]
			continue
		}

		diffs++
		stree = getMakeSeqNode(stree, haplotypeBlock[i])
	}

	atomic.AddUint32(&stree.seen, 1)

	if diffs == endIdx-startIdx+1 {
		return stree, true
	}

	return stree, false
}

func getMakeSeqNode(stree *sequenceTree, val byte) (next *sequenceTree) {
	stree.RLock()

	if stree.next != nil && stree.next[val] != nil {
		next = stree.next[val]
		stree.RUnlock()

		return next
	}

	stree.RUnlock()

	next = new(sequenceTree)
	next.val = val

	stree.Lock()
	stree.next = make(map[byte]*sequenceTree)
	stree.next[val] = next
	next.parent = stree
	stree.Unlock()

	return next
}

func getMakeNode(btree *node, val byte) *node {
	// tree is read-only, except when we need to grow it
	// since most nodes will be re-visited, ok to first
	// check a node without locks, in case in our race
	// we select after creation
	n := getNodeNoLock(btree, val)

	if n != nil {
		return n
	}

	// if not, check again, this time with a read lock
	n = getNode(btree, val)

	// truly null, so lock and create a new node
	if n == nil {
		return newNode(btree, val)
	}

	return n
}

// note that byte and uint8 are equivalent
// just using byte to suggest it's a value meant to be written asis
func newNode(btree *node, val byte) *node {
	nNode := new(node)
	nNode.val = val

	// lock for as little time as needed
	btree.Lock()
	// allow us to walk back up the tree
	// not sure if lock must be before this
	// dont think so, since we don't care what btree's
	// state is here, as long as it isn't deleted
	nNode.parent = btree

	if val == diffVal {
		btree.left = nNode
	} else {
		btree.right = nNode
	}

	btree.Unlock()

	return nNode
}

// mutex unlocks in function scope;
// we can get rid of the function if it really leads to noticeable performance overhead
func getNode(btree *node, val byte) *node {
	var n *node

	btree.RLock()
	// Manually manage locks, beceause defer takes longer
	n = getNodeNoLock(btree, val)

	btree.RUnlock()

	return n
}

func getNodeNoLock(btree *node, val byte) *node {
	var n *node

	// Manually manage locks, beceause defer takes longer
	if val == diffVal {
		n = btree.left
	} else {
		n = btree.right
	}

	return n
}

// func getSeqNodeNoLock(btree *sequenceTree, val byte) *node {
// 	var n *node

// 	// Manually manage locks, beceause defer takes longer
// 	if val == diffVal {
// 		n = btree.left
// 	} else {
// 		n = btree.right
// 	}

// 	return n
// }

// func getSeqNodeNoLock(btree *node, val byte) *node {
// 	var n *sequenceTree

// 	if
// 	// Manually manage locks, beceause defer takes longer
// 	if val == diffVal {
// 		n = btree.left
// 	} else {
// 		n = btree.right
// 	}

// 	return n
// }
