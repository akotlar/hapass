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

type sequenceTree struct {
	// This is a read-only structure, except for next
	// which is a shared Map
	// each node gets the allele it represents
	val byte
	// how many times this particular sequence has recurred
	// Updated using atomic counter
	seen uint32
	// doubly linked
	parent *sequenceTree
	// this node maps to some other alleles as well in the chain
	// next   sync.Map //keys are bytes, values are *sequenceTree
	sync.Mutex
	bases []byte
	next  []*sequenceTree
}

type node struct {
	sync.RWMutex
	left   *node
	right  *node
	parent *node
	val    byte
	// Each diff tree has one sequence tree, with N leaves, representing N
	// alleles
	seq *sequenceTree
}

type sites struct {
	chr string
	// start uint32

	// The final combination; will walk up the tree using parentNode to actually get the full word, to avoid storing it
	// SO FUCKING COOL!
	// In the case of a single position "haplotype" (start == stop), this is the allele of 2N chromosomes (N diploid samples)
	// In the case of a multiple position haplotype (start < stop), this is whether sample N had this haplotype (1, else 0))
	// finalNode node

	// can't pass node by value; will complain about passing locks
	wordLeaves []*node
	sequences  []*sequenceTree //[]byte
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

	var outBase string
	for idx, hPath := range hapsFiles {
		dotIdx := strings.IndexByte(hPath, '.')

		if dotIdx == -1 {
			log.Printf("A file with no extension: %s; calling output by same name\n", hPath)
			outBase = hPath
		} else {
			// since dotIdx is the index of the first '.', and we want 1 byte before
			outBase = hPath[0:dotIdx]
			log.Printf("Writing files with base: %s\n", outBase)
		}

		wmFh, err := os.Create(outBase + ".tped")

		if err != nil {
			log.Fatal(err)
		}

		wfFh, err := os.Create(outBase + ".fam")

		if err != nil {
			log.Fatal(err)
		}

		wMap := bufio.NewWriter(wmFh)
		wFam := bufio.NewWriter(wfFh)

		sPath := sampleFiles[idx]

		samples, err := getSamples(sPath)

		if err != nil {
			log.Fatal(err)
		}

		if len(samples) == 0 {
			log.Fatalf("No samples found for file %s", sPath)
		}

		wFunc := createWriteFunc(wMap, wFam, idx, samples)

		fmt.Println("Reading", hPath, sPath)
		read(config, hPath, sPath, diffRoot, seqRoot, wFunc)

		wFam.Flush()
		wFam.Flush()
	}
}

// Generates a write function that must be called from a single thread
// for a single file
func createWriteFunc(mapWriter, famWriter *bufio.Writer, baseChr int, samples []string) (accumulator func(words sites)) {
	log.Println("CALLED WITH THIS MANY SAMPLES", len(samples), baseChr)

	for _, famLine := range samples {
		famWriter.WriteString(famLine)
		famWriter.WriteString("\n")
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
	var id int
	var idStr string

	// Write the haplotype to the tped file
	accumulator = func(alts sites) {
		// Build the allele slice backwards, since we operate
		// from leaf of the binary tree
		// Every other (from index 0) entry is an allele diff (1 has, 0 no)
		// Every other (from index 1) entry is a byte('\t)
		// We track the number of values in the alt allele
		// to ensure that we overwrite every allele entry

		for idx, wordNode := range alts.wordLeaves {
			id++

			// this is the last node
			diffLeaf := wordNode
			idStr = strconv.Itoa(id)

			mapWriter.WriteString(alts.chr)
			mapWriter.WriteString("\t")
			//if using the compressed version mWriter.Write(decodeSequence(alts.sequences[idx]))
			mapWriter.Write(walkSequenceTree(alts.sequences[idx]))
			mapWriter.WriteString("\t")
			// centimorgans... dummy
			mapWriter.WriteString("0")
			mapWriter.WriteString("\t")
			// starting position; not writing ending position since no space for it
			mapWriter.WriteString(strconv.Itoa(int(alts.starts[idx])))

			i := len(word) - 1
			seen := 0
			for {
				if wordNode.parent == nil {
					break
				}

				word[i] = diffLeaf.val

				wordNode = wordNode.parent
				i -= 2
				seen++
			}

			// i should be -2 since decrement runs one extra time
			if i != -2 {
				log.Fatal("Word length must equal sample length: ", nAlleles, " got: ", i, "seen: ", seen)
			}

			mapWriter.WriteString("\t")
			mapWriter.Write(word)
			mapWriter.WriteString("\n")
		}
	}

	return accumulator
}

func read(config *config, hPath string, sPath string, btree *node, stree *sequenceTree, resultFunc func(alts sites)) {
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
	results := make(chan sites, 100)

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
	var pos int
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
			log.Printf("About to make new: chr: %s, lastSentPos: %d, pos: %d, cutoff: %d, sites: %d\n",
				record[0], lastStart, pos, subChrCutoff, len(records))

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
		log.Printf("Sending remaining sites: lastSentPos: %s, pos: %d, cutoff: %d, sites: %d\n",
			lastStart, pos, subChrCutoff, len(records))
		genotypeQueue <- records
	}

	close(genotypeQueue)
}

// Process the haplotype blocks into haplotypes that are defined
// by slices of each sample's column vector over some range within that block
// then build a binary tree (shared across all threads) that diffs
// against every other sample's vector slice, tracking diffs to avoid overwork
func processGenotypes(genotypeQueue chan [][]string, root *node, stree *sequenceTree, results chan<- sites, wg *sync.WaitGroup) {
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

		// 5 non-sample columns
		rowLen := len(genotypes[0])
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
		var starts []uint32
		var sequences []*sequenceTree
		var newAlleles int
		var pos1 int
		var pos2 int

		// Haps format, has 2 alleles
		// not nec ref and alt
		alleles := make([]byte, 2)
		// var ref byte
		// var alt byte
		for end, row := range genotypes {
			chr = row[0]

			pos, err := strconv.Atoi(row[2])

			if end == 0 {
				pos1 = pos
			} else if end == nSites-1 {
				pos2 = pos
			}

			if err != nil {
				log.Fatal("couldn't read position")
			}

			alleles[0] = row[3][0]

			// In generating alleles we assume first allele is ref
			// and 2nd is alt.
			// This is of course not strictly Haps-required
			// In non-indel cases, it won't matter
			// In indel cases, things may get weird
			if len(row[4]) == len(row[3]) {
				alleles[1] = row[4][0]
			} else if len(row[4]) > len(row[3]) {
				// Insertion if len "alt" is greater than len "ref"
				// Looks to be padded, example: T TGGCTG
				// So ref is still the first byte of row[3]
				alleles[1] = '+'
			} else {
				// Deletion if len "alt" is < len "ref"
				// Have not seen this to determine padding, may be too variable
				// to handle consistently anyway
				// So for now, ref is row[3][0]
				alleles[1] = '-'
			}

			if len(row) != rowLen {
				log.Fatalf("Row length different than expect: chr:%s, pos:%d", chr, pos)
			}

			for i := 5; i < len(row); i++ {
				// Build up the sequence
				// We do this for all sites, even after only 1 sample has every haplotype,
				// because in the backwards pass unique combinations may still result
				sampleGenos[i-5][end] = row[i][0]

				// We should have no null characters
				if row[i][0] == 0 {
					log.Println("Found NULL characters", row)
					os.Exit(0)
				}
			}

			// iterate overall of the samples, setting each sample's column as the "seqence"
			// then for every other sample, calculate a diff; that diff is a binary tree "word"

			// This is done in a backward pass, starting from rowIdx:
			// Because allows us to build sampleGenos up inline
			// ---------x (rowIdx)
			// --------xx
			// -------xxx
			// ------xxxx
			// We stop when we hit max entropy (all alleles are unique)

			// Each loop gives a set of haplotypes of length end - start + 1
		WALK:
			for start := end; start >= 0; start-- {
				seenHaps := make([]bool, len(sampleGenos))
				// var haps int
				var skipCount int

				// haplotypeBlock[start .. end] is a haplotype
				for hIdx, haplotypeBlock := range sampleGenos {
					if seenHaps[hIdx] == true {
						skipCount++
						continue
					}

					seenHaps[hIdx] = true

					// Build from value-less root; golang pointers passed by value
					// (are not references), so ok to modify
					btree := root
				CHECK:
					for j := 0; j < len(sampleGenos); j++ {
						// if this is the currently examined haplotype,
						// it's of course self-identical
						if j == hIdx {
							btree = getMakeNode(btree, sameVal)

							continue CHECK
						}

						// If we've seen this haplotype before, and we've made it
						// past the seenHaps[hIdx], we know this sample's hap
						// must be different from the one being examined
						// This is the corollary to the seenHaps check above
						if seenHaps[j] == true {
							// haps++
							btree = getMakeNode(btree, diffVal)

							continue CHECK
						}

						// Check the haplotype[start : rowIdx + 1] against the sample
						// If identical across entire length, split right
						// else, split left
					DIFF:
						for i := start; i <= end; i++ {
							if sampleGenos[j][i] == 0 {
								log.Fatalf("Can't have null genotype, found at chr: %s, pos: %d, sampleIdx: %d, rowIdx: %d", chr, pos, j, i)
							}

							if haplotypeBlock[i] == sampleGenos[j][i] {
								continue DIFF
							}

							// log.Printf("Diff: %d , %d @%d %v, %v", sampleGenos[j][i], haplotypeBlock[i], i, sampleGenos[j], haplotypeBlock)

							btree = getMakeNode(btree, diffVal)

							// since we found a difference, and have made a split
							// continue to the next sample
							continue CHECK
						}

						// Past DIFF loop, so haplotypes are identical
						seenHaps[j] = true
						btree = getMakeNode(btree, sameVal)
					}

					if btree.seq == nil {
						btree.Lock()

						if btree.seq == nil {
							btree.seq = new(sequenceTree)
						}

						btree.Unlock()
					}

					seq, newAllele := getMakeSequence(btree.seq, haplotypeBlock, alleles, start, end)

					if newAllele {
						newAlleles++
					}

					sequences = append(sequences, seq)

					wordLeaves = append(wordLeaves, btree)
					starts = append(starts, uint32(pos))
				}

				haps := nSamples - skipCount
				// log.Printf("haps: %d, chr: %s, pos: %d, len: %d, startIdx: %d, endIdx %d", haps, chr, pos, end-start+1, start, end)

				if haps >= nSamples {
					log.Println("Reached max entropy", haps, nSamples, chr, start, end)
					break WALK
				}
			}
		}

		if len(wordLeaves) == 0 {
			log.Printf("Found no diff words! chr: %s, pos: %d - %d;  nSites: %d, ", chr, pos1, pos2, nSites)
		} else {
			results <- sites{
				chr:        chr,
				wordLeaves: wordLeaves,
				starts:     starts,
				sequences:  sequences,
			}
		}

		numUpdates := atomic.LoadUint64(&updates)
		log.Printf("\n\nDone with %d sites on chr %s pos:%d - %d. Found alleles: %d, diff tree updates: %d", nSites, chr, pos1, pos2, newAlleles, numUpdates)
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

func hasSeq(diffTrees []*node, desired *node) int {
	var i int
	if len(diffTrees) > 0 {
		for i = 0; i < len(diffTrees); i++ {
			if diffTrees[i] == desired {
				return i
			}
		}

		return -1
	}

	return -1
}

// Instead of using shared map, use 2 arrays (slices)...100x faster
func getMakeSequence(stree *sequenceTree, haplotypeBlock, alleles []byte, startIdx, endIdx int) (*sequenceTree, bool) {
	var diffs int
	var n *sequenceTree

	// Use slices instead of hash because hash has concurrency issues
	// and number of bases is usually <= 4
	for i := startIdx; i <= endIdx; i++ {
		// log.Println(stree)
		idx := bytes.IndexByte(stree.bases, haplotypeBlock[i])

		if idx == -1 {
			stree.Lock()

			// Check again in case things have changed
			// Cool enough, go src does this too (look at previous Map implementation)
			idx = bytes.IndexByte(stree.bases, haplotypeBlock[i])

			if idx == -1 {
				diffs++
				stree.bases = append(stree.bases, haplotypeBlock[i])

				n = new(sequenceTree)
				if haplotypeBlock[i] == '1' {
					n.val = alleles[1]
				} else if haplotypeBlock[i] == '0' {
					n.val = alleles[0]
				} else {
					log.Fatal("Haplotype must be ASCII 49 (1) or 48 (0), found: ", haplotypeBlock[i])
				}

				n.parent = stree

				stree.next = append(stree.next, n)

				stree.Unlock()
				stree = n

				continue
			}

			n = stree.next[idx]
			stree.Unlock()

			stree = n

			continue
		}

		stree = stree.next[idx]
	}

	if diffs == 0 {
		atomic.AddUint32(&stree.seen, 1)
		// log.Println("Shared:  for ", haplotypeBlock[startIdx:endIdx+1])
		return stree, true
	}

	// log.Println("New: for ", haplotypeBlock[startIdx:endIdx+1])

	return stree, false
}

// No locks needed; sequence tree is write-once per node
// and arrays can be read concurrently
func walkSequenceTree(seq *sequenceTree) []byte {
	var s []byte

	node := seq
	for {
		if node.parent == nil {
			break
		}

		node = node.parent

		s = append(s, node.val)
	}

	return s
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

// func getMakeSeqNode(stree *sequenceTree, val byte) (*{}interface, bool) {
// 	next, ok := stree.next.Load(val)

// 	if ok {
// 		return next, true
// 	}

// 	stree.RUnlock()

// 	next = new(sequenceTree)
// 	next.val = val

// 	stree.Lock()
// 	stree.next = make(map[byte]*sequenceTree)
// 	stree.next[val] = next
// 	next.parent = stree
// 	stree.Unlock()

// 	return nil, next
// }

// an inverted version; stores a *node in a global sequenceTree
// where *node is the leaf of a given diff tree/word
// func getMakeSequenceAndAddNode(stree *sequenceTree, diffTree *node, haplotypeBlock []byte, startIdx, endIdx int) (*sequenceTree, bool, bool) {
// 	stree, seen := getMakeSequence(stree, haplotypeBlock, startIdx, endIdx)

// 	i := hasSeq(stree.diffTrees, diffTree)

// 	if i > -1 {
// 		atomic.AddUint32(&stree.diffSeen[i], 1)

// 		return stree, seen, true
// 	}

// 	// check again
// 	stree.Lock()
// 	i = hasSeq(stree.diffTrees, diffTree)

// 	if i > -1 {
// 		atomic.AddUint32(&stree.diffSeen[i], 1)

// 		stree.Unlock()
// 		return stree, seen, true
// 	}

// 	if stree.diffTrees == nil {
// 		stree.diffTrees = []*node{diffTree}
// 		stree.diffSeen = []uint32{1}

// 		stree.Unlock()
// 		return stree, seen, false
// 	}

// 	stree.diffTrees = append(stree.diffTrees, diffTree)
// 	stree.diffSeen = append(stree.diffSeen, 1)

// 	stree.Unlock()

// 	return stree, seen, false
// }

// A map veresion; attractive, but slower.
// func getMakeSequenceMap(stree *sequenceTreeMap, haplotypeBlock []byte, startIdx, endIdx int) (*sequenceTreeMap, bool) {
// 	var diffs int

// 	// var next *sequenceTree
// 	// var ok bool
// 	for i := startIdx; i <= endIdx; i++ {
// 		// log.Println(stree)
// 		next, ok := stree.next.Load(haplotypeBlock[i])

// 		if ok {
// 			stree = next.(*sequenceTreeMap)
// 			continue
// 		}

// 		diffs++
// 		n := new(sequenceTreeMap)
// 		n.val = haplotypeBlock[i]
// 		n.parent = stree

// 		stree.next.Store(haplotypeBlock[i], n)

// 		stree = n
// 	}

// 	if diffs == 0 {
// 		atomic.AddUint32(&stree.seen, 1)
// 		// log.Println("Shared:  for ", haplotypeBlock[startIdx:endIdx+1])
// 		return stree, true
// 	}

// 	// log.Println("New: for ", haplotypeBlock[startIdx:endIdx+1])

// 	return stree, false
// }
