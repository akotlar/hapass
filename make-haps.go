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
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	// "unsafe"
)

type config struct {
	hPath      string
	sPath      string
	maxThreads int
}

type node struct {
	sync.RWMutex
	left   *node
	right  *node
	parent *node
	val    uint8
}

type allele struct {
	chr   string
	start uint16
	stop  uint16
	id    uint16

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
	flag.IntVar(&config.maxThreads, "threads", 8, "Num threads")
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

func createWriteFunc(mWriter *bufio.Writer, pWriter *bufio.Writer, baseChr int, samples []string) (accumulator func(word *allele), completer func()) {
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
	accumulator = func(alt *allele) {
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

	return accumulator, completer
}

func read(config *config, hPath string, sPath string, btree *node, resultFunc func(alt *allele)) {
	reader, err := getReader(hPath)

	if err != nil {
		log.Fatal(err)
	}

	// May be better to encapsulate genotypeQueue into readFile
	// store the genotypes tensor shape [ [ [], [] ] ]
	// sample: [geno1 array, geno2 array]
	// we use buffered channels
	// https://stackoverflow.com/questions/45366954/golang-channel-wierd-deadlock-related-with-make?rq=1
	genotypeQueue := make(chan [][]string)

	// The haplotypes, with sample counts, for a given sub chr
	results := make(chan *allele, 50)
	// var wg sync.WaitGroup

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
				genotypeQueue <- records
			}

			records = [][]string{}

			// fmt.Println("Made new", row[0], lastStart, pos)
		}

		// prefixes[subChr] = append(prefixes[subChr], record[:5])
		records = append(records, record)
	}

	if len(records) > 0 {
		genotypeQueue <- records
	}

	close(genotypeQueue)
}

// This function may go away; it takes an nSites by nAlleles matrix
// Note that nAlleles == 2nSamples
// Then, from index 0, figure out haplotype "words", where a word is the across-sample (i.e 2nSamples) diff
// This diff is stored in unbalanced binary tree
// The tree root is the allele chosen to diff against. Effectively it is a "1"
// The tree's first split is made by looking at the next allele (this one is from the same sample)
// If that allele is equivalent to the root allele, the tree "splits" right (i.e. gets a new Node assign to its .rNode)
// If that allele is not equivalent to the root allele, the tree "splits" left, by assigning .lNode a new Node
// And so on
// When we've done this 2N times, we store the last node (the leaf), to a new object, an allele instance of Allele
// I failed to mention: the tree is doubly linked: each node points to the parent as well as the child
// Then, when it is time to write out genotypes, we iterates over alleles <[]Allele>
// And walk up the tree to write out 2N values; for each "right split" we write a 1, for each "left split" we write a 0
// The idea here is that these trees, by being the diff against a "root" allele, will tend to recur,
// because the forces that shape the identity by state characteristic of a single locus, is effectively equivalent
// to that across nodes, i.e "haplotypes" in a phased system
// The stopping condition here is when there are 2N diffs at a single "locus" (start M to end N > M), meaning all haplotypes
// are unique, and therefore we have no power to detect association
func processGenotypes(genotypeQueue chan [][]string, root *node, results chan<- *allele, wg *sync.WaitGroup) {
	defer wg.Done()
	// Each sample is of the length of the full genotype

	// A new item in the queue is in fact the entire cluster, block, of <= blockSize (default 1e6)
	// Genotypes contains a list of rows, of 2N + 5 length, exact length as .tped file
	// The first 4 columns are chr, id, pos, allele sequence,
	// int because more costly to convert int16 or uint16 to string; cannot use strconv.Itoa or strconf.FormatInt
	for genotypes := range genotypeQueue {
		// subChr++
		// root := new(node)
		// var alleles []*allele

		// lSiteIdx := len(genotypes) - 1
		nSites := len(genotypes)

		// log.Println("launched with this many sites", nSites)
		// log.Println(genotypes[0])
		nSamples := len(genotypes[0]) - 5
		// lSampIdx := nSamples - 1

		// these are the 2N samples, i.e all chromosomes
		// which matches nicely to the .tped format of 2 columns per sample, in .fam order
		// sampleGenos are column vectors: the accumulations of genotypes across all rows
		// we typically won't build the maximum vector, since at some point only 1 sample will exist
		// per "haplotype"
		// Is sample-wise
		sampleGenos := make([][]byte, nSamples, nSamples)

		// Just build a vector of sites for
		for i := range sampleGenos {
			sampleGenos[i] = make([]byte, nSites, nSites)
		}

		// We begin constructing our haplotype words from 0 .. nSites - 1, stopping when no more shared alleles found
		// Note that the longer our allele, the more entropy it contains, and therefore when we reach a point
		// where there are no shared alleles, we have hit maximum entropy, and there is no sense in continuing, except maybe
		// if we're cold
		// lastStart := 0
		// var start int
		// var stop int

		// At each row, we will have to track # of unique haplotypes for the forward and back passes
		// when uniqueHaps == nSamples
		var chr string
		var maxEntropy int
		for rowIdx, row := range genotypes {
			chr = row[0]

			// We do
			// Single site pass:
			// ---------x
			// Backward pass   :
			// --------xx
			// -------xxx
			// ------xxxx
			// etc
			// stop backward pass when all haplotypes (unique) are private
			// Once we've stopped the backward pass, backward passes longer than the
			// length of that all-private allele are skipped/break

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

			// We grow the tree starting from the rootSample
			// Then, at each sample, we calculate the diff between it and rootSample
			// if their haplotypes are identical, we split to the right and set node.val = 1
			// if their haplotypes differ, we split left, and set node.val = 0

			// Haplotype generation loop; 2 steps:
			// 1) Since each sample column represents a haplotype
			// iterate over each sample, setting that column as the reference
			// 2) Iterate over all but that "reference" sample, computing the diff "word"
			// That word is a tree, with all of the pairwise differences
			// Notice that the word always starts with a "1"
			var seenHap []byte
			var hasShared bool
			var hasDiff bool
			var curr *node
		POINT:
			for rootIdx, rootSample := range sampleGenos {
				// biallelic only
				if len(seenHap) >= 2 {
					break
				}

				for _, alt := range seenHap {
					if alt == rootSample[rowIdx] {
						continue POINT
					}
				}

				seenHap = append(seenHap, rootSample[rowIdx])

				// we start from the top
				curr = root

				// Build the word, checking for uniqueness
				hasShared = false
				hasDiff = false
				for sIdx, sample := range sampleGenos {
					// Handle the point pass, with start == rowIdx stop == rowIdx
					// Update currPoint to keep moving down the tree
					curr = updateTree(curr, rowIdx, rowIdx, rootSample, sample, rootIdx == sIdx)

					if sIdx == rootIdx {
						continue
					}

					if curr.val == 1 {
						if !hasShared {
							hasShared = true
						}

						continue
					}

					// currPoint.val == 0 here
					if !hasDiff {
						hasDiff = true
					}
				}

				// Remove alleles where all samples homogenous and we cannot detect association
				if hasShared && hasDiff {
					results <- makeAllele(chr, rowIdx, rowIdx, curr)
				}
			}

			// Backward pass; 2 .. rowIdx, 1 .. rowIdx,  0 .. rowIdx
			// We may want to put this into a funciton, since so similar to above
			// but I save a bit of time by not doing so since backward pass
			// has slightly different needs (max entropy)
			var rootAllele string
			var seen map[string]bool
			for start := rowIdx - 1; start >= 0; start-- {
				var uniqueCount int

				if maxEntropy > 0 && rowIdx-start >= maxEntropy {
					log.Println("HIT maxEntropy")
					continue
				}

				seen = make(map[string]bool)
				for rootIdx, rootSample := range sampleGenos {
					// while kind of a slow step, saves a huge amount of computation
					// since most haplotypes will be shared
					rootAllele = string(rootSample[start : rowIdx+1])
					if seen[rootAllele] == true {
						continue
					}

					seen[rootAllele] = true

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

							continue
						}

						// currPoint.val == 0
						if !hasDiff {
							hasDiff = true
						}
					}

					// to avoid needing to count # of shared alleles
					if hasShared && hasDiff {
						results <- makeAllele(chr, rowIdx, rowIdx, curr)
						continue
					}

					if !hasShared {
						uniqueCount++
					}
				}

				if uniqueCount == len(seen) {
					log.Println(len(seen), uniqueCount)
					log.Println("ALL UNIQUE", uniqueCount, len(seen))
					maxEntropy = rowIdx - start
				}
			}
		}

		// For this to work, uncomment some stuff above
		numUpdates := atomic.LoadUint64(&updates)
		log.Println("Done with this many tree updates: ", numUpdates, " for this many sites", nSites)
	}
}

func reverse(s string) (result string) {
	for _, v := range s {
		result = string(v) + result
	}
	return
}

func makeAllele(chr string, start, stop int, finalNode *node) *allele {
	allele := new(allele)
	allele.chr = chr
	// allele.subChr = subChr
	allele.start = uint16(start)
	allele.stop = uint16(stop)
	// allele.id = id
	allele.finalNode = finalNode

	return allele
}

// can further simplify this, not sure of the overheda, by putting into function the
// code to check if present or updae
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
	// if not, make one, and assign it
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
	// log.Println(btree.right, nNode)

	btree.Unlock()

	// log.Println(btree.right, nNode)

	return nNode
}

// mutex unlocks in function scope;
// we can get rid of the function if it really leads to noticeable performance overhead
func getNode(btree *node, val uint8) *node {
	btree.RLock()
	// defer btree.RUnlock()
	// Maybe can do without an assignment here, using defer btree.RUnlock(),
	// but that may be no faster, since defer has overhead,
	// and not clear to me yet when the lock will release (before return, or after?,
	// and if after, is the returend btree.left or btree.right subject to a race condition?
	var n *node
	if val == 0 {
		n = btree.left
	} else {
		n = btree.right
	}

	btree.RUnlock()

	return n
}
