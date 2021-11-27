package main

import ("os"
        "fmt"
	"strconv"
	"compress/gzip"
	"bufio"
)

type Edge struct {
	nodes map[string]int64
	innodes int64
}

// Just print help messages
func printHelp() {
	fmt.Printf("ringring: circRNA finding program\n")
	fmt.Printf("Parameters:\n")
	fmt.Printf("-r input read\n")
	fmt.Printf("-k kmer size (default 21)\n")
}

// The sequence contains only DNA nucleotides?
func allNucs(seq string) (bool) {
	for i:= 0; i < len(seq); i++ {
		if seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G' && seq[i] != 'C' {
			return false
		}
	}
	return true
}

// Reverse complementer of a sequence
func revcomp(seq string) (string) {
	var rcseq string

	for i:= len(seq) - 1; i >= 0; i-- {
		if seq[i] == 'A' {
			rcseq += "T"
		}
		if seq[i] == 'T' {
			rcseq += "A"
		}
		if seq[i] == 'G' {
			rcseq += "C"
		}
		if seq[i] == 'C' {
			rcseq += "G"
		}
	}

	return rcseq
}

func readsToGraph(filename string, kmer int64) (map[string]Edge) {
	f, e := os.Open(filename)
	if e != nil {
		fmt.Printf("Cannot open file\n")
		return nil
	}

	r, e := gzip.NewReader(f)
	if e != nil {
		fmt.Printf("Not a gzip file\n")
		return nil
	}

	scanner := bufio.NewScanner(r)

	var nextisstuff bool = false
	var i int64
	var graph map[string]Edge
	graph = make(map[string]Edge)
	for scanner.Scan(){
		line := scanner.Text()

		if nextisstuff {
			// Processing read sequence
			if allNucs(line) {
				for i = 0; i <= int64(len(line)) - kmer; i++ {
					seq := line[i:i+kmer]
					rcseq := revcomp(seq)

					_, exist := graph[seq[0:kmer-1]]
					if !exist {
						var actual Edge
						actual.nodes = make(map[string]int64)
						actual.nodes[seq[1:kmer]] = 0
						actual.innodes = 0
						graph[seq[0:kmer-1]] = actual
					} else {
						_, exist = graph[seq[0:kmer-1]].nodes[seq[1:kmer]]
						if !exist {
							graph[seq[0:kmer-1]].nodes[seq[1:kmer]] = 0
						}
					}
					graph[seq[0:kmer-1]].nodes[seq[1:kmer]] += 1

					_, exist = graph[rcseq[0:kmer-1]]
					if !exist {
						var actual Edge
						actual.nodes = make(map[string]int64)
						actual.nodes[rcseq[1:kmer]] = 0
						actual.innodes = 0
						graph[rcseq[0:kmer-1]] = actual
					} else {
						if !exist {
							graph[rcseq[0:kmer-1]].nodes[rcseq[1:kmer]] = 0
						}
					}
					graph[rcseq[0:kmer-1]].nodes[rcseq[1:kmer]] += 1
				}
			}
		}

		if line[0] == '>' || line[0] == '@' {
			nextisstuff = true
		} else {
			nextisstuff = false
		}
	}

	f.Close()

	return graph
}

func calcInNodes(graph map[string]Edge) {
	for _,v := range graph {
		for node, _ := range v.nodes {
			edge,exist := graph[node]
			if exist {
				edge.innodes += 1
				graph[node] = edge
			}
		}
	}
}

func walk(node string, graph map[string]Edge, startnode string, contig string) (bool){

	for k, _ := range graph[node].nodes {
		if k == startnode {
			return true
		}
		contig += k[len(k)-1:]
		walk(k, graph, startnode, contig)
	}
	fmt.Println(contig)
	return false
}

func main() {

	var ifilename string
	var kmer int64 = 21

	// Parameter parsing
	for i:= 0; i < len(os.Args); i++ {
		if os.Args[i] == "-r" {
			ifilename = os.Args[i+1]
		}
		if os.Args[i] == "-h" {
			printHelp();
			return;
		}
		if os.Args[i] == "-k" {
			kmer,_ = strconv.ParseInt(os.Args[i+1], 10, 64)
		}
	}

	graph := readsToGraph(ifilename, kmer)
	calcInNodes(graph)
	for k,_ := range graph {
		//contig := k
		fmt.Println(k, graph[k].innodes, len(graph[k].nodes))
		//walk(k, graph, k, contig)
	}
}
