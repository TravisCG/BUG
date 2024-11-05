package main

import ("os"
        "fmt"
	"strconv"
	"compress/gzip"
	"bufio"
	"strings"
)

type Edge struct {
	nodes map[string]int64
	innodes int64
	visited bool
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
					//rcseq := revcomp(seq)

					_, exist := graph[seq[0:kmer-1]]
					if !exist {
						var actual Edge
						actual.nodes = make(map[string]int64)
						actual.nodes[seq[1:kmer]] = 0
						actual.innodes = 0
						actual.visited = false
						graph[seq[0:kmer-1]] = actual
					} else {
						_, exist = graph[seq[0:kmer-1]].nodes[seq[1:kmer]]
						if !exist {
							graph[seq[0:kmer-1]].nodes[seq[1:kmer]] = 0
						}
					}
					graph[seq[0:kmer-1]].nodes[seq[1:kmer]] += 1

					/*_, exist = graph[rcseq[0:kmer-1]]
					if !exist {
						var actual Edge
						actual.nodes = make(map[string]int64)
						actual.nodes[rcseq[1:kmer]] = 0
						actual.innodes = 0
						actual.visited = false
						graph[rcseq[0:kmer-1]] = actual
					} else {
						if !exist {
							graph[rcseq[0:kmer-1]].nodes[rcseq[1:kmer]] = 0
						}
					}
					graph[rcseq[0:kmer-1]].nodes[rcseq[1:kmer]] += 1*/
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

func isCyclic(node string, graph map[string]Edge, startnode string, depth int) (bool, string){
	edge := graph[node]
	edge.visited = true

	if depth == 0 {
		return false, ""
	}

	for k, _ := range graph[node].nodes {
		_,exist := graph[k]
		if exist {
			if k == startnode {
				return true, k
			}
			edge = graph[k]
			if edge.visited {
				return false, ""
			}
			r, cyck := isCyclic(k, graph, startnode, depth - 1)
			if r == true {
				return true, k + cyck[19:20]
			}
		}
	}
	return false, ""
}

func revertVisitedPath(node string, graph map[string]Edge, depth int) {
	edge := graph[node]
	edge.visited = false
	if depth == 0 {
		return
	}

	for k,_ := range graph[node].nodes {
		_, exist := graph[k]
		if exist {
			edge = graph[k]
			if edge.visited {
				edge.visited = false
				revertVisitedPath(k, graph, depth - 1)
			}
		}
	}
}

func graphToSNF(graph map[string]Edge) {
	for k, v := range graph {
		for node, conn:= range v.nodes {
			fmt.Println(k, conn, node)
		}
	}
}

func covHist(graph map[string]Edge) {
	for k, v := range graph {
		for _, conn := range v.nodes {
			fmt.Println(k, conn)
		}
	}
}

func edgeFilter(graph map[string]Edge, minconn int64) {
	for _, edge := range graph {
		for target, conn := range edge.nodes {
			if conn < minconn {
				delete(edge.nodes, target)
			}
		}
	}
}

func walk(kmerm1 string, graph map[string]Edge) (string){
	contig := strings.Builder{}
	contig.WriteString(kmerm1)
	k := kmerm1
	for true {
		if len(graph[k].nodes) == 1 {
			for k1,_ := range graph[k].nodes{
				k = k1
				contig.WriteString(string(k1[len(k1)-1]))
			}
		} else {
			break
		}
	}
	return contig.String()
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
	edgeFilter(graph, 3)
	calcInNodes(graph)
	for kmerm1, _ := range graph {
		r, seq := isCyclic(kmerm1, graph, kmerm1, 4000)
		if r {
			fmt.Println(seq)
		} else {
			revertVisitedPath(kmerm1, graph, 4000)
		}
	}
//	graphToSNF(graph)
//	nodecount := 1
//	for kmerm1, node := range graph {
//		if node.innodes == 0 {
//			fmt.Println(">node_",nodecount)
//			contig := walk(kmerm1, graph)
//			fmt.Println(contig)
//			nodecount += 1
//		}
//	}
}
