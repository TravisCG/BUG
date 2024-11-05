package main

import ("fmt"
	"os"
	"math/rand"
	"strconv"
	"strings"
)

func printHelp() {
	fmt.Println("Simulate circRNA reads")
	fmt.Println("Usage:")
	fmt.Println("-c               Number of circular RNA")
	fmt.Println("-n               Number of linear RNA")
	fmt.Println("-r               Percentage of read errors")
	fmt.Println("-e               Circular/linear RNA expression ratio")
	fmt.Println("-l               Minimal length of RNA")
	fmt.Println("-L               Maximal length of RNA")
	fmt.Println("-o               Original RNA sequence output")
	fmt.Println("-s               Seed")
}

// Return a linear sequence
func genLinear(length int64) (string) {
	var i int64
	seq := strings.Builder{}
	nuc := [4]string{"A","T","G","C"}
	for i=0; i < length; i++ {
		seq.WriteString(nuc[rand.Intn(4)])
	}
	return seq.String()
}

// Return a circular sequence. Two times the same linear sequence
// When we generate reads, some of them will overlap
func genCircular(length int64) (string) {
	s := genLinear(length)
	seq := strings.Builder{}
	seq.WriteString(s)
	seq.WriteString(s)
	return seq.String()
}

func genReads(seq string, readlen int, preaderr float64) {
	for i:= 0; i < len(seq) - readlen + 1; i++ {
		fmt.Println("@read")
		fmt.Println(seq[i:i+readlen])
		fmt.Println("+read")
		for j:= 0; j < readlen; j++ {
			fmt.Print("I")
		}
		fmt.Println()
	}
}

func main() {
	var preaderr float64
	var numcircrna int
	var numlinrna int
	var minlen float64 = 200.0
	var maxlen float64 = 400.0
	//var rnaratio float64
	var rnaoutname string = "rnaout.fasta"
	var seed int64

	if len(os.Args) == 1 {
		printHelp()
		return
	}
	for i:= 0; i < len(os.Args); i++ {
		if os.Args[i] == "-c" {
			numcircrna,_ = strconv.Atoi(os.Args[i+1])
		}
		if os.Args[i] == "-n" {
			numlinrna,_ = strconv.Atoi(os.Args[i+1])
		}
		if os.Args[i] == "-r" {
			preaderr,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-e" {
			//rnaratio,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-s" {
			seed,_ = strconv.ParseInt(os.Args[i+1], 10, 64)
		}
		if os.Args[i] == "-l" {
			minlen,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-L" {
			maxlen,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-o" {
			rnaoutname = os.Args[i+1]
		}
		if os.Args[i] == "-h" {
			printHelp()
			return
		}
	}

	outfile, err := os.Create(rnaoutname)
	if err != nil {
		fmt.Println("Cannot open output file")
		return
	}

	rand.Seed(seed)
	for i:= 0; i < numlinrna; i++ {
		seq := genLinear(int64(minlen + (maxlen - minlen) * rand.Float64()))
		outfile.WriteString(">linear" + strconv.Itoa(i) + "\n")
		outfile.WriteString(seq + "\n")
		genReads(seq, 100, preaderr)
	}

	for i:=0; i < numcircrna; i++ {
		seq := genCircular(int64(minlen + (maxlen - minlen) * rand.Float64()))
		outfile.WriteString(">circular" + strconv.Itoa(i) + "\n")
		outfile.WriteString(seq[0:len(seq)/2] + "\n")
		genReads(seq, 100, preaderr)
	}

	outfile.Close()
}
