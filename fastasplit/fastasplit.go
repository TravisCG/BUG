package main

import ("fmt"
        "bufio"
	"os"
	"strconv"
)

func printhelp() {
	fmt.Printf("FastaSplit\n")
	fmt.Printf("Split fasta file in equal proportion\n")
	fmt.Printf("-h            This help\n")
	fmt.Printf("-f            Input fasta name\n")
	fmt.Printf("-r            Split into equal records\n")
	fmt.Printf("-b            Split into roughly equal size\n")
	fmt.Printf("-n            Number of parts\n")
	fmt.Printf("-p            Output prefix\n")
}

func main() {
	var records bool = true
	var parts int64 = 2
	var prefix string = "fastasplit"
	var fastaname string

	for i := 0; i < len(os.Args); i++{
		if os.Args[i] == "-h" {
			printhelp();
			return
		}
		if os.Args[i] == "-r" {
			records = true
		}
		if os.Args[i] == "-b" {
			records = false
		}
		if os.Args[i] == "-n" {
			parts,_ = strconv.ParseInt(os.Args[i+1], 10, 64)
		}
		if os.Args[i] == "-p" {
			prefix = os.Args[i+1]
		}
		if os.Args[i] == "-f" {
			fastaname = os.Args[i+1]
		}
	}

	// Read fasta file, count records, determine size of the records
	file, err := os.Open(fastaname)
	if err != nil {
		fmt.Printf("Cannot open fasta file")
		return
	}

	fastarecords := make([]int64, 0)
	var count int64 = 0

	fasta := bufio.NewScanner(file)
	for fasta.Scan() {
		line := fasta.Text()
		if line[0] == '>' {
			fastarecords = append(fastarecords, 0)
			count++ // count records
		} else {
			fastarecords[count-1] += int64(len(line)) // count size
		}
	}

	// Decide which records goes which output
	// FIXME right now is record split
	outrecs := make([][]int64, parts)
	var s int
	if count % parts == 0 {
		s = count / parts
	} else {
		s = count / parts + 1
	}

	var k int = 1
	for i:= int64(0); i < parts; i++ {
		outrecs[i] = make([]int64, s)
		for j:= 0; j < s; j++ {
			outrecs[i][j] = k
			k++
		}
	}

	// Write out results
	for i := int64(0); i < parts; i++ {
		file, err := os.Open(prefix + "." + string(i) + ".fasta", os.O_CREATE | os.O_WRONLY)
		if err != nil {
			fmt.Printf("Cannot write out results")
			return
		}
		file.Close()
	}

	// Write out records
	if records {
		fmt.Printf(prefix)
	}
	file.Close()
	return
}
