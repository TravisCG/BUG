package main

import (
	"fmt"
	"bufio"
	"os"
	"compress/gzip"
	"strconv"
	"math/rand"
)

func printHelp() {
	fmt.Printf("fastasubsample: create a small subset of a fastafile\n")
	fmt.Printf("Usage:\n")
	fmt.Printf("fastasubsample options inputfile\n")
	fmt.Printf("-h       This help message (what a suprise)\n")
	fmt.Printf("-p       Percentage of the input file to put to the output. 0 < p < 1.0 (default: 0.1)\n")
	fmt.Printf("-g       Input file is gzipped\n")
	fmt.Printf("-o       Output file name instead of STDOUT\n")
}

func main() {
	var isgzipped bool = false
	var percentage float64 = 0.1
	var output *os.File = os.Stdout
	var err error

	if (len(os.Args) < 2) || (os.Args[1] == "-h") {
		printHelp()
		return
	}

	for i:= 1; i < len(os.Args) - 1; i++ {
		if os.Args[i] == "-p" {
			percentage, err = strconv.ParseFloat(os.Args[i+1], 64)
			if err != nil {
				percentage = 0.1
			}
		}
		if os.Args[i] == "-g" {
			isgzipped = true;
		}
		if os.Args[i] == "-o" {
			output, err = os.Create(os.Args[i+1])
			if err != nil {
				fmt.Fprintf(os.Stderr, "Cannot write output file")
				return
			}
		}
	}

	fastaname := os.Args[len(os.Args) - 1]

	file, err := os.Open(fastaname)
	if err != nil {
		file.Close()
		fmt.Println(os.Stderr, "Cannot open file:", err)
		return
	}

	var fasta *bufio.Scanner

	if isgzipped {
		reader, err := gzip.NewReader(file)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Cannot decompress the gzipped file")
			file.Close()
			return
		}
		fasta = bufio.NewScanner(reader)
	} else {
		reader := bufio.NewReader(file)
		fasta = bufio.NewScanner(reader)
	}

	var show bool = false

	for fasta.Scan() {
		line := fasta.Text()

		if line[0] == '>' {
			if rand.Float64() < percentage {
				show = true
			} else {
				show = false
			}
		}

		if show {
			fmt.Fprintf(output, "%s\n", line)
		}
	}

	file.Close()
}
