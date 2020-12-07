package main

import (
	"fmt"
	"regexp"
	"bufio"
	"os"
)

// How to use the program
func printHelp() {
	fmt.Printf("FastaGrep: find and extract fasta records by regular expression\n\n")
	fmt.Printf("Usage:\n")
	fmt.Printf("fastagrep -vh 'regular_expression' input_fasta\n")
	fmt.Printf("regular_expression: standard regexp\n")
	fmt.Printf("input_fasta: Fasta formatted input sequence\n")
	fmt.Printf("Results will be go to the standard output\n")
}

func main() {
	var show bool     = false
	var revmatch bool = false
	var number bool   = false

	// Argument handling
	if len(os.Args) < 3 {
		printHelp()
		return
	}

	for i:=1; i < len(os.Args) - 2; i++{
		if os.Args[i] == "-v" {
			revmatch = true
		}
		if os.Args[i] == "-h" || os.Args[i] == "-help" {
			printHelp()
			return
		}
		if os.Args[i] == "-n" {
			number = true
		}
	}

	regexpstr := os.Args[len(os.Args) - 2]
	fastaname := os.Args[len(os.Args) - 1]

	// regular expression compiling
	re, err := regexp.Compile(regexpstr)
	if err != nil {
		fmt.Println(os.Stderr, "Cannot understand regular expression:", err)
		return
	}

	// file reading operation
	file, err := os.Open(fastaname)
	if err != nil {
		file.Close()
		fmt.Println(os.Stderr, "Cannot open file:", err)
		return
	}

	fasta := bufio.NewScanner(file)
	reccount := 0
	for fasta.Scan() {
		// file processing
		line := fasta.Text()
		if line[0] == '>' {
			show = false
			reccount += 1
			match := re.Match([]byte(line))
			if (match && !revmatch) || (!match && revmatch) {
				if number {
					fmt.Println(reccount)
				} else {
					show = true
				}
			}
		}
		if show {
			fmt.Println(line)
		}
	}
	if err:= fasta.Err(); err != nil {
		fmt.Println(os.Stderr, "reading file error:", err)
	}
}
