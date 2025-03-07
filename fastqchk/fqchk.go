package main

import ("os"
        "bufio"
	"compress/gzip"
	"fmt"
	"io"
	//"strings"
)

// Check the header part
func CheckHeader(line string) bool {
	if line[0] != '@'{
		return false
	}
	return true
}

// Check the sequence
func CheckSeq(line string) (bool, int) {
	for i := 0; i < len(line); i++ {
		if line[i] != 'A' && line[i] != 'T' && line[i] != 'G' && line[i] != 'C' && line[i] != 'N' {
			fmt.Println("This is not a nucleotide:", line[i], i, line)
			return false, 0
		}
	}
	return true, len(line)
}

// Check the separator
func CheckSep(line string) bool {
	if line[0] != '+' {
		return false
	}
	return true
}

// Check the quality
func CheckQual(line string, sl int) bool {
	if sl != len(line){
		return false
	}
	for i := 0; i < len(line); i++{
		if( int(line[i]) < 33 || int(line[i]) > 126 ) {
			return false
		}
	}
	return true
}

func main() {
	var fastq *bufio.Scanner
	var good bool
	var seqlen int

	file, e := os.Open(os.Args[1])
	if e != nil {
		fmt.Println("Cannot open Fastq file")
		return
	}

	reader, e := gzip.NewReader(file)
	if e != nil {
		file.Seek(0, io.SeekStart) // gzip already read data
		fastq = bufio.NewScanner(file)
	} else {
		fastq = bufio.NewScanner(reader)
	}

	counter := 0
	for fastq.Scan() {
		line := fastq.Text()
		if counter % 4 == 0{
			good = CheckHeader(line)
			if !good {
				fmt.Println("Problem with the header in line", counter + 1)
				break
			}
		} else if counter % 4 == 1 {
			good,seqlen = CheckSeq(line)
			if !good {
				fmt.Println("Problem with the sequene in line", counter + 1)
				break
			}
		} else if counter % 4 == 2 {
			good = CheckSep(line)
			if !good {
				fmt.Println("Problem with the separator in line", counter + 1)
				break
			}
		} else if counter % 4 == 3 {
			good = CheckQual(line, seqlen)
			if !good {
				fmt.Println("Problem with the quality in line", counter + 1)
				break
			}
		}

		counter++
	}

	file.Close()
}
