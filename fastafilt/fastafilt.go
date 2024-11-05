package main

import ("fmt"
	"bufio"
	"compress/gzip"
	"os"
	"strconv"
	"strings"
)

func printHelp() {
	fmt.Println("fastafilt:")
	fmt.Println("-fasta: fasta file name")
	fmt.Println("-lmin: length filter minimum")
	fmt.Println("-lmax: length filter maximum")
	fmt.Println("-gmin: GC minimum")
	fmt.Println("-gmax: GC maximum")
	fmt.Println("-nmin: Minumum N content")
	fmt.Println("-nmax: Maximum N content")
}

func decide(header string, seq strings.Builder, lmin int64, lmax int64, gmin float64, gmax float64, nmin float64, nmax float64) {
	var ncount int64 = 0
	var gccount int64 = 0

	s := seq.String()
	if int64(len(s)) < lmin {
		return
	}
	if int64(len(s)) > lmax && lmax != -1{
		return
	}
	for i:=0; i < len(s); i++ {
		if s[i] == 'N' || s[i] == 'n' {
			ncount++
		}
		if s[i] == 'G' || s[i] == 'C' || s[i] == 'g' || s[i] == 'c' {
			gccount++
		}
	}

	fmt.Println(header)
	fmt.Println(s)
}

func main() {
	var lmin int64 = 0
	var lmax int64 = -1
	var gmin float64
	var gmax float64
	var nmin float64
	var nmax float64
	var fastaname string
	var fasta *bufio.Scanner
	var file *os.File
	var err error
	var gzipped bool = false
	var seq strings.Builder
	var header string

	for i:=1; i < len(os.Args); i++ {
		if os.Args[i] == "-h" {
			printHelp()
			return
		}
		if os.Args[i] == "-lmin" {
			lmin,_ = strconv.ParseInt(os.Args[i+1], 10, 64)
		}
		if os.Args[i] == "-lmax" {
			lmax,_ = strconv.ParseInt(os.Args[i+1], 10, 64)
		}
		if os.Args[i] == "-gmin" {
			gmin,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-gmax" {
			gmax,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-nmin" {
			nmin,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-nmax" {
			nmax,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-fasta" {
			fastaname = os.Args[i+1]
		}
		if os.Args[i] == "-g" {
			gzipped = true
		}
	}

	file, err = os.Open(fastaname)
	if err != nil {
		fmt.Println(os.Stderr, "Cannot open file")
		return
	}

	if gzipped {
		reader, err := gzip.NewReader(file)
		if err != nil {
			fmt.Println(os.Stderr, "Cannot decompress file")
			file.Close()
			return
		}
		fasta = bufio.NewScanner(reader)
	} else {
		reader := bufio.NewReader(file)
		fasta = bufio.NewScanner(reader)
	}


	for fasta.Scan() {
		line := fasta.Text()
		if line[0] == '>' {
			if seq.Len() > 0 {
				decide(header, seq, lmin, lmax, gmin, gmax, nmin, nmax)
			}
			header = line
			seq = strings.Builder{}
		} else {
			seq.WriteString(line)
		}
	}
	decide(header, seq, lmin, lmax, gmin, gmax, nmin, nmax)

	file.Close()
}
