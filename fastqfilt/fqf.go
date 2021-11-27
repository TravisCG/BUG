package main

import ("os"
        "bufio"
	"compress/gzip"
	"fmt"
	"strings"
)

func printHelp() {
	fmt.Println("Fastq filter\n")
	fmt.Println("-f gzipped fastq file")
	fmt.Println("-i id file every line is a fastq read id")
	fmt.Println("-h guess what should it be")
}

func main() {
	var fastqname string
	var idfilename string
	var allow map[string]bool

	for i:= 0; i < len(os.Args); i++ {
		if os.Args[i] == "-f" {
			fastqname = os.Args[i+1]
		}
		if os.Args[i] == "-i" {
			idfilename = os.Args[i+1]
		}
		if os.Args[i] == "-h" {
			printHelp()
			return
		}
	}

	allow = make(map[string]bool)

	f, e := os.Open(idfilename)
	if e != nil {
		fmt.Println("id file missing")
		return
	}
	scan := bufio.NewScanner(f)
	for scan.Scan() {
		allow[scan.Text()] = true
	}
	f.Close()

	f, e = os.Open(fastqname)
	if e != nil {
		fmt.Println("fastq file missing")
		return
	}
	reader,e := gzip.NewReader(f)
	if e != nil {
		fmt.Println("Not a gzip file")
		return
	}
	fastq := bufio.NewScanner(reader)
	count := 1
	var output bool = false
	for fastq.Scan() {
		if count % 4 == 1 {
			fields := strings.Fields(fastq.Text())
			_, exist := allow[fields[0][1:]]
			if exist {
				output = true
			} else {
				output = false
			}
			count = 1
		}

		if output {
			fmt.Println(fastq.Text())
		}
		count += 1
	}
	f.Close()
}
