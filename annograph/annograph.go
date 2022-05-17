package main

import ("os"
	"bufio"
	"fmt"
	"strings"
	"strconv"
	"regexp"
)

type Region struct {
	contig string;
	start  int64;
	end    int64;
	strand string;
}

func PrintHelp() {
	fmt.Println("AnnoGraph")
	fmt.Println("Using GFF annotation to scaffols fragmented genome assembly")
	fmt.Println("Usage:\n")
	fmt.Println("-h:        you already know")
	fmt.Println("-f:        genome in fasta format")
	fmt.Println("-g:        annotation in GFF3 format")
}

func ParseFasta(filename string) (map[string]int64){
	file, err := os.Open(filename)
	if err != nil {
		fmt.Println("Cannot open ",filename)
		return nil
	}
	scan := bufio.NewScanner(file)
	fasta := make(map[string]int64)
	var header string
	for scan.Scan() {
		line := scan.Text()
		if line[0] == '>' {
			// header
			header = line[1:len(line)]
			fasta[header] = 0
		} else {
			// sequence
			fasta[header] += int64(len(line))
		}
	}
	file.Close()
	return fasta
}

func ParseGFF3(filename string) (map[string][]Region){
	graph := make(map[string][]Region)

	file, err := os.Open(filename)
	if err != nil {
		fmt.Println("Cannot open", filename)
		return nil
	}
	re, err := regexp.Compile("ID=([^.]+)")
	if err != nil {
		fmt.Println("Cannot compile regexp")
		return nil
	}
	scan := bufio.NewScanner(file)
	for scan.Scan() {
		line := scan.Text()
		if line[0] == '#' {
			continue
		}
		cols := strings.Split(line,"\t")
		if cols[2] != "gene" {
			continue
		}
		contig   := cols[0]
		start,_  := strconv.ParseInt(cols[3], 10, 32)
		end,_    := strconv.ParseInt(cols[4], 10, 32)
		strand   := cols[6]
		info     := cols[8]
		match    := re.FindAllStringSubmatch(info, -1)
		genename := match[0][1]
		region   := Region{contig: contig, start: start, end: end, strand:strand,}
		_,exists := graph[genename]
		if !exists {
			graph[genename] = make([]Region, 0)
		}
		graph[genename] = append(graph[genename], region)
	}
	file.Close()
	return graph
}

func main() {
	var genomename string = ""
	var annoname   string = ""

	for i:=0; i < len(os.Args); i++ {
		if os.Args[i] == "-f" {
			genomename = os.Args[i+1]
		}
		if os.Args[i] == "-g" {
			annoname = os.Args[i+1]
		}
		if os.Args[i] == "-h" {
			PrintHelp()
			return
		}
	}

	if genomename == "" {
		fmt.Println("No genome specified")
		return
	}

	if annoname == "" {
		fmt.Println("No annotation specified")
		return
	}

	reflen := ParseFasta(genomename)
	graph  := ParseGFF3(annoname)

	for gene, regions := range(graph) {
		if len(regions) == 2 {
			from := regions[0]
			to   := regions[1]
			fromlen, e1 := reflen[from.contig]
			tolen, e2 := reflen[to.contig]
			if !e1 || !e2 {
				continue
			}
			dfrom := fromlen - from.end
			dto   := tolen - to.end
			if (dfrom < 500 && to.start < 500) || (dto < 500 && from.start < 500) {
				fmt.Println(from.contig, fromlen, from.start, from.end, gene, to.contig, tolen, to.start, to.end)
			}
		}
	}
}
