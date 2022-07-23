package main

import (
	"fmt"
	"bufio"
	"compress/gzip"
	"os"
	"strings"
	"strconv"
	"regexp"
)

func printHelp() {
	fmt.Println("recombination position finder")
	fmt.Println("-h             This help")
	fmt.Println("-vcf           VCF file name")
	fmt.Println("-m             Mother id in the VCF file")
	fmt.Println("-f             Father id in the VCF file")
	fmt.Println("-w             Window size")
	fmt.Println("All other individual in the VCF file threated as offsprings")
}

func openvcf(filename string) (*bufio.Scanner) {
	file, e := os.Open(filename)
	if e != nil {
		fmt.Println("Cannot open VCF file")
		return nil
	}
	reader, e := gzip.NewReader(file)
	if e != nil {
		fmt.Println("VCF file not gzipped")
		file.Close()
		return nil
	}
	vcf := bufio.NewScanner(reader)
	return vcf
}

// Parse the genotype 0: hom ref, 1: het, 2: hom alt, 3: hard to find out what it is
func getGT(fields string) (int) {
	gt := strings.Split(fields, ":")[0]
	if len(gt) == 1 {
		return 3
	}
	if gt[0] == '0' && gt[2] == '0' {
		return 0
	}
	if (gt[0] == '0' && gt[2] == '1') || (gt[0] == '1' && gt[2] == '0') {
		return 1
	}
	if gt[0] == '1' && gt[2] == '1' {
		return 2
	}
	return 3
}

func compareClasses(prev []int, actual []int) bool {
	var same int = 0
	var comparable int = 0
	for i:=0; i < len(prev); i++ {
		if prev[i] != -1 && actual[i] != -1 {
			comparable++
		}
		if prev[i] == actual[i] && prev[i] != -1{
			same++
		}
	}
	if same > 0 && same < comparable {
		//fmt.Println("Recombination!", prev, actual, same, comparable)
		return true
	}
	return false
}

func clustering(m []int, hetcol int, homcol int) []int {
	classes := make([]int, len(m))
	for j:=0; j < len(m); j++ {
		if j == hetcol || j == homcol {
			classes[j] = -1 // -1 means do not compare the clusters
		} else {
			classes[j] = m[j] - (m[homcol] / 2)
		}
		if classes[j] == 2 {
			classes[j] = -1 // it is impossible to be a hom alt sibling when the hom parent is ref. It is a sequencing error, set it -1 (no compare)
		}
	}
	return classes

}

func processMatrix(m [][]int, positions []int, wsize int, contig string, hetcol int, homcol int) {

	if len(m) < 1 {
		return
	}
	classes := make([]int, len(m[0]))
	prevclasses := make([]int, len(m[0]))
	for i:=0; i < len(m); i++ {
		if m[i][hetcol] == 1 && (m[i][homcol] == 0 || m[i][homcol] == 2) {
			classes = clustering(m[i], hetcol, homcol)
			if i > 0 {
				if compareClasses(prevclasses, classes) {
					fmt.Println(contig, positions[i])
				}
			}
			prevclasses = classes
		}
	}
}

func processContig(m [][]int, positions []int, wsize int, contig string, mother int, father int) {
	processMatrix(m, positions, wsize, contig, mother, father)
	processMatrix(m, positions, wsize, contig, father, mother)
}

func main() {
	var vcfname string = ""
	var motherid string = ""
	var fatherid string = ""
	var windowsize int = 16

	for i:= 0; i < len(os.Args); i++ {
		if os.Args[i] == "-h" {
			printHelp()
			return
		}
		if os.Args[i] == "-vcf" {
			vcfname = os.Args[i+1]
		}
		if os.Args[i] == "-m" {
			motherid = os.Args[i+1]
		}
		if os.Args[i] == "-f" {
			fatherid = os.Args[i+1]
		}
		if os.Args[i] == "-w" {
			windowsize,_ = strconv.Atoi(os.Args[i+1])
		}
	}

	var err bool = false
	if motherid == "" {
		fmt.Println("No mother id")
		err = true
	}
	if fatherid == "" {
		fmt.Println("No father id")
		err = true
	}
	if vcfname == "" {
		fmt.Println("VCF file not specified")
		err = true
	}
	if err {
		return
	}

	vcf := openvcf(vcfname)
	if vcf == nil {
		return
	}

	var fathercol int
	var mothercol int
	var ctg string
	var prevctg string = ""
	var matrix [][]int // all windows per contig for every individual
	var row []int // sum of het genotypes in a window per every individual
	var poslist []int
	var contiglen map[string]int = make(map[string]int)

	re,rerr := regexp.Compile("##contig=<ID=([^,]+),length=([0-9]+)>")
	if rerr != nil {
		fmt.Println("Cannot compile regexp")
		return
	}
	for vcf.Scan() {
		line := vcf.Text()
		if line[1] == '#' {
			if line[0:8] == "##contig" {
				match := re.FindAllStringSubmatch(line, -1)
				//contiglen[match[0][1]],_ = strconv.ParseInt(match[0][2],10,64)
				contiglen[match[0][1]],_ = strconv.Atoi(match[0][2])
			}
			continue
		}
		cols := strings.Split(line, "\t")
		if line[0] == '#' {
			for i:= 0; i < len(cols); i++ {
				if cols[i] == motherid {
					mothercol = i
				}
				if cols[i] == fatherid {
					fathercol = i
				}
			}
			row = make([]int, len(cols) - 9)
			continue
		}
		ctg = cols[0]
		pos,err := strconv.Atoi(cols[1])
		if err != nil {
			fmt.Println("Invalid VCF file. Position is not integer:", cols[1])
			return
		}
		if len(cols[3]) != 1 || len(cols[4]) != 1 {
			// Skip indels
			continue
		}
		if ctg != prevctg {
			if prevctg != "" {
				// The contig is not the first one
				processContig(matrix, poslist, windowsize, prevctg, mothercol - 9, fathercol - 9)
			}
			// new contig, new matrix
			matrix = make([][]int,0)
			poslist = make([]int, 0)
		}
		if getGT(cols[mothercol]) != 3 && getGT(cols[fathercol]) != 3 {
		//if (getGT(cols[mothercol]) == 1 && getGT(cols[fathercol]) == 0) {//|| (getGT(cols[mothercol]) == 0 && getGT(cols[fathercol]) == 1) {
			// if the mother is het and the father is hom then
			// parse and store the offsprings genotype
			row = make([]int, len(cols) - 9)
			for i:=9; i < len(cols); i++ {
				row[i-9] = getGT(cols[i])
			}
			matrix = append(matrix, row)
			poslist = append(poslist, pos)
		}
		prevctg = ctg
	}
	processContig(matrix, poslist, windowsize, prevctg, mothercol - 9, fathercol - 9)
}
