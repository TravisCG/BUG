package main

import (
	"fmt"
	"bufio"
	"compress/gzip"
	"os"
	"strings"
	"strconv"
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
func getGT(fields string) (int64) {
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

func viterbi(Y []float64, emp [][]float64, trp [][]float64) {
	f
}

func processMatrix(m [][]int64, positions []int64, wsize int, contig string, mothercol int, fathercol int, targetcol int) {
	var eclass int

	emp := make([][]float64, 4)
	emp[0] = make([]float64, 9)
	emp[1] = make([]float64, 9)
	emp[2] = make([]float64, 9)
	emp[3] = make([]float64, 9)

	emp[0][0] = 0.44
	emp[0][1] = 0.22
	emp[0][2] = 0.0
	emp[0][3] = 0.22
	emp[0][4] = 0.11
	emp[0][5] = 0.0
	emp[0][6] = 0.0
	emp[0][7] = 0.0
	emp[0][8] = 0.0

	emp[1][0] = 0.0
	emp[1][1] = 0.0
	emp[1][2] = 0.0
	emp[1][3] = 0.22
	emp[1][4] = 0.11
	emp[1][5] = 0.0
	emp[1][6] = 0.44
	emp[1][7] = 0.22
	emp[1][8] = 0.0

	emp[2][0] = 0.0
	emp[2][1] = 0.22
	emp[2][2] = 0.44
	emp[2][3] = 0.0
	emp[2][4] = 0.11
	emp[2][5] = 0.22
	emp[2][6] = 0.0
	emp[2][7] = 0.0
	emp[2][8] = 0.0

	emp[3][0] = 0.0
	emp[3][1] = 0.0
	emp[3][2] = 0.0
	emp[3][3] = 0.0
	emp[3][4] = 0.11
	emp[3][5] = 0.22
	emp[3][6] = 0.0
	emp[3][7] = 0.22
	emp[3][8] = 0.44

	trp := make([][]float64, 4)
	trp[0] = make([]float64, 4)
	trp[1] = make([]float64, 4)
	trp[2] = make([]float64, 4)
	trp[3] = make([]float64, 4)

	trp[0][0] = 0.25
	trp[0][1] = 0.25
	trp[0][2] = 0.25
	trp[0][3] = 0.25

	trp[1][0] = 0.25
	trp[1][1] = 0.25
	trp[1][2] = 0.25
	trp[1][3] = 0.25

	trp[2][0] = 0.25
	trp[2][1] = 0.25
	trp[2][2] = 0.25
	trp[2][3] = 0.25

	trp[3][0] = 0.25
	trp[3][1] = 0.25
	trp[3][2] = 0.25
	trp[3][3] = 0.25

	for i:= 0; i < len(m) - wsize; i+=wsize {
		Y := make([]float64, wsize)
		for j:=0; j < wsize; j++ {
			eclass = int(m[i+j][fathercol] * 3 + m[i+j][mothercol])
			Y[j] = eclass
		}
		viterbi(Y, trp, emp)
	}
}

func main() {
	var vcfname string = ""
	var motherid string = ""
	var fatherid string = ""
	var windowsize int = 16
	var targetid string = ""

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
		if os.Args[i] == "-t" {
			targetid = os.Args[i+1]
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
	if targetid == "" {
		fmt.Println("Target sybling not specified")
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
	var targetcol int
	var ctg string
	var prevctg string = ""
	var matrix [][]int64 // all windows per contig for every individual
	var row []int64 // sum of het genotypes in a window per every individual
	var poslist []int64

	for vcf.Scan() {
		line := vcf.Text()
		if line[1] == '#' {
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
				if cols[i] == targetid {
					targetcol = i
				}
			}
			row = make([]int64, len(cols) - 9)
			continue
		}
		ctg = cols[0]
		pos,err := strconv.ParseInt(cols[1], 10, 64)
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
				processMatrix(matrix, poslist, windowsize, prevctg, mothercol - 9, fathercol - 9, targetcol - 9)
			}
			// new contig, new matrix
			matrix = make([][]int64,0)
			poslist = make([]int64, 0)
		}
		if (getGT(cols[mothercol]) == 1 && getGT(cols[fathercol]) == 0) {//|| (getGT(cols[mothercol]) == 0 && getGT(cols[fathercol]) == 1) {
			// if the mother is het and the father is hom then
			// parse and store the offsprings genotype
			row = make([]int64, len(cols) - 9)
			for i:=9; i < len(cols); i++ {
				row[i-9] = getGT(cols[i])
			}
			matrix = append(matrix, row)
			poslist = append(poslist, pos)
		}
		prevctg = ctg
	}
	processMatrix(matrix, poslist, windowsize, prevctg, mothercol - 9, fathercol - 9, targetcol - 9)
}
