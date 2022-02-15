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

// Parse the genotype 0: hom 1: het
func getGT(fields string) (int64) {
	gt := strings.Split(fields, ":")[0]
	if len(gt) < 3 {
		return 0
	}
	if gt[0] == gt[2] && (gt[0] == '0' || gt[0] == '.'){
		return 0
	}
	if (gt[0] == '0' || gt[0] == '.') && gt[2] != '0' && gt[2] != '.' {
		return 1
	}
	return 2
}

func main() {
	var vcfname string = ""
	var motherid string = ""
	var fatherid string = ""

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
	var outfiles []*os.File

	for vcf.Scan() {
		line := vcf.Text()
		if line[1] == '#' {
			continue
		}
		cols := strings.Split(line, "\t")
		if line[0] == '#' {
			outfiles = make([]*os.File, len(cols)-9)
			for i:= 9; i < len(cols); i++ {
				if cols[i] == motherid {
					mothercol = i
				}
				if cols[i] == fatherid {
					fathercol = i
				}
				outfiles[i-9],_ = os.Create(cols[i] + ".bed")
			}
			continue
		}
		ctg = cols[0]
		pos,err := strconv.ParseInt(cols[1], 10, 64)
		if err != nil {
			fmt.Println("Invalid VCF file. Position is not integer:", cols[1])
			return
		}
		if (getGT(cols[mothercol]) == 1 && getGT(cols[fathercol]) == 0) {
			// if the mother is het and the father is hom then
			// parse and store the offsprings genotype
			for i:=9; i < len(cols); i++ {
				var colour string = "0,255,0"
				if getGT(cols[i]) == 0 {
					colour = "255,0,0"
				}
				if getGT(cols[i]) == 1 {
					colour = "0,0,255"
				}
				outfiles[i-9].WriteString(ctg + "\t" + cols[1] + "\t" + strconv.FormatInt(pos + 1, 10) + "\t.\t0\t+\t" + cols[1] + "\t" + strconv.FormatInt(pos + 1, 10) + "\t" + colour + "\n")
			}
		}
	}

	for i:=0; i < len(outfiles); i++ {
		outfiles[i].Close()
	}
}
