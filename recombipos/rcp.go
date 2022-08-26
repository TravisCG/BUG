package main

import (
	"fmt"
	"bufio"
	"compress/gzip"
	"os"
	"strings"
	"strconv"
	"regexp"
	"math"
)

type Recombi struct {
	position int
	siblings []int
}

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

func compareClasses(prev string, actual string) []int {
	var same int = 0
	var comparable int = 0
	var recombi []int
	recombi = make([]int, 0)
	for i:=0; i < len(prev); i++ {
		if prev[i] == 'x' || actual[i] == 'x' {
			continue
		}
		comparable++
		if prev[i] == actual[i]{
			same++
		}
		if prev[i] != actual[i]{
			recombi = append(recombi, i)
		}
	}
	if same > 0 && same < comparable {
		return recombi
	}
	return nil
}

func clustering(m []int, hetcol int, homcol int) string {
	classes := strings.Builder{}
	for j:=0; j < len(m); j++ {
		if j == hetcol || j == homcol {
			classes.WriteString("x")
		} else {
			clnum := m[j] - (m[homcol] / 2)
			if clnum == -1 {
				classes.WriteString("x")
			} else {
				classes.WriteString(strconv.Itoa(clnum))
			}
		}
	}
	return classes.String()
}

func normclasses(classes string, prev string) string {
	diff := 0
	anticlass := strings.Builder{}
	for j:=0; j < len(classes); j++ {
		if classes[j] != prev[j] {
			diff++
		}
		if classes[j] != 'x' {
			num,_ := strconv.Atoi(string(classes[j]))
			i := int(math.Abs(float64(num - 1)))
			anticlass.WriteString(strconv.Itoa(i))
		} else {
			anticlass.WriteString("x")
		}
	}
	if diff > len(classes) / 2 {
		return anticlass.String()
	}
	return classes
}

func getPosition(m [][]int, positions []int, startpos int, wsize int, prevclasses string, classes string, hetcol int, homcol int) int {
	var i int
	for i=0; i < wsize; i++ {
		if m[startpos + i][hetcol] == 1 && (m[startpos+i][homcol] == 0 || m[startpos + i][homcol] == 2){
			actclass := clustering(m[startpos + i], hetcol, homcol)
			if actclass == classes {
				break
			}
		}
	}
	return positions[startpos + i]
}

func processMatrix(m [][]int, positions []int, wsize int, hetcol int, homcol int) ([]Recombi, string){
	var ret []Recombi
	var classes string
	var prevclasses string
	var firstclasses string
	var first bool = true

	if len(m) == 0 {
		return ret, ""
	}

	ret = make([]Recombi, 0)

	for i:=0; i < len(m) - wsize; i++ {
		// go through the window
		classabu := make(map[string]int)
		for j:=0; j < wsize; j++ {
			if m[i+j][hetcol] == 1 && (m[i+j][homcol] == 0 || m[i+j][homcol] == 2){
				classes = clustering(m[i+j], hetcol, homcol)
				_,exists := classabu[classes]
				if !exists {
					classabu[classes] = 0
				}
				classabu[classes] += 1
			}
		}
		// get most abundant class
		max := 0
		maxclass := ""
		for k,v := range classabu {
			if v > max {
				maxclass = k
			}
		}
		if maxclass == "" {
			continue
		}
		if first {
			firstclasses = classes
			first = false
		} else {
			// clustering
			classes = normclasses(classes, prevclasses)
			r := compareClasses(prevclasses, classes)
			// recombination
			if r != nil {
				var actual Recombi
				actual.position = getPosition(m, positions, i, wsize, prevclasses, classes, hetcol, homcol)
				actual.siblings = r
				ret = append(ret, actual)
			}
		}
		prevclasses = classes
	}
	return ret, firstclasses
}

func writeBED(rp []Recombi, firstclasses string, contig string, contiglen int, sibnames []string, parent int, otherparent int){
	var prevpos []int
	var out []*os.File
	colors := []string{"255,0,0", "0,255,0"}
	prevpos = make([]int, len(sibnames))
	flipflop := make([]int, len(sibnames))
	out = make([]*os.File, len(sibnames))

	for i:=0; i < len(sibnames); i++ {
		if firstclasses[i] != 'x' && i != otherparent && i != parent{
			prevpos[i] = 0
			out[i],_ = os.OpenFile(sibnames[i] + "." + sibnames[parent] + ".bed", os.O_CREATE | os.O_APPEND | os.O_WRONLY, 0755)
			flipflop[i],_ = strconv.Atoi(string(firstclasses[i]))
		}
	}

	for i:=0; i < len(rp); i++ {
		for s:=0; s < len(rp[i].siblings); s++ {
			sib := rp[i].siblings[s]
			if sib == otherparent || sib == parent {
				continue
			}

			out[sib].WriteString(contig + "\t" + strconv.Itoa(prevpos[sib]) + "\t" + strconv.Itoa(rp[i].position - 1) + "\t0\t0\t+\t" + strconv.Itoa(prevpos[sib]) + "\t" + strconv.Itoa(rp[i].position - 1) + "\t" + colors[flipflop[sib]] + "\n")
			prevpos[sib] = rp[i].position
			flipflop[sib] = int(math.Abs(float64(flipflop[sib] - 1)))
		}
	}

	for i:=0; i < len(sibnames); i++ {
		if i == otherparent || i == parent {
			continue
		}
		if firstclasses[i] != 'x' {
			out[i].WriteString(contig + "\t" + strconv.Itoa(prevpos[i]) + "\t" + strconv.Itoa(contiglen-1) + "\t0\t0\t+\t" + strconv.Itoa(prevpos[i]) + "\t" + strconv.Itoa(contiglen-1) + "\t" + colors[flipflop[i]] + "\n")
			out[i].Close()
		}
	}
}

func processContig(m [][]int, positions []int, wsize int, contig string, contiglen int, mother int, father int, sibnames []string) {
	var recpos []Recombi
	var firstclasses string
	recpos,firstclasses = processMatrix(m, positions, wsize, mother, father)
	writeBED(recpos, firstclasses, contig, contiglen, sibnames, mother, father)
	recpos,firstclasses = processMatrix(m, positions, wsize, father, mother)
	writeBED(recpos, firstclasses, contig, contiglen, sibnames, father, mother)
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
	var sibnames []string = make([]string, 0)

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
				if i > 8 {
					sibnames = append(sibnames, cols[i])
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
				processContig(matrix, poslist, windowsize, prevctg, contiglen[prevctg], mothercol - 9, fathercol - 9, sibnames)
			}
			// new contig, new matrix
			matrix = make([][]int,0)
			poslist = make([]int, 0)
		}
		if getGT(cols[mothercol]) != 3 && getGT(cols[fathercol]) != 3 {
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
	processContig(matrix, poslist, windowsize, prevctg, contiglen[prevctg], mothercol - 9, fathercol - 9, sibnames)
}
