package main

import (
	"fmt"
	"bufio"
	"compress/gzip"
	"os"
	"io"
	"strings"
	"strconv"
	"regexp"
	"math"
)

type Recombi struct {
	position int
	siblings []int
}

/* Genotype information. The two haplotypes
   does not have any order, because the VCF
   also has no any */
type Genotype struct {
	h1 int
	h2 int
}

func printHelp() {
	fmt.Println("recombination position finder")
	fmt.Println("-h             This help")
	fmt.Println("-vcf           VCF file name")
	fmt.Println("-m             Mother id in the VCF file")
	fmt.Println("-f             Father id in the VCF file")
	fmt.Println("-w             Window size")
	fmt.Println("-r             Reference sequence in fasta format")
	fmt.Println("-s             Comma separated list of offsprings (sybling option)")
	fmt.Println("All other individual in the VCF file threated as offsprings (except -s is specified)")
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

func refRead(filename string) (map[string]string) {
	var header string
	var seq strings.Builder
	fastastore := make(map[string]string)

	file, e := os.Open(filename)
	if e != nil {
		fmt.Println("Cannot open reference fasta")
		file.Close()
		return nil
	}
	fasta := bufio.NewReader(file)
	for {
		line,_,err := fasta.ReadLine()
		if err != nil {
			if err == io.EOF {
				break
			}
		}
		if line[0] == '>' {
			if seq.Len() > 0{
				fastastore[header] = seq.String()
			}
			header = strings.Split(string(line[1:]), " ")[0]
			seq = strings.Builder{}
		} else {
			seq.WriteString(string(line))
		}
	}
	fastastore[header] = seq.String()
	file.Close()
	return fastastore
}

// Parse the genotype
func getGT(fields string) (Genotype) {
	var r Genotype
	rex := regexp.MustCompile("[/|]")
	gtfield := strings.Split(fields, ":")[0]
	gt := rex.Split(gtfield,-1)

	if len(gt) == 1 {
		r.h1 = -1
		r.h2 = -1
	} else {
		r.h1,_ = strconv.Atoi(gt[0])
		r.h2,_ = strconv.Atoi(gt[1])
	}

	return r
}

// Genotype is homozygous
func gtIsHom(gt Genotype) bool {
	if gt.h1 == gt.h2 {
		return true
	}
	return false
}

// The parents could have a child with the given genotype?
func isParentsChild(parent1 Genotype, parent2 Genotype, child Genotype) bool {
	if( parent1.h1 == child.h1 && parent2.h2 == child.h2){
		return true
	}
	if( parent1.h1 == child.h1 && parent2.h1 == child.h2){
		return true
	}
	if( parent1.h2 == child.h1 && parent2.h1 == child.h2){
		return true
	}
	if( parent1.h2 == child.h1 && parent2.h2 == child.h2){
		return true
	}
	if( parent2.h1 == child.h1 && parent1.h1 == child.h2){
		return true
	}
	if( parent2.h1 == child.h1 && parent1.h2 == child.h2){
		return true
	}
	if( parent2.h2 == child.h2 && parent1.h1 == child.h2){
		return true
	}
	if( parent2.h2 == child.h2 && parent1.h2 == child.h2){
		return true
	}
	return false
}

// Both alleles are reference
func gtIsRef(gt Genotype) bool {
	if gt.h1 == 0 && gt.h2 == 0{
		return true
	}
	return false
}

// Genotype is homozygous and not contains reference allele
func gtIsHomAlt(gt Genotype) bool {
	if gt.h1 > 0 && gt.h2 > 0 && gt.h1 == gt.h2 {
		return true
	}
	return false
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

func clustering(m []Genotype, hetcol int, homcol int) string {
	classes := strings.Builder{}
	for j:=0; j < len(m); j++ {
		if j == hetcol || j == homcol {
			classes.WriteString("x")
		} else {
			if(gtIsHom(m[j])){
				classes.WriteString("0")
			} else {
				classes.WriteString("1")
			}
		}
	}
	return classes.String()
}

func antiClass(classes string) string {
	anticlass := strings.Builder{}
	for j:=0; j < len(classes); j++ {
		if classes[j] != 'x' {
			num,_ := strconv.Atoi(string(classes[j]))
			i := int(math.Abs(float64(num - 1)))  //TODO 2 become 1. It can couse problem?
			anticlass.WriteString(strconv.Itoa(i))
		} else {
			anticlass.WriteString("x")
		}
	}
	return anticlass.String()
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
			i := int(math.Abs(float64(num - 1))) //TODO 2 become 1. Possible bug?
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

func getPosition(m [][]Genotype, positions []int, startpos int, wsize int, prevclasses string, classes string, hetcol int, homcol int) int {
	var i int
	for i=0; i < wsize; i++ {
		actclass := clustering(m[startpos + i], hetcol, homcol)
		anticlass := antiClass(actclass)
		if actclass == classes || anticlass == classes {
			break
		}
	}
	return positions[startpos + i]
}

/* Create a filtered matrix
   We keep variations whose can be used for recombination detection
*/
func filtMatrix(m [][]Genotype, positions []int, nucs [][]string, hetcol int, homcol int, contig string, wsize int) (fm [][]Genotype, fp []int) {
	fm = make([][]Genotype, 0)
	fp = make([]int, 0)

	for i:=0; i < len(m); i++ {
		if len(nucs[i]) > 2 || len(nucs[i][1]) > 1 || len(nucs[i][0]) > 1{
			continue
		}
		if gtIsHom(m[i][hetcol]) == false && gtIsHom(m[i][homcol]) == true {
			errorflag := false
			for j:=0; j < len(m[i]); j++ {
				if isParentsChild(m[i][hetcol], m[i][homcol], m[i][j]) == false {
					fmt.Println("Sequencing error:", contig, positions[i], j)
					errorflag = true
				}
			}
			if errorflag {
				continue
			}
			fm = append(fm, m[i])
			fp = append(fp, positions[i])
			if len(fp) > 1 && fp[len(fp) - 1] - fp[len(fp) - 2] > wsize {
				fmt.Println("No enough variation in the given window:", contig, fp[len(fp)-1], fp[len(fp)-2])
			}
		}
	}
	return fm, fp
}

// Calculate Shannon entropy
func ShannonEnt(classabu map[string]int) float64 {
	var sum int = 0
	var shannon float64 = 0.0
	var p float64

	for _,v := range classabu {
		sum += v
	}

	for _,v := range classabu {
		p = float64(v) / float64(sum)
		shannon += p * math.Log(p) * -1.0
	}

	return shannon
}

func getMaxClass(m [][]Genotype, startpos int, wsize int, hetcol int, homcol int) (string, float64) {
	var classes string
	classabu := make(map[string]int)
	for j:=0; j < wsize; j++ {
		classes = clustering(m[startpos + j], hetcol, homcol)
		_,exists := classabu[classes]
		if !exists {
			classabu[classes] = 0
		}
		classabu[classes]++
	}
	ent := ShannonEnt(classabu)
	max := 0
	maxclass := ""
	for k,v := range classabu {
		if v > max {
			max = v
			maxclass = k
		}
	}
	return maxclass, ent
}

func processMatrix(m [][]Genotype, positions []int, wsize int, hetcol int, homcol int, contig string) ([]Recombi, string){
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
		maxclass, entropy := getMaxClass(m, i, wsize, hetcol, homcol)
		if maxclass == "" || entropy > 2.0 {
			continue
		}
		if first {
			firstclasses = maxclass
			first = false
			prevclasses = maxclass
		} else {
			// clustering
			classes = normclasses(maxclass, prevclasses)
			r := compareClasses(prevclasses, classes)
			// recombination
			if r != nil {
				var actual Recombi
				actual.position = getPosition(m, positions, i, wsize, prevclasses, classes, hetcol, homcol)
				actual.siblings = r
			
				if len(ret) == 0 || actual.position > ret[len(ret)-1].position{
					ret = append(ret, actual)
				}
			}
			prevclasses = classes
		}
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
	// if there is no recombination at all, we do not put any output
	if len(firstclasses) == 0 {
		return
	}

	for i:=0; i < len(sibnames); i++ {
		if firstclasses[i] != 'x' && i != otherparent && i != parent{
			prevpos[i] = 0
			out[i],_ = os.OpenFile(sibnames[i] + "." + sibnames[parent] + ".bed", os.O_CREATE | os.O_APPEND | os.O_WRONLY, 0644)
			flipflop[i],_ = strconv.Atoi(string(firstclasses[i]))
			if flipflop[i] > 1 {
				flipflop[i] = 1
			}
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
			out[i].WriteString(contig + "\t")
			out[i].WriteString(strconv.Itoa(prevpos[i]) + "\t")
			out[i].WriteString(strconv.Itoa(contiglen-1) + "\t0\t0\t+\t")
			out[i].WriteString(strconv.Itoa(prevpos[i]) + "\t")
			out[i].WriteString(strconv.Itoa(contiglen-1) + "\t")
			out[i].WriteString(colors[flipflop[i]] + "\n")
			out[i].Close()
		}
	}
}

func hapDecision(flipflop []int, nucs []string, genotype []Genotype, parentcol int) (string, string) {
	hap1 := make(map[Genotype]int, 0)
	hap2 := make(map[Genotype]int, 0)

	// collect all the genotypes per groups
	for i:=0; i < len(flipflop); i++ {
		if flipflop[i] == 2 {
			// Parent column maybe. I will use this information, but not right now
		} else {
			if flipflop[i] == 0 {
				// first haplotype
				_,e := hap1[genotype[i]]
				if !e {
					hap1[genotype[i]] = 0
				}
				hap1[genotype[i]]++
			} else {
				// second haplotype
				_,e := hap2[genotype[i]]
				if !e {
					hap2[genotype[i]] = 0
				}
				hap2[genotype[i]]++
			}
		}
	}

	// voting for the most abundant genotype in the group
	max := 0
	var hmax1 Genotype //:= 0
	hetnum1 := 0
	for k,v:= range(hap1){
		if v > max{
			max = v
			hmax1 = k
		}
		if gtIsHom(k) == false {
			hetnum1 += 1
		}
	}
	max = 0
	var hmax2 Genotype //:= 0
	hetnum2 := 0
	for k,v := range(hap2){
		if v > max {
			max = v
			hmax2 = k
		}
		if gtIsHom(k) == false {
			hetnum2 += 1
		}
	}

	// decision
	var hap1str string
	var hap2str string

	if gtIsRef(hmax1) && gtIsRef(hmax2) {
		if hetnum1 > hetnum2 {
			hap1str = nucs[genotype[parentcol].h1]
			hap2str = nucs[0]
		} else {
			hap1str = nucs[0]
			hap2str = nucs[genotype[parentcol].h2]
		}
	} else if gtIsHom(hmax1) == false && gtIsRef(hmax2) {
		if genotype[parentcol].h1 != 0 {
			hap1str = nucs[genotype[parentcol].h1]
		} else {
			hap1str = nucs[genotype[parentcol].h2]
		}
		hap2str = nucs[0]
	} else if gtIsHomAlt(hmax1) && gtIsRef(hmax2) {
		hap1str = nucs[hmax1.h1]
		hap2str = nucs[0]
	} else if gtIsRef(hmax1) && gtIsHom(hmax2) == false {
		hap1str = nucs[0]
		if genotype[parentcol].h1 != 0 {
			hap2str = nucs[genotype[parentcol].h1]
		} else {
			hap2str = nucs[genotype[parentcol].h2]
		}
	} else if gtIsHom(hmax1) == false && gtIsHom(hmax2) == false {
		hap1str = nucs[0]
		hap2str = nucs[hmax2.h1] //FIXME there is no way to find out the correct haplotype, right now it is just a dummy thing
	} else if gtIsHomAlt(hmax1) && gtIsHom(hmax2) == false {
		hap1str = nucs[hmax1.h1]
		hap2str = nucs[0]
	} else if gtIsRef(hmax1) && gtIsHomAlt(hmax2) {
		hap1str = nucs[0]
		hap2str = nucs[hmax2.h1]
	} else if gtIsHom(hmax1) == false && gtIsHomAlt(hmax2) {
		hap1str = nucs[0]
		hap2str = nucs[hmax2.h1]
	} else if gtIsHomAlt(hmax1) && gtIsHomAlt(hmax2) {
		hap1str = nucs[hmax1.h1]
		hap2str = nucs[hmax2.h1]
	}
	return hap1str, hap2str
}

func writeHAP(m [][]Genotype, positions []int, firstclasses string, nucs [][]string, recpos []Recombi, parentcol int, outfilename string, contig string, fastarecord string) {
	hap1 := strings.Builder{}
	hap2 := strings.Builder{}
	flipflop := make([]int, len(m[0]))
	var hapstart int = 0
	var rpindex int = 0

	if len(firstclasses) == 0 {
		return
	}
	for i:=0; i < len(m[0]); i++ {
		if firstclasses[i] != 'x' {
			flipflop[i],_ = strconv.Atoi(string(firstclasses[i]))
			if flipflop[i] > 1 {
				flipflop[i] = 1
			}
		} else {
			flipflop[i] = 2
		}
	}
	for i:=0; i < len(m); i++ {
		pos := positions[i] - 1
		// change the siblings groups if recombination occured
		if rpindex < len(recpos) && recpos[rpindex].position < pos {
			for s:=0; s < len(recpos[rpindex].siblings); s++ {
				sib := recpos[rpindex].siblings[s]
				flipflop[sib] = int(math.Abs(float64(flipflop[sib] - 1)))
			}
			rpindex++
		}
		// appending non-variable sequences from reference
		if pos > hapstart || pos == hapstart {
			hap1.WriteString(fastarecord[hapstart:pos])
			hap2.WriteString(fastarecord[hapstart:pos])
		} else {
			continue
		}
		// decide which variation goes to which haplotype
		if gtIsHom(m[i][parentcol]) == false {
			// heterozygous parent. THe decision can be difficult
			n1, n2 := hapDecision(flipflop, nucs[i], m[i], parentcol)
			hap1.WriteString(n1)
			hap2.WriteString(n2)
			pos = pos + len(nucs[i][0])
		} else {
			hap1.WriteString(nucs[i][ m[i][parentcol].h1 ])
			hap2.WriteString(nucs[i][ m[i][parentcol].h2 ])
			pos = pos + len(nucs[i][0])
		}
		hapstart = pos
	}
	hap1.WriteString(fastarecord[hapstart:])
	hap2.WriteString(fastarecord[hapstart:])

	out,_ := os.OpenFile(outfilename, os.O_APPEND | os.O_CREATE | os.O_WRONLY, 0644)
	out.WriteString(">" + contig + "_hap1\n")
	out.WriteString(hap1.String() + "\n")
	out.WriteString(">" + contig + "_hap2\n")
	out.WriteString(hap2.String() + "\n")
	out.Close()
}

func processContig(m [][]Genotype, positions []int, wsize int, contig string, contiglen int, mother int, father int, sibnames []string, nucs [][]string, fastarecord string) {
	var recpos []Recombi
	var firstclasses string
	var filtm [][]Genotype
	var filtp []int

	filtm,filtp = filtMatrix(m, positions, nucs, mother, father, contig, wsize)
	recpos,firstclasses = processMatrix(filtm, filtp, wsize, mother, father, contig)
	writeBED(recpos, firstclasses, contig, contiglen, sibnames, mother, father)
	writeHAP(m, positions, firstclasses, nucs, recpos, mother, "motherhaplotype.fasta", contig, fastarecord)

	filtm,filtp = filtMatrix(m, positions, nucs, father, mother, contig, wsize)
	recpos,firstclasses = processMatrix(filtm, filtp, wsize, father, mother, contig)
	writeBED(recpos, firstclasses, contig, contiglen, sibnames, father, mother)
	writeHAP(m, positions, firstclasses, nucs, recpos, father, "fatherhaplotype.fasta", contig, fastarecord)
}

func main() {
	var vcfname string = ""
	var refname string = ""
	var motherid string = ""
	var fatherid string = ""
	var windowsize int = 16
	var allowedSyblings []string

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
		if os.Args[i] == "-r" {
			refname = os.Args[i+1]
		}
		if os.Args[i] == "-s" {
			allowedSyblings = strings.Split(os.Args[i+1], ",")
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
	if refname == "" {
		fmt.Println("Reference fasta file missing")
		err = true
	}
	if err {
		return
	}

	fasta := refRead(refname)

	vcf := openvcf(vcfname)
	if vcf == nil {
		return
	}

	var fathercol int
	var mothercol int
	var remappedf int
	var remappedm int
	var ctg string
	var prevctg string = ""
	var matrix [][]Genotype // all windows per contig for every individual
	var row []Genotype // sum of het genotypes in a window per every individual
	var poslist []int
	var nucs [][]string
	var contiglen map[string]int = make(map[string]int)
	var samplenames []string = make([]string, 0)
	var allowedindex []int = make([]int, 0)

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
					remappedm = i - 9
				}
				if cols[i] == fatherid {
					fathercol = i
					remappedf = i - 9
				}
				if i > 8 {
					if len(allowedSyblings) > 0 {
						if i == fathercol || i == mothercol {
							allowedindex = append(allowedindex, i)
							samplenames = append(samplenames, cols[i])
							if i == fathercol {
								remappedf = len(samplenames) - 1
							} else {
								remappedm = len(samplenames) - 1
							}
						}
						for _, sa := range(allowedSyblings){
							if sa == cols[i] {
								allowedindex = append(allowedindex, i)
								samplenames = append(samplenames, sa)
							}
						}
					} else {
						samplenames = append(samplenames, cols[i])
						allowedindex = append(allowedindex, i)
					}
				}
			}
			row = make([]Genotype, len(samplenames))
			continue
		}
		ctg = cols[0]
		pos,err := strconv.Atoi(cols[1])
		if err != nil {
			fmt.Println("Invalid VCF file. Position is not integer:", cols[1])
			return
		}
		if ctg != prevctg {
			if prevctg != "" {
				// The contig is not the first one
				seq, exists := fasta[prevctg]
				if exists {
					processContig(matrix, poslist, windowsize, prevctg, contiglen[prevctg], remappedm, remappedf, samplenames, nucs, seq)
				} else {
					fmt.Println(prevctg,"not found in the reference file")
				}
			}
			// new contig, new matrix
			matrix = make([][]Genotype,0)
			poslist = make([]int, 0)
			nucs = make([][]string, 0) // store alternative nucleotides
		}
		// parse and store the offsprings genotype
		row = make([]Genotype, len(samplenames))
		for i:=9; i < len(cols); i++ {
			for ii,index := range(allowedindex){
				if index == i {
					row[ii] = getGT(cols[i])
				}
			}
		}
		if cols[4] != "*" {
			matrix = append(matrix, row)
			poslist = append(poslist, pos)
			actnucs := make([]string, 0)
			actnucs = append(actnucs, cols[3])
			actnucs = append(actnucs, strings.Split(cols[4], ",")...)
			nucs = append(nucs, actnucs)
		}

		prevctg = ctg
	}
	processContig(matrix, poslist, windowsize, prevctg, contiglen[prevctg], remappedm, remappedf, samplenames, nucs, fasta[prevctg])
}
