package main

import ("fmt"
	"math/rand"
	"strconv"
	"strings"
	"os"
)

type ChildHap struct {
	mother int
	father int
}

// How can we use the program
func printHelp() {
	fmt.Println("simhap: simulating haplotypes\n")
	fmt.Println("-l: length of the sequence")
	fmt.Println("-m: mother output")
	fmt.Println("-f: father output")
	fmt.Println("-vcf: vcf output")
	fmt.Println("-r: reference output name")
	fmt.Println("-c: number of childs")
	fmt.Println("-p: mutation probability (0.001)")
	fmt.Println("-s: sequencing error probability (0.001)")
}

// This writes out one sequence
func saveRef(filename string, seq string){
	f,_ := os.OpenFile(filename, os.O_CREATE | os.O_WRONLY, 0660)
	f.WriteString(">reference\n")
	f.WriteString(seq)
	f.WriteString("\n")
	f.Close()
}

// This writes out diploid sequence
func saveHap(filename string, hap1 string, hap2 string){
	f,_ := os.OpenFile(filename, os.O_CREATE | os.O_WRONLY, 0660)
	f.WriteString(">hap1\n")
	f.WriteString(hap1)
	f.WriteString("\n>hap2\n")
	f.WriteString(hap2)
	f.WriteString("\n")
	f.Close()
}

func mutate(ref string, nuc [4]string) ([2]string) {
	alt := nuc[rand.Intn(4)]
	var hap [2]string

	// finding alternative nucleotide which is not reference
	for ;alt == ref;{
		alt = nuc[rand.Intn(4)]
	}

	// which haplotype affected?
	switch rand.Intn(3) {
	case 0:
		hap[0] = alt
		hap[1] = alt
	case 1:
		hap[0] = ref
		hap[1] = alt
	case 2:
		hap[0] = alt
		hap[1] = ref
	}

	return hap
}

func vcfHeader(vcf *os.File, chrlen int64, childnum int64) {
	var i int64 = 0

        vcf.WriteString("##fileformat=VCFv4.2\n")
        vcf.WriteString("##contig=<ID=reference,length=" + strconv.FormatInt(chrlen, 10) + ">\n")
        vcf.WriteString("##INFO=<ID=RCM,Number=1,Type=Integer,Description=\"Recombination in this position, mother genome\">\n")
        vcf.WriteString("##INFO=<ID=RCF,Number=1,Type=Integer,Description=\"Recombination in this position, father genome\">\n")
        vcf.WriteString("##INFO=<ID=RMS,Number=1,Type=String,Description=\"Recombination in sample, mother genome\">\n")
        vcf.WriteString("##INFO=<ID=RFS,Number=1,Type=String,Description=\"Recombination in sample, father genome\">\n")
        vcf.WriteString("##INFO=<ID=FF,Number=1,Type=Integer,Description=\"False Flag. 0: no read error, 1: read error\">\n")
        vcf.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmother\tfather")
	for i=0; i < childnum; i++ {
		vcf.WriteString("\tchild" + strconv.FormatInt(i, 10))
	}
	vcf.WriteString("\n")
}

func vcfRecord(vcf *os.File, pos int64, ref string, mh [2]string, fh [2]string, childs []ChildHap, rmut float64, seprob float64) {
	var alts map[string]int // all the alternative nucleotides
	var alts2 []string      // all alternative nucleotides
	var index int = 1 // zero is reference
	var samples []string
	var info string
	var altstr string
	var mmutflag bool = false
	var fmutflag bool = false
	var rmsm string = ""
	var rmsf string = ""
	nuc:= [4]string{"A","T","G","C"}

	//FIXME It is a bit ugly. Maybe an array would be nicer
	alts = make(map[string]int)
	alts[ref] = 0
	alts2 = append(alts2, ref)
	if ref != mh[0] {
		alts[mh[0]] = index
		alts2 = append(alts2, mh[0])
		index++
	}
	if ref != mh[1] && mh[0] != mh[1] {
		alts[mh[1]] = index
		alts2 = append(alts2, mh[1])
		index++
	}
	a := strconv.Itoa(alts[mh[0]])
	b := strconv.Itoa(alts[mh[1]])
	samples = append(samples, a + "/" + b + ":0") //TODO there is no sequencing error in parents
	_,ok := alts[fh[0]]
	if ref != fh[0] && !ok {
		alts[fh[0]] = index
		alts2 = append(alts2, fh[0])
		index++
	}
	_,ok = alts[fh[1]]
	if ref != fh[1] && !ok {
		alts[fh[1]] = index
		alts2 = append(alts2, fh[1])
		index++
	}
	a = strconv.Itoa(alts[fh[0]])
	b = strconv.Itoa(alts[fh[1]])
	samples = append(samples, a + "/" + b + ":0") //TODO there is no sequencing error in parents

	for i:=0; i < len(childs); i++ {
		falseflag := ":0"
		// recombination
		if rand.Float64() < rmut {
			mmutflag = true
			rmsm = rmsm + ";RMS=child" + strconv.Itoa(i)
			childs[i].mother = (childs[i].mother + 1) % 2
		}
		if rand.Float64() < rmut {
			fmutflag = true
			rmsf = rmsf + ";RFS=child" + strconv.Itoa(i)
			childs[i].father = (childs[i].father + 1) % 2
		}

		a := strconv.Itoa(alts[mh[childs[i].mother]])
		b := strconv.Itoa(alts[fh[childs[i].father]])

		// sequencing error
		if rand.Float64() < seprob {
			falseflag = ":1"
			nuc1 := nuc[rand.Intn(4)]
			nuc2 := ref
			for nuc1 == ref {
				nuc1 = nuc[rand.Intn(4)]
			}
			alts[nuc1] = index
			alts2 = append(alts2, nuc1)
			index++
			if rand.Float64() < 0.5 {
				nuc2 = nuc1
			}
			a = strconv.Itoa(alts[nuc1])
			b = strconv.Itoa(alts[nuc2])
		}
		samples = append(samples, a + "/" + b + falseflag)
	}

	if mmutflag {
		info = "RCM=1" + rmsm
	} else {
		info = "RCM=0"
	}
	if fmutflag {
		info = info + ";RCF=1" + rmsf
	} else {
		info = info + ";RCF=0"
	}

	if index == 1 && mmutflag == false && fmutflag == false {
		return
	}

	altstr = strings.Join(alts2[1:], ",")
	vcf.WriteString("reference\t" + strconv.FormatInt(pos + 1, 10) + "\t.\t" + ref + "\t" + altstr + "\t100\tPASS\t" + info + "\tGT:FF")

	for i:=0; i < len(samples); i++ {
		vcf.WriteString("\t" + samples[i])
	}
	vcf.WriteString("\n")
}

func main() {
	var seqlen int64 = 2000
	var motheroutput string = "mother.out.fasta"
	var fatheroutput string = "father.out.fasta"
	var refoutput string = "reference.fasta"
	nuc:= [4]string{"A","T","G","C"}
	var reference strings.Builder
	var motherhap [2]strings.Builder
	var fatherhap [2]strings.Builder
	var mutprob float64 = 0.001
	var recmutprob float64 = 0.001
	var seprob float64 = 0.001 // Sequencing error probability
	var childnum int64 = 1
	var childs []ChildHap
	var vcfname string = "out.vcf"

	// Argument handling
	for i:=0; i < len(os.Args); i++ {
		if os.Args[i] == "-l" {
			seqlen,_ = strconv.ParseInt(os.Args[i+1], 10, 64)
		}
		if os.Args[i] == "-m" {
			motheroutput = os.Args[i+1]
		}
		if os.Args[i] == "-f" {
			fatheroutput = os.Args[i+1]
		}
		if os.Args[i] == "-r" {
			refoutput = os.Args[i+1]
		}
		if os.Args[i] == "-p" {
			mutprob,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-c" {
			childnum,_ = strconv.ParseInt(os.Args[i+1], 10, 64)
		}
		if os.Args[i] == "-vcf" {
			vcfname = os.Args[i+1]
		}
		if os.Args[i] == "-s" {
			seprob,_ = strconv.ParseFloat(os.Args[i+1], 64)
		}
		if os.Args[i] == "-h" {
			printHelp()
			return
		}
	}

	var i int64

	// define childs
	childs = make([]ChildHap, childnum)
	for i=0; i < childnum; i++ {
		childs[i].mother = rand.Intn(2)
		childs[i].father = rand.Intn(2)
	}

	vcf,_ := os.OpenFile(vcfname, os.O_CREATE | os.O_WRONLY, 0660)
	vcfHeader(vcf, seqlen, childnum)

	// building sequences
	for i=0; i < seqlen; i++ {
		var mh [2]string
		var fh [2]string
		ref := nuc[rand.Intn(4)] // this is the reference nucleotide
		mh[0] = ref //mother haplotype 1
		mh[1] = ref //mother haplotype 2
		fh[0] = ref //father haplotype 1
		fh[1] = ref //father haplotype 2

		// Time for some mutation
		if rand.Float64() < mutprob {
			mh = mutate(ref, nuc)
		}
		if rand.Float64() < mutprob {
			fh = mutate(ref, nuc)
		}

		// save VCF record
		vcfRecord(vcf, i, ref, mh, fh, childs, recmutprob, seprob)

		// store nucleotides
		motherhap[0].WriteString(mh[0])
		motherhap[1].WriteString(mh[1])
		fatherhap[0].WriteString(fh[0])
		fatherhap[1].WriteString(fh[1])
		reference.WriteString(ref)
	}

	// write sequences to file
	saveRef(refoutput, reference.String())
	saveHap(motheroutput, motherhap[0].String(), motherhap[1].String())
	saveHap(fatheroutput, fatherhap[0].String(), fatherhap[1].String())
	vcf.Close()
}
