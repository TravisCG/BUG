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

func vcfRecord(vcf *os.File, pos int64, ref string, mh [2]string, fh [2]string, childs []ChildHap) {
	var alts map[string]int // all the alternative nucleotides
	var index int = 1

	//FIXME It is a bit ugly. Maybe an array would be nicer
	alts = make(map[string]int)
	if ref != mh[0] {
		alts[mh[0]] = index
		index++
	}
	if ref != mh[1] && mh[0] != mh[1] {
		alts[mh[1]] = index
		index++
	}
	_,ok := alts[fh[0]]
	if ref != fh[0] && !ok {
		alts[fh[0]] = index
		index++
	}
	_,ok = alts[fh[1]]
	if ref != fh[1] && !ok {
		alts[fh[1]] = index
		index++
	}

	for i:=0; i < len(childs); i++ {
		_,ok = alts[childs]
	}

	vcf.WriteString("reference\t" + strconv.FormatInt(pos + 1, 10) + "\t" + ref + "\t")
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
	var childnum int64 = 1
	var childs []ChildHap
	var vcfname string

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
	}

	var i int64

	// define childs
	childs = make([]ChildHap, childnum)
	for i=0; i < childnum; i++ {
		childs[i].mother = rand.Intn(2)
		childs[i].father = rand.Intn(2)
	}

	vcf,_ := os.OpenFile(vcfname, os.O_CREATE | os.O_WRONLY, 0660)

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
		vcfRecord(vcf, i, ref, mh, fh, childs)

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
