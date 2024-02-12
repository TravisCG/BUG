package main

import ("fmt"
	"math/rand"
	"strconv"
	"strings"
	"os"
)

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

func mutate(ref string, nuc [4]string) (string, string) {
	alt := nuc[rand.Intn(4)]
	var hap1 string
	var hap2 string

	// finding alternative nucleotide which is not reference
	for ;alt == ref;{
		alt = nuc[rand.Intn(4)]
	}

	// which haplotype affected?
	switch rand.Intn(3) {
	case 0:
		hap1 = alt
		hap2 = alt
	case 1:
		hap1 = ref
		hap2 = alt
	case 2:
		hap1 = alt
		hap2 = ref
	}

	return hap1, hap2
}

func main() {
	var seqlen int64 = 2000
	var motheroutput string = "mother.out.fasta"
	var fatheroutput string = "father.out.fasta"
	var refoutput string = "reference.fasta"
	nuc:= [4]string{"A","T","G","C"}
	var reference strings.Builder
	var motherhap1 strings.Builder
	var motherhap2 strings.Builder
	var fatherhap1 strings.Builder
	var fatherhap2 strings.Builder
	var mutprob float64 = 0.001

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
	}

	// building sequences
	var i int64
	for i=0; i < seqlen; i++ {
		ref := nuc[rand.Intn(4)] // this is the reference nucleotide
		mh1 := ref //mother haplotype 1
		mh2 := ref //mother haplotype 2
		fh1 := ref //father haplotype 1
		fh2 := ref //father haplotype 2

		// Time for some mutation
		if rand.Float64() < mutprob {
			mh1, mh2 = mutate(ref, nuc)
		}
		if rand.Float64() < mutprob {
			fh1, fh2 = mutate(ref, nuc)
		}

		// store nucleotides
		motherhap1.WriteString(mh1)
		motherhap2.WriteString(mh2)
		fatherhap1.WriteString(fh1)
		fatherhap2.WriteString(fh2)
		reference.WriteString(ref)
	}

	// write sequences to file
	saveRef(refoutput, reference.String())
	saveHap(motheroutput, motherhap1.String(), motherhap2.String())
	saveHap(fatheroutput, fatherhap1.String(), fatherhap2.String())

}
