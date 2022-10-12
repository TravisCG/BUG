package main

import ("fmt"
	"math/rand"
	"strconv"
	"strings"
	"os"
)
type ParentChr struct {
	mother int
	father int
}

func printHelp() {
	fmt.Println("-f: output reference name in fasta")
	fmt.Println("-s: seed (default 100)")
	fmt.Println("VCF will be seen in stdout")
}

func printgt(whichchr int){
	switch whichchr {
	case 0:
		fmt.Print("0/0:0")
	case 1,2:
		fmt.Print("0/1:0")
	case 3:
		fmt.Print("1/1:0")
	}
}

func childgt(inherited ParentChr, mutatedm int, mutatedf int, falsecal float64) (string){
	var frommother bool = false
	var fromfather bool = false
	var ret int = 0
	retstr := [3]string{"0/0", "0/1", "1/1"}
	var falseflag string = ":0"

	if mutatedm == 3 {
		frommother = true
	}
	if mutatedf == 3 {
		fromfather = true
	}
	if inherited.mother == 0 && mutatedm == 1 {
		frommother = true
	}
	if inherited.mother == 1 && mutatedm == 2 {
		frommother = true
	}
	if inherited.father == 0 && mutatedf == 1 {
		fromfather = true
	}
	if inherited.father == 1 && mutatedf == 2 {
		fromfather = true
	}

	if frommother && fromfather {
		ret = 2
	} else if frommother || fromfather {
		ret = 1
	}

	if rand.Float64() < falsecal {
		falseret := ret
		for ;ret == falseret; {
			falseret = rand.Intn(3)
		}
		ret = falseret
		falseflag = ":1"
	}

	return retstr[ret] + falseflag
}

func main() {
	var chrlen int = 2000000
	var mutprob float64 = 0.01
	var numchild = 8
	nuc := [4]string{ "A", "T", "G", "C" }
	var childchr []ParentChr
	var seed int64 = 100
	var recombimut float64 = 0.0001
	var falsecal float64 = 0.001
	var fasta strings.Builder
	var refname string

	for i:=0; i < len(os.Args); i++ {
		if os.Args[i] == "-f" {
			refname = os.Args[i+1]
		}
		if os.Args[i] == "-s" {
			seed,_ = strconv.ParseInt(os.Args[i+1], 10, 64)
		}
		if os.Args[i] == "-h" {
			printHelp()
			return
		}
	}

	fasta.WriteString(">onechr\n")

	rand.Seed(seed)
	fmt.Println("##fileformat=VCFv4.2")
	fmt.Println("##contig=<ID=onechr,length=" + strconv.Itoa(chrlen) + ">")
	fmt.Println("##INFO=<ID=RCM,Number=1,Type=Integer,Description=\"Recombination in this position, mother genome\">")
	fmt.Println("##INFO=<ID=RCF,Number=1,Type=Integer,Description=\"Recombination in this position, father genome\">")
	fmt.Println("##INFO=<ID=RMS,Number=1,Type=String,Description=\"Recombination in sample, mother genome\">")
	fmt.Println("##INFO=<ID=RFS,Number=1,Type=String,Description=\"Recombination in sample, father genome\">")
	fmt.Println("##INFO=<ID=FF,Number=1,Type=Integer,Description=\"False Flag. 0: no read error, 1: read error\">")
	fmt.Print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmother\tfather")
	childchr = make([]ParentChr, numchild)
	for i:= 0; i < numchild; i++ {
		fmt.Print("\tchild" + strconv.Itoa(i))
		childchr[i].mother = rand.Intn(2) // which homologue chromosome inherited from mother?
		childchr[i].father = rand.Intn(2)
	}
	fmt.Println()

	var count int = 1
	for i:= 0; i < chrlen; i++ {
		ref := nuc[rand.Intn(4)]
		if rand.Float64() < mutprob {
			ref  = nuc[rand.Intn(4)]
			alt := nuc[rand.Intn(4)]
			whichmotherchr := rand.Intn(4) // 0 neither chr, 1 first, 2 second, 3 both chr contain mutation
			whichfatherchr := rand.Intn(4) // which father haplotype contains the mutation?
			for ;alt == ref; {
				alt = nuc[rand.Intn(4)]
			}
			// recombination
			rcinfom := "RCM=0"
			rcinfof := ";RCF=0"
			for j:=0; j < numchild; j++ {
				if rand.Float64() < recombimut {
					if rcinfom == "RCM=0" {
						rcinfom = "RCM=1"
					}
					rcinfom = rcinfom + ";RMS=child" + strconv.Itoa(j)
					childchr[j].mother = (childchr[j].mother + 1) % 2
				}
				if rand.Float64() < recombimut {
					if rcinfof == ";RCF=0" {
						rcinfof = ";RCF=1"
					}
					rcinfof = rcinfof + ";RFS=child" + strconv.Itoa(j)
					childchr[j].father = (childchr[j].father + 1) % 2
				}
			}
			// put VCF record to the standard output
			fmt.Print("onechr\t")
			fmt.Print(strconv.Itoa(i+1) + "\t")
			fmt.Print("id" + strconv.Itoa(count) + "\t")
			fmt.Print(ref + "\t")
			fmt.Print(alt + "\t")
			fmt.Print("100\tPASS\t"+rcinfom+rcinfof+"\tGT:FF\t")
			printgt(whichmotherchr)
			fmt.Print("\t")
			printgt(whichfatherchr)
			for j:=0; j < numchild; j++ {
				fmt.Print("\t" + childgt(childchr[j], whichmotherchr, whichfatherchr, falsecal))
			}
			fmt.Println()
			count++
		}
		fasta.WriteString(ref)
	}
	// save reference fasta file
	f,_ := os.OpenFile(refname, os.O_CREATE | os.O_WRONLY, 0660)
	f.WriteString(fasta.String() + "\n")
	f.Close()

}
