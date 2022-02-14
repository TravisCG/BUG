package main

import ("fmt"
	"math/rand"
	"strconv"
)

func childgt(inherited int, mutated int) (string){
	if inherited == mutated {
		return "0/1"
	}
	return "0/0"
}

func main() {
	var chrlen int = 2000000
	var mutprob float64 = 0.01
	var numchild = 8
	nuc := [4]string{ "A", "T", "G", "C" }
	var childchr []int
	var seed int64 = 100
	var recombimut float64 = 0.0001

	rand.Seed(seed)
	fmt.Println("##fileformat=VCFv4.2")
	fmt.Println("##contig=<ID=onechr,length=" + strconv.Itoa(chrlen) + ">")
	fmt.Println("##INFO=<ID=RCP,Number=1,Type=Integer,Description=\"Recombination in this position\">")
	fmt.Println("##INFO=<ID=RCS,Number=1,Type=String,Description=\"Recombination in sample\">")
	fmt.Print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmother\tfather")
	childchr = make([]int, numchild)
	for i:= 0; i < numchild; i++ {
		fmt.Print("\tchild" + strconv.Itoa(i))
		childchr[i] = rand.Intn(2) // which homologue chromosome inherited from mother?
	}
	fmt.Println()

	var count int = 1
	for i:= 0; i < chrlen; i++ {
		if rand.Float64() < mutprob {
			ref := nuc[rand.Intn(4)]
			alt := nuc[rand.Intn(4)]
			whichmotherchr := rand.Intn(2) // which mother haplotype contains the mutation?
			for ;alt == ref; {
				alt = nuc[rand.Intn(4)]
			}
			// recombination
			rcinfo := "RCP=0"
			for j:=0; j < numchild; j++ {
				if rand.Float64() < recombimut {
					if rcinfo == "RCP=0" {
						rcinfo = "RCP=1"
					}
					rcinfo = rcinfo + ";RCS=child" + strconv.Itoa(j)
					childchr[j] = (childchr[j] + 1) % 2
				}
			}
			// put VCF record to the standard output
			fmt.Print("onechr\t")
			fmt.Print(strconv.Itoa(i) + "\t")
			fmt.Print("id" + strconv.Itoa(count) + "\t")
			fmt.Print(ref + "\t")
			fmt.Print(alt + "\t")
			fmt.Print("100\tPASS\t"+rcinfo+"\tGT\t")
			fmt.Print("0/1\t")
			fmt.Print("0/0")
			for j:=0; j < numchild; j++ {
				fmt.Print("\t" + childgt(childchr[j], whichmotherchr))
			}
			fmt.Println()
			count++
		}
	}
}
