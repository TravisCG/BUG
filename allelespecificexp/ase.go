package main

import (
	"os"
	"fmt"
	"strings"
)

func same(sampleA []float64, sampleB []float64, limit float64) bool {
	t, v := welch(sampleA, sampleB)
	p := ttestpvalue(t,v)

	if p < limit {
		return false
	}
	return true
}

func  makedecision(het float64, hethom float64, hom float64) string {
	if het > 0.05 && hethom < 0.05 && hom < 0.05 {
		return "both parents are heterozygous"
	}
	if het < 0.05 && hethom > 0.05 && hom < 0.05 {
		return "one parent is homozygous, the other one is heterozygous"
	}
	if het < 0.05 && hethom < 0.05 && hom > 0.05 {
		return "both parents are homozygous"
	}
	return "Results are not clear"
}

func main() {

	expname := ""
	var groupids []string
	var groups map[string][]float64

	for i:= 0; i < len(os.Args); i++ {
		if os.Args[i] == "-exp" {
			expname = os.Args[i+1]
		}
		if os.Args[i] == "-group" {
			groupids = strings.Split(os.Args[i+1], " ")
		}
	}

	expmatrix := expreader(expname)
	fmt.Println("gene_id\toffspring_with_parental_phenotype\tp_het\tp_het_hom\tp_hom\tdecision")
	for geneid, expvalues := range expmatrix {
		groups = make(map[string][]float64)
		for i := 0; i < len(groupids); i++ {
			groups[groupids[i]] = append(groups[groupids[i]], expvalues[i])
		}

		parpheno := 0 // number of parental phenotype
		n        := 0 // number of offsprings
		for groupid, values := range groups {
			if groupid == "F" || groupid == "M" {
				continue
			}
			if same(groups["F"], values, 0.05) || same(groups["M"], values, 0.05) {
				parpheno += 1
			}
			n += 1
		}

		hetero := binomtest(float64(n), float64(parpheno), 0.5) // both parents are heterozygous
		hethom := binomtest(float64(n), float64(parpheno), 1.0) // one parent is homozygous the other one is heterozygous
		hom    := binomtest(float64(n), float64(parpheno), 0.0) // both parents are homozygous
		deci   := makedecision(hetero, hethom, hom)
		fmt.Printf("%s\t%d\t%f\t%f\t%f\t%s\n", geneid, parpheno, hetero, hethom, hom, deci)
	}
}
