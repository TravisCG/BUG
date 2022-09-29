package main

import (
	"os"
	"github.com/biogo/hts/bam"
	"io"
	"fmt"
	"github.com/biogo/hts/sam"
	"math"
)

type Window struct {
	Clipped []int64;
	Total   []int64;
}

func bincoef(k int64, n int64) float64 {
	var s float64 = 1.0
	var i int64
	if n-k < k {
		k = n - k
	}
	for i=0; i < k; i++ {
		s *= float64(n-i) / float64(i+1)
	}
	return s
}

func binomTest(k int64, n int64, p float64) float64 {
	var s float64 = 0.0
	var i int64
	for i=k; i <= n; i++ {
		s += bincoef(i,n) * math.Pow(p,float64(i)) * math.Pow(1.0 - p, float64(n - i))
	}
	return s
}

func findProblems(storage map[string]Window, mincoverage int64, limit float64) {
	for refName,window := range storage {
		l := len(window.Clipped)
		for i:=0; i < l; i++ {
			if window.Total[i] >= mincoverage {
				p := binomTest(window.Clipped[i], window.Total[i], 0.5)
				if p < limit {
					fmt.Println(refName, window.Clipped[i], window.Total[i],p)
				}
			}
		}
	}
}

func main(){
	var windowsize int = 100
	var minclipping int = 200
	var mincoverage int64 = 16
	var limit float64 = 0.05
	var storage map[string]Window

	storage = make(map[string]Window)

	file, e := os.Open(os.Args[1])
	if e != nil {
		fmt.Println("Cannot open BAM file")
		return
	}
	b, e := bam.NewReader(file, 0)
	if e != nil {
		fmt.Println("Cannot open BAM file")
		file.Close()
		return
	}

	for {
		rec, e := b.Read()
		if e == io.EOF {
			break
		}
		_,exists := storage[rec.Ref.Name()]
		if !exists {
			var actual Window
			actual.Clipped = make([]int64, (rec.Ref.Len() / windowsize) + 1)
			actual.Total = make([]int64, (rec.Ref.Len() / windowsize) + 1)
			storage[rec.Ref.Name()] = actual
		}

		if len(rec.Cigar) > 0 && rec.Flags & 2048 == 0 && rec.Flags & 256 == 0 && (rec.Cigar[0].Type() == sam.CigarSoftClipped || rec.Cigar[0].Type() == sam.CigarHardClipped) {
			if rec.Cigar[0].Len() > minclipping {
				//fmt.Println(rec.Ref.Name(), rec.Pos, rec.Cigar[0].Len())
				storage[rec.Ref.Name()].Clipped[rec.Pos / windowsize]++
			}
		}
		storage[rec.Ref.Name()].Total[rec.Pos / windowsize]++
	}

	findProblems(storage, mincoverage, limit)
}
