package main

import (
	"os"
	"github.com/biogo/hts/bam"
	"io"
	"fmt"
	"github.com/biogo/hts/sam"
)

type Window struct {
	Clipped []int64;
	Total   []int64;
}

func findProblems(storage map[string]Window) {
	for refName,window := range storage {
		l := len(window.Clipped)
		for i:=0; i < l; i++ {
			fmt.Println(refName, window.Clipped[i], window.Total[i])
		}
	}
}

func main(){
	var windowsize int = 100
	var minclipping int = 200
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

		if len(rec.Cigar) > 0 && (rec.Cigar[0].Type() == sam.CigarSoftClipped || rec.Cigar[0].Type() == sam.CigarHardClipped) {
			if rec.Cigar[0].Len() > minclipping {
				//fmt.Println(rec.Ref.Name(), rec.Pos, rec.Cigar[0].Len())
				storage[rec.Ref.Name()].Clipped[rec.Pos / windowsize]++
			}
		}
		storage[rec.Ref.Name()].Total[rec.Pos / windowsize]++
	}

	findProblems(storage)
}
