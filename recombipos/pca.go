package main

import (
	"math"
)

func znorm(m [][]int64) ([][]float64) {
	sum := make([]float64, len(m[0]))
	avg := make([]float64, len(m[0]))
	std := make([]float64, len(m[0]))

	for i:= 0; i < len(m); i++ {
		for j:= 0; j < len(m[i]); j++ {
			sum[j] += float64(m[i][j])
		}
	}

	for j:= 0; j < len(m[0]); j++ {
		avg[j] = sum[j] / float64(len(m))
	}

	ss := make([]float64, len(m[0]))
	for i:= 0; i < len(m); i++ {
		for j:= 0; j < len(m[i]); j++ {
			ss[j] += (float64(m[i][j]) - avg[j]) * (float64(m[i][j]) - avg[j])
		}
	}

	for j:= 0; j < len(m[0]); j++ {
		std[j] = math.Sqrt(ss[j] / float64(len(m)-1))
	}

	ret := make([][]float64, 0)
	for i:=0; i < len(m); i++ {
		row := make([]float64, len(m[i]))
		for j:=0; j < len(m[i]); j++ {
			row[j] = (float64(m[i][j]) - avg[j]) / std[j]
		}
		ret = append(ret, row)
	}

	return ret
}

func cov(m [][]float64, x int, y int) (float64) {
	var sum float64 = 0
	for i:=0; i < len(m); i++ {
		for j:=0; j < len(m); j++ {
			sum += 0.5 * (m[i][x] - m[j][x]) * (m[i][y] - m[j][y])
		}
	}
	return sum / (float64(len(m)) * float64(len(m)-1))
}

func covmatrix(m [][]float64) ([][]float64){
	d := len(m[0])

	ret := make([][]float64, d)
	for i:=0; i < d; i++ {
		row := make([]float64, d)
		for j:=0; j < d; j++ {
			row[j] = cov(m,i,j)
		}
		ret = append(ret, row)
	}
	return ret
}

/*
func pca(m [][]int64){
	n := znorm(m)
}*/
