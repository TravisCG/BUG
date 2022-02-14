package main

import (
	"testing"
)

func TestNorm(t *testing.T) {
	m := make([][]int64, 0)
	for i:=0; i < 10; i++ {
		r := make([]int64, 5)
		r[0] = int64(i)
		r[1] = int64(i+1)
		r[2] = int64(i+2)
		r[3] = int64(i+3)
		r[4] = int64(i+4)
		m = append(m, r)
	}
	n := znorm(m)

	if n[0][0] > -1.4863 {
		t.Errorf("Normalization is incorrect")
	}
}

func TestCov(t *testing.T) {
	m := make([][]float64, 0)
	for i:=0; i < 10; i++ {
		r := make([]float64, 5)
		for j:=0; j < 5; j++ {
			r[j] = float64(i)+float64(j)
		}
		m = append(m,r)
	}
	m[0][0] = 0.7585539
	m[1][0] = 1.1057340
	m[2][0] = 0.8465542
	m[3][0] = 0.3123003
	m[4][0] = -1.5446309
	m[5][0] = -0.5354912
	m[6][0] = -1.7155467
	m[7][0] = 0.1410860
	m[8][0] = 0.4613748
	m[9][0] = 1.3486318

	c := cov(m, 0, 0)
	if c < 1.1267 || c > 1.12674 {
		t.Errorf("Covariance calculation is incorrect")
	}
}
