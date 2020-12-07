package main

import (
	"math"
)

func mean(values []float64) float64 {
	n := float64(len(values))
	s := 0.0
	for i:=0; i < len(values); i++ {
		s += values[i]
	}
	return s / n
}

func variance(values []float64) float64 {
	m := mean(values)
	s := 0.0
	n := float64(len(values))
	for i:=0; i < len(values); i++ {
		s += (values[i] - m) * (values[i] - m)
	}
	return s / (n - 1.0)
}

func sd(values []float64) float64 {
	return math.Sqrt(variance(values))
}

func welch(samplea []float64, sampleb []float64) (float64,float64) {
	va := variance(samplea)
	vb := variance(sampleb)
	na := float64(len(samplea))
	nb := float64(len(sampleb))
	sdelta := math.Sqrt( (va / na) + (vb / nb) )
	if sdelta == 0.0 {
		return 0.0,0.0
	}
	t := (mean(samplea) - mean(sampleb)) / sdelta

	df1 := ((va / na) + (vb / nb)) * ((va / na) + (vb / nb))
	df2 := ((va * va) / (na * na * (na - 1))) + ((vb * vb) / (nb * nb * (nb - 1)))
	df := df1 / df2

	return t, df
}

func tpdf(t float64, v float64) float64 {
	return math.Gamma((v + 1.0) / 2.0) / (math.Sqrt(v * math.Pi ) * math.Gamma(v / 2.0)) * math.Pow(1.0 + (t * t / v), -(v + 1) / 2.0)
}

func ttestpvalue(t float64, v float64) float64 {
	s     := 0.0
	limit := 7.0
	step  := 0.001

	if t < 0.0 {
		t = -t
	}
	for i:= t; i < limit; i+= step {
		s += tpdf(i, v) * step
	}
	return 2.0 * s
}

func binomial(n float64, k float64, p float64) float64 {
	return math.Exp(lnfactorial(n) - (lnfactorial(k) + lnfactorial(n - k))) * math.Pow(p, k) * math.Pow(1.0 - p, n - k)
}

func binomtest(n float64, k float64, p float64) float64 {
	s     := 0.0
	start := k
	end   := n

	if k < n * p {
		start = 1.0
		end   = k
	}

	for i := start; i <= end; i++ {
		s = s + binomial(n, i, p)
	}
	if s > 0.5 {
		s = 0.5
	}
	return 2.0 * s
}
