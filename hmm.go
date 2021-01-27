package main

import (
	"fmt"
	"math/rand"
	"time"
)

type HMM struct {
	obslen int
	statelen int
	B [][]float64  // emission probabilities
	A [][]float64  // transmission probabilities
	Y []string     // observations
	X []string     // hidden states
}

func initHMM(h *HMM, statelen int, obslen int) {
	h.obslen = obslen
	h.statelen = statelen
	h.Y = make([]string, obslen, obslen)
	h.X = make([]string, statelen, statelen)
	h.A = make([][]float64, statelen, statelen)
	h.B = make([][]float64, statelen, statelen)

	for i := 0; i < statelen; i++ {
		h.A[i] = make([]float64, statelen, statelen)
		h.B[i] = make([]float64, obslen, obslen)
	}
}

func rndselect(probs []float64) int {
	state := 0
	p := rand.Float64()
	limit := probs[0]
	for j:= 1; j < len(probs); j++ {
		if p > limit && p < limit + probs[j] {
			state = j
			break
		}
		limit = limit + probs[j]
	}
	return state
}

func generate(h HMM, start []float64, length int) ([]int, []int) {
	rand.Seed(time.Now().UnixNano())
	prevstate := rndselect(start)
	obs := rndselect(h.B[prevstate])

	hidden := make([]int, length, length)
	observ := make([]int, length, length)

	for i:= 0; i < length; i++ {
		prevstate = rndselect(h.A[prevstate])
		obs = rndselect(h.B[prevstate])
		hidden[i] = prevstate
		observ[i] = obs
	}
	return hidden, observ
}

func argmax(T [][]float64, i int, j int, h HMM, obs int) int {
	maxvalue := 0.0
	maxindex := 0;
	for k := 0; k < h.statelen; k++ {
		v := T[k][i] * h.A[k][j] * h.B[j][obs]
		if v > maxvalue {
			maxvalue = v;
			maxindex = k
		}
	}
	return maxindex
}

func viterbi(observations []int, h HMM, start []float64) []int {

	T := make([][]float64, h.statelen, h.statelen)
	for i := 0; i < h.statelen; i++ {
		T[i] = make([]float64, len(observations), len(observations))
		T[i][0] = h.B[i][observations[0]] * start[i]
	}

	for i := 1; i < len(observations); i++ {
		for j := 0; j < h.statelen; j++ {
			k := argmax(T, i-1, j, h, observations[i])
			T[j][i] = T[k][i-1] * h.A[k][j] * h.B[j][observations[i]]
		}
	}
	// backtrack
	path := make([]int, len(observations), len(observations))
	for i:= len(observations)-1; i >= 0; i-- {
		maxvalue := 0.0
		maxindex := 0
		for k:= 0; k < h.statelen; k++ {
			if T[k][i] > maxvalue {
				maxvalue = T[k][i]
				maxindex = k
			}
		}
		path[i] = maxindex
	}
	return path
}

func main() {
	var weather HMM

	initHMM(&weather, 2, 3)
	weather.Y[0] = "walk"
	weather.Y[1] = "shop"
	weather.Y[2] = "clean"

	weather.X[0] = "rainy"
	weather.X[1] = "sunny"

	weather.A[0][0] = 0.7
	weather.A[0][1] = 0.3
	weather.A[1][0] = 0.4
	weather.A[1][1] = 0.6

	weather.B[0][0] = 0.1
	weather.B[0][1] = 0.4
	weather.B[0][2] = 0.5
	weather.B[1][0] = 0.6
	weather.B[1][1] = 0.3
	weather.B[1][2] = 0.1

	startp := []float64 {0.6, 0.4}

	samplesize := 20
	hidden, observ := generate(weather, startp, samplesize)
	solution := viterbi(observ, weather, startp)

	for i:= 0; i < samplesize; i++ {
		fmt.Println(weather.X[hidden[i]], weather.Y[observ[i]], weather.X[solution[i]])
	}
}
