package main

import (
	"math"
)

func lnfactorial(num float64) float64 {

	switch {
	case num <= 1.0:
		return 0.0
	case num == 2.0:
		return math.Log(2.0)
	case num == 3.0:
		return math.Log(6.0)
	case num == 4.0:
		return math.Log(24.0)
	case num == 5.0:
		return math.Log(120.0)
	case num == 6.0:
		return math.Log(720.0)
	case num == 7.0:
		return math.Log(5040.0)
	case num == 8.0:
		return math.Log(40320.0)
	case num == 9.0:
		return math.Log(362880)
	case num > 9.0 && num < 120.0:
		return lnsum(num)
	default:
		return stirling(num)
	}
}

func lnsum(num float64) float64 {
	var sum float64 = 0.0

	for i := 2.0; i <= num; i++{
		sum += math.Log(i)
	}

	return sum
}

// stirling's approximation for log factorial
func stirling(n float64) float64 {
	return (n * math.Log(n)) - n + (0.5 * math.Log(2.0 * math.Pi * n))
}

