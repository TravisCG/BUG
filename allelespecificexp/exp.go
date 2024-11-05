package main

import (
	"fmt"
	"os"
	"bufio"
	"strings"
	"strconv"
)

func expreader(expname string) map[string][]float64 {
	matrix := make(map[string][]float64)

	// Read the expression file
	expfile, err := os.Open(expname)
	if err != nil {
		fmt.Println("Cannot open expression file")
		return nil
	}

	exp := bufio.NewScanner(expfile)
	exp.Scan() // Read header
	for exp.Scan() {
		line := exp.Text()
		columns := strings.Split(line, "\t")
		geneid := columns[0]
		_, exists := matrix[geneid]
		if !exists {
			matrix[geneid] = make([]float64, 0)
		}
		for i:= 1; i < len(columns); i++ {
			value, err := strconv.ParseFloat(columns[i], 64)
			if err != nil {
				value = 0.0
				println("Cannot convert:", columns[i])
			}
			matrix[geneid] = append(matrix[geneid], value)
		}
	}
	return matrix
}
