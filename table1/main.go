// Copyright (c) 2026 G. K. Bagger, A. R. Booker, B. Kerr, K. J. McGown,
// V. Starichkova and T. Trudgian.
// Released under the MIT License; see LICENSE.
//
// Table 1 of the paper.  Two modes:
//
//	go run .                  certify each row of Table1.csv at f = f_0(q_1)
//	go run . construct [n]    reconstruct the table from scratch, scanning
//	                          test points a*10^b with a in [1,10) in steps
//	                          of 1/n (default n = 100, as in the paper)
//
// The companion script certify_table1_range.py checks the other half of the
// claim, namely that each row continues to hold for all f up to 10^22.
package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"math/big"
	"os"
	"strconv"
)

func main() {
	mode := "certify"
	if len(os.Args) > 1 {
		mode = os.Args[1]
	}
	switch mode {
	case "certify":
		certify()
	case "construct":
		steps := 100
		if len(os.Args) > 2 {
			var err error
			if steps, err = strconv.Atoi(os.Args[2]); err != nil || steps < 1 {
				log.Fatalln("usage: go run . construct [steps (positive int)]")
			}
		}
		construct(steps)
	default:
		log.Fatalln("usage: go run . [certify | construct [steps]]")
	}
}

// certify reads Table1.csv and checks condition (10) of Theorem 3.3, together
// with the hypotheses of that theorem and condition (9), at f = f_0(q_1).
func certify() {
	file, err := os.Open("Table1.csv")
	if err != nil {
		log.Fatalln(err)
	}
	defer file.Close()

	cases, err := csv.NewReader(file).ReadAll()
	if err != nil {
		log.Fatalln(err)
	}

	allok := true
	for i, record := range cases {
		if i == 0 {
			continue // header
		}
		val1, err := strconv.Atoi(record[0])
		if err != nil {
			log.Fatalln(err)
		}
		q1 := uint64(val1)
		f, err := strconv.ParseFloat(record[1], 64)
		if err != nil {
			log.Fatalln(err)
		}
		lambda, err := strconv.ParseFloat(record[2], 64)
		if err != nil {
			log.Fatalln(err)
		}
		vars := VariableStatev2New(q1, math.Sqrt(f), lambda, 3)
		ok := vars.ValidateStatev2()
		allok = allok && ok
		fmt.Printf("q1=%3v  f_0=%.2e  lambda=%.2f  h=%v  ok=%v\n", q1, f, lambda, vars.h, ok)
	}
	fmt.Println("\nALL ROWS PASS AT f_0:", allok)
	if !allok {
		os.Exit(1)
	}
}

// construct reproduces the (q_1, f_0) column of Table 1: for each test point
// it reports the largest prime q_1 admissible there, printing a line whenever
// that value increases.
func construct(steps int) {
	oldq1 := uint64(3)
	for b := 14; b <= 20; b++ {
		for i := 0; i < 9*steps; i++ {
			if b == 20 && i >= steps {
				return // test points stop at 2*10^20
			}
			a := 1.0 + float64(i)/float64(steps)
			sqrtF := math.Sqrt(math.Ceil(a * math.Pow10(b)))
			currentq1 := SmallestFreeWinVaryingLambda(sqrtF, 3, 1000, steps)

			if currentq1 > oldq1 {
				oldq1 = currentq1
				if big.NewInt(int64(currentq1)).ProbablyPrime(2) {
					fmt.Printf("q1=%3v  f_0=%.2f*10^%v\n", currentq1, a, b)
				}
			}
		}
	}
}
