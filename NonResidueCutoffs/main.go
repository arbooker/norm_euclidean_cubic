package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"strconv"
)

func main() {
	file, err := os.Open("Table1.csv")
	if err != nil {
		fmt.Println(err)
		return
	}
	defer file.Close()

	reader := csv.NewReader(file)
	cases, err := reader.ReadAll()
	if err != nil {
		fmt.Println(err)
		return
	}

	for i, record := range cases {
		if i == 0 {
			continue
		}
		val1, err := strconv.Atoi(record[0])
		if err != nil {
			fmt.Println(err)
			return
		}
		q1 := uint64(val1)
		f, err := strconv.ParseFloat(record[1], 64)
		if err != nil {
			fmt.Println(err)
			return
		}
		lambda, err := strconv.ParseFloat(record[2], 64)
		if err != nil {
			fmt.Println(err)
			return
		}
		fmt.Printf("[%v %.1e %v]", q1, f, lambda)
		vars := VariableStatev2New(q1, math.Sqrt(f), lambda, 3)
		fmt.Println(vars.ValidateStatev2())
	}
}

/*
func main() {
	if len(os.Args) != 2 {
		log.Fatalln("Wrong arguments. Usage: go run . <steps (int)>")
	}
	steps, err := strconv.Atoi(os.Args[1])
	if err != nil {
		log.Fatalln("Something went wrong with strconv")
	}
	oldq1 := uint64(7)
	for b := 14; b <= 20; b++ {
		for i := 0; i < 9*steps; i++ {
			if b == 20 && i >= steps {
				os.Exit(0)
			}
			a := 1.0 + float64(i)/float64(steps)
			sqrtF := math.Sqrt(math.Ceil(a * math.Pow10(b)))
			currentq1 := nonresiduecutoffs.SmallestFreeWinVaryingLambda(sqrtF, 3, 1000, steps)

			if currentq1 > oldq1 {
				oldq1 = currentq1
				if big.NewInt(int64(currentq1)).ProbablyPrime(2) {
					fmt.Printf("%v, %.2f*10^%v \n", currentq1, a, b)
				}
			}
		}
	}
}
*/
