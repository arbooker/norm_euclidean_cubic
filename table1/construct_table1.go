// Copyright (c) 2026 G. K. Bagger, A. R. Booker, B. Kerr, K. J. McGown,
// V. Starichkova and T. Trudgian.
// Released under the MIT License; see LICENSE.
//
// Construction of Table 1: the search over lambda and over test points f_0
// that determines, for each prime q_1, the least threshold at which the
// hypotheses of Theorem 3.3 and condition (10) can be met.

package main

import (
	"errors"
	"log"
	"math"
	"math/big"
)

type VariableState struct {
	H                        *big.Float
	q1, q2, h, u, sigma, phi uint64
	X, E, sqrtF              float64
	Case, r                  int
}

func VariableStateNew(q1, q2 uint64, sqrtF, lambda float64, r int) VariableState {
	bigSqrtF := big.NewFloat(sqrtF)
	F := big.NewFloat(1).Mul(bigSqrtF, bigSqrtF)
	H := big.NewFloat(1).Quo(F, big.NewFloat(float64(q1*q2)))

	var Case int
	var u, sigma, phi uint64
	if q2 < (q1 << 1) {
		H.Quo(H, big.NewFloat(3))
	} else {
		H.Quo(H, big.NewFloat(2))
	}
	h := hSet(sqrtF, lambda, r)

	if q2 < h {
		Case = 1
		u = q1 * q2
		sigma = (q1 + 1) * (q2 + 1)
		phi = (q1 - 1) * (q2 - 1)
	} else if q1 < h {
		Case = 2
		u = q1
		sigma = (q1 + 1)
		phi = (q1 - 1)
	} else {
		Case = 3
		u = 1
		sigma = 1
		phi = 1
	}

	X, _ := big.NewFloat(1).Quo(H, big.NewFloat(float64(h))).Float64()
	E := 1 - (float64(sigma)/4.0+(float64(phi)/float64(u))+(float64(phi)/X))*
		(float64(sigma)/X)*
		math.Pi*
		math.Pi/6.0

	return VariableState{H, q1, q2, h, u, sigma, phi, X, E, sqrtF, Case, r}
}

func (vars VariableState) validateE() bool {
	return vars.E > 0
}

func (vars VariableState) validateX() bool {
	return vars.X > float64(vars.u)
}

func (vars VariableState) validateQ1Q2() bool {
	return vars.q1 < vars.q2
}

func (vars VariableState) dSet() (float64, error) {
	switch vars.r {
	case 1:
		fallthrough
	case 2:
		return 1, nil
	case 3:
		return 1 + float64(1)/(float64(6*vars.h)), nil
	case 4:
		return 1 + float64(2)/(float64(3*vars.h)), nil
	case 5:
		return 1 + float64(5)/(float64(3*vars.h)), nil
	case 6:
		return 1 + float64(10)/(float64(3*vars.h)) + float64(5)/(float64(36*vars.h*vars.h)), nil
	default:
		return 0, errors.New("Invalid choice of r")
	}
}

func hSet(sqrtF, lambda float64, r int) uint64 {
	return uint64(math.Ceil(lambda * math.Pow(sqrtF, 1.0/float64(r))))
}

func (vars VariableState) WSet() (float64, error) {
	d, err := vars.dSet()
	if err != nil {
		return 0, err
	}

	factorial := 1.0 / float64(vars.h)
	for i := 2; i <= vars.r; i++ {
		factorial *= float64(i) / float64(vars.h)
	}
	return float64(2*vars.r-1) + vars.sqrtF*d*factorial, nil
}

func (vars VariableState) validateH() bool {
	if !vars.validateE() || !vars.validateX() || !vars.validateQ1Q2() {
		return false
	} else {
		W, err := vars.WSet()
		if err != nil {
			log.Printf("error: %v\n", err)
			return false
		}
		HSquared := big.NewFloat(1).Mul(vars.H, vars.H)
		RHS := big.NewFloat(1)
		RHS.Mul(RHS, big.NewFloat(1.0/vars.E))
		RHS.Mul(RHS, big.NewFloat(math.Pi*math.Pi/6.0))
		RHS.Mul(RHS, big.NewFloat(float64(vars.sigma)/float64(vars.phi)))
		RHS.Mul(RHS, big.NewFloat(float64(vars.u)))
		RHS.Mul(RHS, big.NewFloat(float64(vars.h)))
		RHS.Mul(RHS, big.NewFloat(vars.sqrtF))
		Temp := big.NewFloat(2.0 * float64(vars.h) / float64(vars.h-uint64(3*(3-vars.Case))))
		Temp.Mul(Temp, Temp)
		Pow := big.NewFloat(1)
		for i := 1; i <= vars.r; i++ {
			Pow.Mul(Pow, Temp)
		}
		RHS.Mul(RHS, Pow)
		RHS.Mul(RHS, big.NewFloat(W))

		return HSquared.Cmp(RHS) == 1

	}
}

func largestPossibleQ2(sqrtF float64) uint64 {
	return uint64(math.Floor(1.821 * math.Sqrt(sqrtF) * math.Pow(2.0*math.Log(sqrtF), 1.5)))
}

func smallestFreeWin(sqrtF, lambda float64, r, maxTest int) uint64 {
	var varsQ2Large, varsQ2Small VariableState
	q2Max := largestPossibleQ2(sqrtF)

	h := hSet(sqrtF, lambda, r) //This is to validateH when q2 ==h

	for q1 := uint64(3); q1 < q2Max; q1 += 2 {
		varsQ2Large = VariableStateNew(q1, q2Max, sqrtF, lambda, r)
		varsQ2Small = VariableStateNew(q1, h, sqrtF, lambda, r)
		if varsQ2Large.validateH() && (q1 >= h || varsQ2Small.validateH()) {
			continue
		} else {
			return q1 - 2
		}
	}
	return uint64(maxTest)
}

func SmallestFreeWinVaryingLambda(sqrtF float64, r, maxTest, steps int) uint64 {
	min := 0.1
	max := 2.0
	length := max - min
	courseness := length / float64(steps)
	q1Boundary := uint64(1)
	var q1Try uint64
	for lambda := min; lambda <= max; lambda += courseness {
		q1Try = smallestFreeWin(sqrtF, lambda, r, maxTest)
		if q1Try > q1Boundary {
			q1Boundary = q1Try
		}
	}
	return q1Boundary
}
