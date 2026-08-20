package main

import (
	"log"
	"math"
	"math/big"
)

func VariableStatev2New(q1 uint64, sqrtF, lambda float64, r int) VariableState {
	bigSqrtF := big.NewFloat(sqrtF)
	F := big.NewFloat(1).Mul(bigSqrtF, bigSqrtF)
	h := hSet(sqrtF, lambda, r)
	Q := largestPossibleQ2(sqrtF)
	H1 := big.NewFloat(1).Quo(F, big.NewFloat(float64(2*q1*Q)))
	H2 := big.NewFloat(1).Mul(bigSqrtF, big.NewFloat(math.Sqrt(float64(h/2))))

	H := big.NewFloat(1)
	if H1.Cmp(H2) < 0 {
		H.Set(H1)
	} else {
		H.Set(H2)
	}

	Case := 2
	u := q1
	sigma := (q1 + 1)
	phi := (q1 - 1)

	X, _ := big.NewFloat(1).Quo(H, big.NewFloat(float64(h))).Float64()
	E := 1 - (float64(sigma)/4.0+(float64(phi)/float64(u))+(float64(phi)/X))*
		(float64(sigma)/X)*
		math.Pi*
		math.Pi/6.0

	return VariableState{H, q1, 0, h, u, sigma, phi, X, E, sqrtF, Case, r}
}

func (vars VariableState) Cond2() bool {
	return vars.sqrtF/(float64(73*vars.q1*2)*math.Log(vars.sqrtF)) > math.Max(float64(vars.h), float64(2*vars.q1))
}

func (vars VariableState) ValidateStatev2() bool {
	if !vars.validateE() || !vars.validateX() || !vars.Cond2() {
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
