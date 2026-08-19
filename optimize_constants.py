#!/usr/bin/env python3
# Copyright (c) 2026 G. K. Bagger, A. R. Booker, B. Kerr,
#                    K. J. McGown, V. Starichkova and T. Trudgian
# Released under the MIT License; see LICENSE.
#
# Joint optimization of the parameters behind Theorem 7.1 and Proposition 4.2,
# and the scan that locates the analytic floor f >= 1e22.
"""Joint optimization of the analytic threshold F:
   thm:main1 (K1, theta, rho) x Prop 4.2 (r, C, M) coupled via Case II.
   Lemma set: R-S two-sided (x>=286), BJS lower bound exact for x<=1e12,
   Dusart pi(x) bound (x>=355991), Ma q2-bound, q1 >= 379 hypothesis."""
from mpmath import mp, mpf, log, sqrt, pi, exp, ceil, findroot

mp.dps = 30
PI26 = pi**2/6
C12 = 12/pi**2
AL = sqrt(mpf(3)/2)
Q1MIN = 379

def E_UV(U, V):
    lU = log(U)
    return mpf(6)/pi**2*(1/U+1/V)*lU**2 + (mpf('2.66')/U+mpf('2.26')/V)*lU + mpf('7.57')/U + mpf('9.82')/V

def eps_at(f, K1, th):
    V = ceil(K1*f**mpf('0.375'))
    U = ceil(AL*V/f**mpf('0.25'))
    T = ceil(f**mpf('0.25')/AL)
    Delta = (2*sqrt(mpf(6)))**mpf('0.25')/sqrt(mpf(K1))
    E = E_UV(U, V)
    E1 = Delta*(C12*th*log(U) + 1 + U**th*E)**mpf('0.25')
    E2 = C12*U**(-th) + E
    return E1+E2

def dusart(x):
    l = log(x)
    return (1 + 1/l + mpf('2.51')/l**2)/l

def errors(f, K1, rho):
    V0 = K1*f**mpf('0.375')           # lower bound for V0; all terms decreasing
    Q = V0**rho
    # BJS: Mertens lower bound is exact for x <= 1e12.  By [Ma, Cor 2],
    # q2 <= 1.821 f^{1/4} (log f)^{3/2}; when that is <= 1e12, delta0(q2)=0.
    # Otherwise sup_{x>=Q} delta0(x) = 1/(2 log^2 max(Q,1e12)).
    ma = mpf('1.821')*f**mpf('0.25')*log(f)**mpf('1.5')
    d0 = mpf(0) if ma <= mpf('1e12') else 1/(2*log(max(Q, mpf('1e12')))**2)
    eta = mpf('2.5e-5')
    return 1/mpf(Q1MIN) + 1/Q + d0 + 1/(2*log(V0)**2) + (1+eta)*dusart(V0)

def rho_min(f, K1, th):
    eps = eps_at(f, K1, th)
    rho = mpf('0.60')
    for _ in range(60):
        rhs = 2*(1-eps)/3 - errors(f, K1, rho)
        if rhs <= 0: return None
        newrho = exp(-rhs)
        if abs(newrho-rho) < mpf('1e-12'): break
        rho = newrho
    return rho

# ---- Prop 4.2: minimal M at optimal C ----
def M_min(F, r, C):
    hx = 1/mpf(2*r)
    h = ceil(C*F**hx)
    lamp = C + 1/F**hx
    fact = {3:6, 4:24, 5:120}[r]
    dr = {3: 1+mpf(1)/(6*h), 4: 1+mpf(2)/(3*h), 5: 1+mpf(5)/(3*h)}[r]
    W = (2*r-1) + fact/C**(2*r)*dr*C**r      # f^{1/2}/h^r <= 1/C^r
    W = (2*r-1) + fact/C**r*dr
    F8 = (mpf(2)*h/(h-6))**(2*r)
    # E_u(X): X = M f^{(r+1)/4r}/h -- take M ~ 30 as safe lower bound for X
    X = 30*F**(mpf(r+1)/(4*r))/h
    E = 1 - PI26*(mpf('1.25')+1/X)/X
    return sqrt(PI26*lamp*F8*W/E)

def prop_opt(F, r):
    best = None
    C = mpf('0.8')
    while C < 3:
        M = M_min(F, r, C)
        if best is None or M < best[1]: best = (C, M)
        C += mpf('0.01')
    return best

# ---- feasibility of floor F ----
def check(F, verbose=False):
    best = None
    for r in (3, 4):
        Hx = mpf(r+1)/(4*r)
        C, M = prop_opt(F, r)
        for K1 in range(350, 1101, 25):
            for th10 in range(40, 56, 2):
                th = mpf(th10)/100
                rho = rho_min(F, K1, th)
                if rho is None: continue
                A = ((1+mpf('2.5e-5'))*K1)**rho
                B = mpf('0.375')*rho
                e2 = 1 - 2*B - Hx
                val = 3*A*A*M / F**e2
                if best is None or val < best[0]:
                    best = (val, r, C, M, K1, th, rho, A, B)
    if verbose and best:
        val, r, C, M, K1, th, rho, A, B = best
        print(f"  F={float(F):.3g}: CaseII ratio={float(val):.4f} r={r} C={float(C):.2f} "
              f"M={float(M):.1f} K1={K1} th={float(th):.2f} rho={float(rho):.5f} A={float(A):.2f} B={float(B):.5f}")
    return best[0] < 1 if best else False, best

import sys
print("scan for minimal feasible floor F:")
for Fs in ['1.0e22','1.2e22','1.4e22','1.6e22','1.8e22','2.0e22','2.3e22','2.6e22','3e22','4.6e22']:
    ok, best = check(mpf(Fs), verbose=True)
    print(f"  --> {'FEASIBLE' if ok else 'infeasible'}")
