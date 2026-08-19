#!/usr/bin/env python3
# Copyright (c) 2026 G. K. Bagger, A. R. Booker, B. Kerr,
#                    K. J. McGown, V. Starichkova and T. Trudgian
# Released under the MIT License; see LICENSE.
#
# Verification of the constants in Sections 6-8 (Theorem 7.1, Proposition 4.2
# and the endgame), and regeneration of the LaTeX table in Theorem 7.1's proof.
"""Verification of the constants in Sections 6-8:
   Theorem 7.1 (thm:main1): K1=500, theta=0.44, rho=0.5805
       => q1 <= q2 <= 36.88 f^{0.21769} for 10^22 <= f <= 10^50, q1 >= 379.
   Prop 4.2 (thm:main2): r=4, C=1.82, m <= 86 f^{5/16} for q1 > 1.82 f^{1/8}.
   Section 8: Cases I/II and 10 q1^2 q2 at f >= 10^22.
   Inputs: Rosser-Schoenfeld Thm 5 (two-sided Mertens, x>=286) and the
   x <= 10^12 lower bound (R-S for x<=10^8, Bordignon-Johnston-Starichkova
   Lemma 16 for the extension); Dusart 1999 pi(x) bound (x >= 355991);
   Ma-McGown-Rhodes-Wanner Cor 2 (q2 <= 1.821 f^{1/4} (log f)^{3/2}).
   Also regenerates the LaTeX table in the proof of Theorem 7.1, with all
   printed entries rounded outward so the table itself is a valid proof."""
from mpmath import mp, mpf, log, sqrt, pi, ceil
import math

mp.dps = 30
PI26 = pi**2/6; C12 = 12/pi**2; AL = sqrt(mpf(3)/2)
K1 = 500; TH = mpf('0.44'); RHO = mpf('0.5805'); ETA = mpf('2.5e-5')
F0 = mpf('1e22')

def E_UV(U, V):
    lU = log(U)
    return mpf(6)/pi**2*(1/U+1/V)*lU**2 + (mpf('2.66')/U+mpf('2.26')/V)*lU + mpf('7.57')/U + mpf('9.82')/V

DELTA = (2*sqrt(mpf(6)))**mpf('0.25')/sqrt(mpf(K1))   # Delta^4 <= 2 sqrt(6)/K1^2

def eps_envelope(fl, fr):
    """upper bound for E1+E2 on [fl,fr]: terms increasing in U at fr, rest at fl"""
    Vmin = K1*fl**mpf('0.375')
    Umin = AL*K1*fl**mpf('0.125')
    Umax = (AL*(K1*fr**mpf('0.375')+1))/fr**mpf('0.25') + 1
    br = C12*TH*log(Umax) + 1 + Umin**TH*E_UV(Umin, Vmin)
    return DELTA*br**mpf('0.25') + C12*Umin**(-TH) + E_UV(Umin, Vmin)

def dusart(x):
    l = log(x)
    return (1 + 1/l + mpf('2.51')/l**2)/l

def err_max(fl, fr):
    """upper bound for the error terms in eq:rho2bound on [fl,fr]"""
    V0 = K1*fl**mpf('0.375')          # V0 >= K1 f^{3/8}; all terms decreasing
    Q = V0**RHO
    ma_right = mpf('1.821')*fr**mpf('0.25')*log(fr)**mpf('1.5')
    d0 = mpf(0) if ma_right <= mpf('1e12') else 1/(2*log(max(Q, mpf('1e12')))**2)
    return 1/mpf(379) + 1/Q + d0 + 1/(2*log(V0)**2) + (1+ETA)*dusart(V0)

CUTS = ['1e22','1.3e22','2e22','6e22','3e23','1e25','1e28','1e33','1e35','1e41','1e50']
TEX = {'1e22':'10^{22}','1.3e22':'1.3\\times10^{22}','2e22':'2\\times10^{22}',
       '6e22':'6\\times10^{22}','3e23':'3\\times10^{23}','1e25':'10^{25}','1e28':'10^{28}',
       '1e33':'10^{33}','1e35':'10^{35}','1e41':'10^{41}','1e50':'10^{50}'}

print(f"Theorem 7.1: need log(1/rho2) > log(1/{float(RHO)}) = {float(log(1/RHO)):.6f}")
ok = True
for i in range(len(CUTS)-1):
    e = eps_envelope(mpf(CUTS[i]), mpf(CUTS[i+1]))
    r = err_max(mpf(CUTS[i]), mpf(CUTS[i+1]))
    eR = math.ceil(float(e)*100000)/100000        # round outward for printing
    rR = math.ceil(float(r)*100000)/100000
    lb = 2*(1-mpf(repr(eR)))/3 - mpf(repr(rR))
    lbR = math.floor(float(lb)*100000)/100000
    ok &= lbR > float(log(1/RHO))
    print(f"$[{TEX[CUTS[i]]},{TEX[CUTS[i+1]]}]$ & ${eR:.5f}$ & ${rR:.5f}$ & ${lbR:.5f}$\\\\")
A = ((1+ETA)*(K1+mpf('1e-8')))**RHO
print(f"intervals {'ALL PASS' if ok else '*** FAIL ***'};  ((1+eta)K1)^rho = {float(A):.4f} <= 36.88: {A <= mpf('36.88')}")
print(f"0.375*rho = {float(mpf('0.375')*RHO):.7f} <= 0.21769: {mpf('0.375')*RHO <= mpf('0.21769')}")
print(f"UV < f check at 1e22: {float(mpf('1.001')*AL*K1**2/sqrt(F0)):.2e} << 1")

# Prop 4.2 with the written-proof roundings
h = 1024; lam = mpf('1.8218'); W = mpf('9.189'); F8 = mpf('268.4'); M = 86
lhs = (1/mpf('0.999996'))*PI26*lam*F8*W/M**2
print(f"\nProp 4.2 (C=1.82, M=86): written-proof bound {float(lhs):.5f} < 0.9994: {lhs < mpf('0.9994')}")
print(f"  checks: f^(1/8)>=562: {F0**mpf('0.125') >= 562}; h=ceil(1.82*562.34..)={int(ceil(mpf('1.82')*F0**mpf('0.125')))};",
      f"lam'={float(mpf('1.82')+1/F0**mpf('0.125')):.5f}<=1.8218;",
      f"W={float(7+24/mpf('1.82')**4*(1+mpf(2)/(3*1024))):.4f}<=9.189;",
      f"F8={float((mpf(2048)/1018)**8):.2f}<=268.4")

# Section 8
B = mpf('0.21769'); Ac = mpf('36.88')
e2 = 1 - 2*B - mpf('0.3125')
print(f"\nSection 8 at f=1e22: CaseII 3*36.88^2*86 = {float(3*Ac**2*M):.0f} <= f^{float(e2):.5f} = {float(F0**e2):.0f}:",
      3*Ac**2*M <= F0**e2)
print(f"  CaseI ratio {float(3*mpf('1.82')*Ac*14*sqrt(log(F0))/F0**(1-(mpf('0.125')+B+mpf(1)/3))):.3e};",
      f"cond3 ratio {float(10*Ac**3/F0**(1-3*B)):.3e}")
