#!/usr/bin/env python3
# Copyright (c) 2026 G. K. Bagger, A. R. Booker, B. Kerr,
#                    K. J. McGown, V. Starichkova and T. Trudgian
# Released under the MIT License; see LICENSE.
#
# Verification of the constants in Sections 3-5 (Theorem 3.3, Lemma 3.9,
# Theorem 4.1, Propositions 4.2 and 5.2).
"""Checks the numerical claims of Sections 3-5:

- Lemma 3.9 (Polya-Vinogradov): for squarefree s with nu prime factors,
  (s,f)=1, and 1 <= M < f with M >= 2^{nu+1} (s/phi(s)) (sqrt(f) log f + 1),
  there is m <= M with (m,s)=1 and chi(m)=omega.  With s=q1q2 (nu=2,
  s/phi(s) <= 3) this gives M0 = 24 (sqrt(f) log f + 1).
- Small-q2 regime: q1 q2 <= sqrt(f)/(73 log f)  ==>  3 q1 q2 M0 <= f.
- Large-q2 regime: q2 > Y := sqrt(f)/(73 q1 log f) (checked >= max(h,2q1)),
  Theorem 3.3 with u=q1, v=q2, ell=1, r=3,
  H = min( f/(c q1 Q(f)), sqrt(f h / 2) ),  Q(f)=1.821 f^{1/4} log^{3/2} f,
  c = 2 (q2 >= 2q1 guaranteed); for q1=2, denominator 3*Q (criterion
  3 q2 m <= f); for q1=3, denominator 6*Q.
- The disjointness hypotheses X >= 2 and 2HX <= f of Theorem 3.3, and the
  u=1 case.

Section A reads the published Table 1 from table1/Table1.csv; the tools in
that directory construct and certify the table itself.
"""
import csv, os
from math import ceil, sqrt, log, pi
from sympy import primerange
PI26 = pi*pi/6

def Q(f): return 1.821 * f**0.25 * log(f)**1.5

def lhs10(f, q1, h, H, r=3, ell=1):
    X = H/h
    u, su, pu = q1, q1+1, q1-1
    if X < max(u, 2.0): return None
    E = 1 - PI26*su*(su/4 + pu/u + pu/X)/X
    if E <= 0: return None
    W = 2*r-1 + sqrt(f)/h**r * 6 * (1 + 1.0/(6*h))
    return (1/E)*PI26*(su/pu)*u*h*sqrt(f)*(2.0*h/(h-3*ell))**(2*r)*W/H**2

LAMBDAS = [0.1 + 0.002*i for i in range(951)]

def check_thm33(f, q1, denom_c=2.0):
    """Return (ok, bestlhs, Y/h margin at best lambda)."""
    Y = sqrt(f)/(73*q1*log(f))
    best = None
    for lam in LAMBDAS:
        h = ceil(lam * f**(1/6.))
        if h <= q1: continue
        if Y < h or Y < 2*q1: continue          # seam: PV zone must reach h and 2q1
        H = min(f/(denom_c*q1*Q(f)), sqrt(f*h/2.0))
        v = lhs10(f, q1, h, H)
        if v is not None and v < 1:
            if best is None or v < best[0]: best = (v, Y/h)
    return best

_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'table1', 'Table1.csv')
with open(_csv) as _fp:
    table = {int(r['q1']): float(r['f']) for r in csv.DictReader(_fp)}
DENOM = {2: 3.0, 3: 6.0}          # c*q_1 in H = f/(c q_1 Q(f))

print("=== A. Published Table 1 entries (table1/Table1.csv) ===")
worst = (0,None)
bad = []
for q1, f0 in table.items():
    dc = DENOM.get(q1, 2.0)
    r0 = check_thm33(f0, q1, denom_c=dc)
    if r0 is None: bad.append((q1,f0)); continue
    # monotonicity: margins at larger f
    for fmul in (3, 100, 1e4):
        r = check_thm33(f0*fmul, q1, denom_c=dc)
        if r is None or r[0] > r0[0]*1.5: bad.append((q1, f0*fmul)); break
    for fbig in (1e22, 1e30, 1e40, 1e50):
        if fbig > f0 and check_thm33(fbig, q1, denom_c=dc) is None: bad.append((q1,fbig)); break
    if r0[0] > worst[0]: worst = (r0[0], q1)
print("all entries pass:" , not bad, " failures:", bad[:5])
print(f"worst (10)-margin at its f0: LHS={worst[0]:.3f} at q1={worst[1]}")

print("=== B. q1 in {2,3,5,7} at f=1e14 (and above) ===")
for q1, dc in ((2,3.0),(3,6.0),(5,2.0),(7,2.0)):
    for f in (1e14, 1e15, 1e18, 1e22, 1e30, 1e50):
        r = check_thm33(f, q1, denom_c=dc)
        assert r is not None, (q1, f)
        if f == 1e14: print(f"q1={q1}: at 1e14 best LHS={r[0]:.4f}, Y/h={r[1]:.1f}", end="  ")
    print("| passes up to 1e50")

print("=== C. displayed constants for the written q1=2,3 cases (lambda=3/4) ===")
for q1, mult in ((2,3.0),(3,6.0)):
    f = 1e14; h = ceil(0.75*f**(1/6.))
    Hi = f/(mult*Q(f)); Hc = sqrt(f*h/2.0)
    X = min(Hi,Hc)/h
    su, pu = q1+1, q1-1
    E = 1 - PI26*su*(su/4 + pu/q1 + pu/X)/X
    W = 5 + (6+1/h)*(0.75)**-3
    prod = 10*2.1**6*19.24
    condi = sqrt(prod*(mult*1.821)**2)
    print(f"q1={q1}: X={X:.3g}, E={E:.6f} (claim >=0.9995), W={W:.3f}<=19.24, "
          f"10*2.1^6*19.24={prod:.1f}<16510, cond (i): {condi:.0f} <= f^(5/12)/log^1.5 f "
          f"(={f**(5/12.)/log(f)**1.5:.0f} at 1e14)")
print(f"cond (ii): f^(1/2) > 2*16510/0.75 = {2*16510/0.75:.0f} i.e. f > {(2*16510/0.75)**2:.3g}")
print(f"(1/0.9995)*pi^2/6*sigma/phi*u = {(1/0.9995)*PI26*3*2:.4f} (q1=2) "
      f"{(1/0.9995)*PI26*2*3:.4f} (q1=3), both < 10")

print("=== D. PV-lemma arithmetic ===")
# threshold: 2^{nu+1}(s/phi)(sqrt f log f + 1) with nu=2, s/phi<=3 -> 24(...)
# regime q1q2 <= sqrt f/(73 log f): 3*q1q2*M0 <= f?
for f in (1e14, 2.25e19, 1e22):
    M0 = 24*(sqrt(f)*log(f)+1)
    lhs = 3*(sqrt(f)/(73*log(f)))*M0
    print(f"f={f:.0e}: 3*(f^0.5/73logf)*M0/f = {lhs/f:.4f} (<=1), M0<f: {M0<f}")
print("s/phi(s) worst = (2/1)(3/2) =", 2*1.5)

print("=== E. Prop 4.2: E_1 and the 2HX <= f cap ===")
f = 1e22
E1new = 1 - PI26*(0.25+2+1/6.2e5)/6.2e5
final = (1/E1new)*PI26*1.8218*268.4*9.189/86**2
print(f"E1 = {E1new:.7f}; final LHS = {final:.5f} < 0.9994: {final<0.9994}")
print(f"2HX/f = 2*86^2/1.82*f^(-1/2) = {2*86**2/1.82/sqrt(1e22):.2e} <= 1 for f>=1e22; "
      f"= 1 when f = {(2*86**2/1.82)**2:.3g}")

print("=== F. Prop 5.2 Case II with min-H over f-grid ===")
def p52(f):
    h = ceil(1.3*f**(1/6.))
    H = min(f/(2*373*Q(f)), sqrt(f*h/2.0))
    return max(lhs10(f,q1,h,H) for q1 in primerange(5,374))
for f in (1e22, 1e26, 1e35, 1e42, 1e50):
    print(f"  f={f:.0e}: max LHS = {p52(f):.3g}")
