#!/usr/bin/env python3
# Copyright (c) 2026 G. K. Bagger, A. R. Booker, B. Kerr,
#                    K. J. McGown, V. Starichkova and T. Trudgian
# Released under the MIT License; see LICENSE.
#
# Table 1 records, for each prime q_1 <= 199, a threshold f_0(q_1) and a
# parameter lambda for which condition (10) of Theorem 3.3 holds.  The Go
# programs in this directory construct the table and check each row at
# f = f_0(q_1).  This script checks the remaining half of the claim: that
# each row continues to hold for every f in [f_0(q_1), 10^22].
#
# The left-hand side of (10) is piecewise smooth in f, with jumps at the
# points where h = ceil(lambda f^{1/6}) increments.  It is decreasing on
# each piece, so it suffices to look at the left endpoint of every interval
# of constancy of h.  We do that exhaustively.
import csv, os
from mpmath import mp, mpf, log, ceil, sqrt, pi
mp.dps = 30

TOP = mpf('1e22')
CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Table1.csv')

def Q(f):
    return mpf('1.821')*f**mpf('0.25')*log(f)**mpf('1.5')

def denom(q1):
    """c*q_1 in H = f/(c q_1 Q(f)); see the proof of Theorem 4.1."""
    return 3 if q1 == 2 else (6 if q1 == 3 else 2*q1)

def L(f, q1, h):
    """LHS of (10) with h supplied explicitly (r=3, u=q1, v=q2, ell=1)."""
    H = min(f/(denom(q1)*Q(f)), sqrt(f*h/2)); X = H/h
    sig = q1+1; phi = q1-1
    E = 1 - (pi**2/6)*sig*(mpf(sig)/4 + mpf(phi)/q1 + mpf(phi)/X)/X
    if E <= 0: return None, None, None
    W = 5 + sqrt(f)/h**3*6*(1 + 1/(6*h))
    return (1/E)*(pi**2/6)*(mpf(sig)/phi)*q1*h*sqrt(f)*(2*h/(h-3))**6*W/H**2, X, E

def cond9(f, q1, h):
    """Condition (9) at f, and uniformly for larger f."""
    if sqrt(f)/(73*q1*log(f)) < max(h, 2*q1): return False
    lam = mpf(h)/f**(mpf(1)/6)
    return f**(mpf(1)/3) >= 73*q1*(lam + f**(-mpf(1)/6))*log(f)

with open(CSV) as fp:
    rows = [(int(r['q1']), mpf(r['f']), mpf(r['lambda'])) for r in csv.DictReader(fp)]

print(" q1     f_0      lambda   L(f_0)    #h-intervals  worst jump ratio  max over endpoints")
certified = True
for q1, f0, lam in rows:
    h0 = int(ceil(lam*f0**(mpf(1)/6)))
    v0, X0, E0 = L(f0, q1, h0)
    if v0 is None or v0 >= 1 or not cond9(f0, q1, h0):
        certified = False
        print(f"{q1:4d}  ROW FAILS AT f_0"); continue
    hmax = int(ceil(lam*TOP**(mpf(1)/6)))
    worst = 0; mx = v0; nj = 0
    for h in range(h0+1, hmax+1):
        fj = (mpf(h-1)/lam)**6      # left endpoint of the interval where ceil(lam f^{1/6}) = h
        if fj > TOP: break
        a, _, _ = L(fj, q1, h)      # value just after the jump
        b, _, _ = L(fj, q1, h-1)    # value just before
        if a is None:
            certified = False; print("  E<=0 at", q1, h); break
        nj += 1
        r = a/b
        if r > worst: worst = r
        if a > mx: mx = a
    if mx > v0: certified = False
    print(f"{q1:4d}  {float(f0):8.2e}  {float(lam):.2f}  {float(v0):8.6f}   {nj:6d}"
          f"        {float(worst):.9f}     {float(mx):8.6f}")

# The jump ratios are occasionally a hair above 1, so the LHS is not globally
# monotone; what is certified -- and all that is needed -- is that no left
# endpoint of an interval of constancy of h carries a value exceeding L(f_0).
print("\nSUP ATTAINED AT f_0 FOR EVERY q1:", certified)
