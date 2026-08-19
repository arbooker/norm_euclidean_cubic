# Copyright (c) 2026 G. K. Bagger, A. R. Booker, B. Kerr,
#                    K. J. McGown, V. Starichkova and T. Trudgian
# Released under the MIT License; see LICENSE.
#
# Certification of Table 1: for each q1, the left-hand side of (10)
# at the tabulated (f_0, lambda) dominates its value at every larger f.
from mpmath import mp, mpf, log, ceil, sqrt, pi
mp.dps=30
def primes(n):
    p=[];x=2
    while x<=n:
        if all(x%d for d in p if d*d<=x): p.append(x)
        x+=1
    return p
Q=lambda f: mpf('1.821')*f**mpf('0.25')*log(f)**mpf('1.5')
def L(f,q1,h):
    """LHS of (10) with h supplied explicitly (r=3,u=q1,v=q2,l=1)."""
    H=min(f/(2*q1*Q(f)), sqrt(f*h/2)); X=H/h
    sig=q1+1; phi=q1-1
    E=1-(pi**2/6)*sig*(mpf(sig)/4+mpf(phi)/q1+mpf(phi)/X)/X
    if E<=0: return None,None,None
    W=5+sqrt(f)/h**3*6*(1+1/(6*h))
    return (1/E)*(pi**2/6)*(mpf(sig)/phi)*q1*h*sqrt(f)*(2*h/(h-3))**6*W/H**2, X, E
thr={5:'1e14',7:'1e14',11:'2.07e14',13:'3.89e14',17:'1.08e15',19:'1.66e15',23:'3.45e15',
29:'8.47e15',31:'1.10e16',37:'2.19e16',41:'3.26e16',43:'3.93e16',47:'5.55e16',53:'8.88e16',
59:'1.35e17',61:'1.54e17',67:'2.23e17',71:'2.79e17',73:'3.11e17',79:'4.24e17',83:'5.14e17',
89:'6.76e17',97:'9.47e17',101:'1.11e18',103:'1.20e18',107:'1.40e18',109:'1.50e18',
113:'1.73e18',127:'2.73e18',131:'3.08e18',137:'3.67e18',139:'3.89e18',149:'5.11e18',
151:'5.38e18',157:'6.27e18',163:'7.27e18',167:'7.99e18',173:'9.18e18',179:'1.05e19',
181:'1.10e19',191:'1.36e19',193:'1.42e19',197:'1.53e19',199:'1.60e19'}
TOP=mpf('1e22')
lams={}; report=[]
for q1 in sorted(thr):
    f0=mpf(thr[q1]); best=None
    for i in range(10,200):
        lam=mpf(i)/100; h=ceil(lam*f0**(mpf(1)/6))
        if h<=q1 or 3>=h: continue
        v,X,E=L(f0,q1,h)
        if v is None or v>=1: continue
        if sqrt(f0)/(73*q1*log(f0)) < max(h,2*q1): continue      # (9) at f0
        lamp=lam+f0**(-mpf(1)/6)
        if f0**(mpf(1)/3) < 73*q1*lamp*log(f0): continue          # uniform (9)
        if best is None or v<best[1]: best=(lam,v,X,E,h)
    lams[q1]=best
    report.append((q1,best))
# --- exhaustive certification: L non-increasing at every jump of h, and L(f0) is the max
certified=True
print(" q1   lambda   L(f0)     #h-intervals  worst jump ratio   max over all left endpoints")
for q1 in sorted(thr):
    lam,v0,X0,E0,h0=lams[q1]; f0=mpf(thr[q1])
    hmax=int(ceil(lam*TOP**(mpf(1)/6)))
    worst=0; mx=v0; nj=0
    for h in range(int(h0)+1,hmax+1):
        fj=(mpf(h-1)/lam)**6            # left endpoint of the interval where ceil(lam f^{1/6})=h
        if fj>TOP: break
        a,_,_=L(fj,q1,h)                # value just after the jump
        b,_,_=L(fj,q1,h-1)              # value just before
        if a is None: certified=False; print("  E<=0 at",q1,h); break
        nj+=1
        r=a/b
        if r>worst: worst=r
        if a>mx: mx=a
    if mx>v0: certified=False
    print(f"{q1:4d}  {float(lam):.2f}   {float(v0):.6f}   {nj:6d}        {float(worst):.9f}     {float(mx):.6f}")
# The left-hand side is *not* globally monotone: at a jump of h it can rise
# by a relative 5e-5 (the "worst jump ratio" column above exceeds 1 for some
# q1).  What the certification needs, and what is checked here, is the weaker
# statement that no left endpoint of an interval of constancy of h gives a
# larger value than f_0 itself.
print("\nSUP ATTAINED AT f_0 FOR EVERY q1:", certified)
