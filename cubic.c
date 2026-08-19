// Enumeration of prime conductors f in (2e14, 1e22) of cyclic cubic
// fields that are not ruled out as norm-Euclidean by Criterion 3.1 of
//
//   G. K. Bagger, A. R. Booker, B. Kerr, K. J. McGown, V. Starichkova
//   and T. Trudgian, The determination of norm-Euclidean cyclic cubic
//   fields.
//
// Copyright (c) 2026 the above authors.
// Released under the MIT License; see LICENSE.
//
// compile with gcc -O2 -march=native -o cubic cubic.c -lm
//
// unified successor to cubic2.c (singly-focused) and cubic2d.c
// (doubly-focused); the post-processing and table machinery is shared,
// and the enumeration kernel is selected per mode:
//
//   small    (2e14, 3.26e16), wheel 2*3*5*7, stride kernel
//   large    (3.26e16, 1e22), wheel of primes <= 41, doubly-focused kernel
//
// In all cases the admissible x mod M for a given y are the sums
// x = x_1' + x_2' - jM of scaled residues taken from two lists.  When
// M << 2X (stride kernel) each pair is stepped through [-X,X] in
// increments of M; when M > 2X (doubly-focused kernel, after Bernstein)
// each admissible class contains at most one candidate, the per-y lists
// are generated on the fly from four precomputed sublists by the same
// CRT identity, and the candidates are found by sorting the two lists
// and merging, with the sums confined to [0,X] u [M-X,M+X] u [2M-X,2M).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <sys/sysinfo.h>
typedef __int128_t int128_t;
typedef __uint128_t uint128_t;

// Table 1 (MixedBound) values of f_0(q_1) for the primes q_1 <= 199;
// candidates with q_1 > 199 (none are expected) fall through to the
// general criterion test below
#define nq1 46
static double thresholds[nq1] = {
// each of these is exactly representable as a double
1e14, 1e14, 1e14, 1e14, // first four values are not used
2.07e14, 3.89e14, 1.08e15, 1.66e15, 3.45e15, 8.47e15, 1.10e16, 2.19e16,
3.26e16, 3.93e16, 5.55e16, 8.88e16, 1.35e17, 1.54e17, 2.23e17, 2.79e17,
3.11e17, 4.24e17, 5.14e17, 6.76e17, 9.47e17, 1.11e18, 1.20e18, 1.40e18,
1.50e18, 1.73e18, 2.73e18, 3.08e18, 3.67e18, 3.89e18, 5.11e18, 5.38e18,
6.27e18, 7.27e18, 7.99e18, 9.18e18, 1.05e19, 1.10e19, 1.36e19, 1.42e19,
1.53e19, 1.60e19
};

#define nq2 100
static struct {
	int q,unused;
	double negqinv;
	uint128_t thresh;
	int8_t **chi;
} t[nq2];
static int8_t *chiptr[nq2];
static uint64_t B1;
static int nq0;

static inline uint32_t submod(uint32_t a,uint32_t b,uint32_t m) {
	return (a < b) ? a-b+m : a-b;
}

static inline uint32_t addmod(uint32_t a,uint32_t b,uint32_t m) {
	uint32_t s = a+b;
	return (s < a || s >= m) ? s-m : s;
}

static inline uint32_t mulmod(uint32_t a,uint32_t b,uint32_t m) {
	return (uint64_t)a*b % m;
}

static uint32_t modpow(uint32_t x,uint64_t e,uint32_t m) {
	uint32_t res;
	if (!e) return 1;
	for (;!(e&1);e>>=1)
		x = mulmod(x,x,m);
	for (res=x,e>>=1;e;e>>=1) {
		x = mulmod(x,x,m);
		if (e & 1) res = mulmod(res,x,m);
	}
	return res;
}

// replace a and b by x and y such that
// ax+by = gcd(a,b), and return the gcd
static int64_t xgcd(int64_t *a,int64_t *b) {
	int64_t r1,s1,t1,r2,s2,t2,q,z;
	r1 = *a, r2 = *b;
	s1 = 1, s2 = 0;
	t1 = 0, t2 = 1;
	while (r2) {
		q = r1 / r2;
		z = r1-q*r2, r1 = r2, r2 = z;
		z = s1-q*s2, s1 = s2, s2 = z;
		z = t1-q*t2, t1 = t2, t2 = z;
	}
	if (r1 < 0)
		r1 = -r1, s1 = -s1, t1 = -t1;
	*a = s1, *b = t1;
	return r1;
}

static uint32_t modinv(uint32_t x,uint32_t m) {
	int64_t b=m,a=x%m;
	xgcd(&a,&b);
	return (a < 0) ? a+m : a;
}

// assumes a < p and p small
static uint32_t modsqrt(uint32_t a,uint32_t p) {
	for (uint32_t x=0;2*x<=p;x++)
		if (mulmod(x,x,p) == a) return x;
	return 0;
}

// element u+v*sqrt(-3) in F_q[sqrt(-3)] where q == 2 mod 3
typedef struct {
	uint32_t u,v;
} fp2_t;

static inline fp2_t fp2_mul(fp2_t x,fp2_t y,uint32_t q) {
	return (fp2_t){
		.u=((uint64_t)x.u*y.u + (uint64_t)3*(q-x.v)*y.v) % q,
		.v=((uint64_t)x.u*y.v + (uint64_t)x.v*y.u) % q
	};
}

static fp2_t fp2_pow(fp2_t x,uint64_t e,uint32_t q) {
	fp2_t res;
	if (!e) return (fp2_t){.u=1,.v=0};
	for (;!(e&1);e>>=1)
		x = fp2_mul(x,x,q);
	for (res=x,e>>=1;e;e>>=1) {
		x = fp2_mul(x,x,q);
		if (e & 1) res = fp2_mul(res,x,q);
	}
	return res;
}

// table of admissible pairs (x,y) mod q
// assumes 245*(q-1)^2 < 2^64 (which is true in any practical case)
static int8_t **chi_table(uint32_t q) {
	int8_t **f = (int8_t **)calloc(q*(sizeof(int8_t *)+q),1);
	int8_t *cbase = (int8_t *)(f+q);
	uint32_t x,y;
	for (y=0;y<q;y++) f[y] = cbase+(uint64_t)y*q;
	if (q == 2)
		f[0][1] = f[1][0] = 1;
	else if (q == 3)
		f[0][2] = f[1][2] = f[2][2] = 1;
	else if (q % 3 == 1) {
		uint32_t s = modsqrt(q-3,q);
		for (y=0;y<q;y++)
		for (x=!y;x<q;x++) {
			uint32_t pmodq = ((uint64_t)x*x+(uint64_t)243*y*y) % q;
			if (!pmodq) continue;
			uint32_t c = modpow(mulmod(submod(x,mulmod(9*y,s,q),q),pmodq,q),(q-1)/3,q);
			// 1 when chi = 1, 2 when chi = omega, -2 when chi = omega^-1
			if (c == 1)
				f[y][x] = 1;
			else if (addmod(c,c,q) == s-1)
				f[y][x] = 2;
			else
				f[y][x] = -2;
		}
	} else
		for (y=0;y<q;y++)
		for (x=!y;x<q;x++) {
			uint32_t v = fp2_pow((fp2_t){.u=x,.v=9*y%q},((uint64_t)q*q-1)/3,q).v;
			// 1 when chi = 1, 2 when chi = omega, -2 when chi = omega^-1
			f[y][x] = v ? (2*v-q)*2 : 1;
		}
	return f;
}

static inline int modq(int i,int64_t x) {
	int r = x+t[i].q*floor(x*t[i].negqinv);
	return (r < 0) ? r+t[i].q : r;
}

static inline int chi(int i,int64_t x) {
	return (int)chiptr[i][modq(i,x)];
}

// assumes u < 2^96
static inline int is_cube(uint128_t u) {
	uint64_t v = (uint64_t)cbrtl((long double)u);
	uint128_t v3 = (uint128_t)(v*v)*v;
	if (v3 < u) {
		do ++v; while ((v3=(uint128_t)(v*v)*v) < u);
		return (v3 == u);
	}
	if (v3 > u) {
		do --v; while ((v3=(uint128_t)(v*v)*v) > u);
		return (v3 == u);
	}
	return 1;
}

static inline uint64_t intsqrt(uint128_t x) {
	uint64_t s = (uint64_t)sqrtl((long double)x);
	uint128_t s2 = (uint128_t)s*s;
	if (s2 > x) {
		do --s; while ((uint128_t)s*s > x);
		return s;
	}
	if (s2 < x) {
		do ++s; while ((uint128_t)s*s <= x);
		return s-1;
	}
	return s;
}

// return odd part of gcd of x and y
// assumes y > 0 and x != 0
static inline int64_t oddgcd(int64_t x,int64_t y) {
	if (x < 0) x = -x;
	x >>= __builtin_ctzl(x);
	y >>= __builtin_ctzl(y);
	while (x != y)
		if (x > y)
			x -= y, x >>= __builtin_ctzl(x);
		else
			y -= x, y >>= __builtin_ctzl(y);
	return x;
}

static inline uint32_t max(uint32_t x,uint32_t y) {
	return (x < y) ? y : x;
}

static const char *sprint128(uint128_t x) {
	static char buf[40];
	char *s=&buf[39];
	*s = 0;
	if (!x)
		*(--s) = '0';
	else do {
		*(--s) = x % 10 + '0';
		x /= 10;
	} while (x);
	return s;
}

static void print(uint128_t p,int64_t x,int64_t y) {
	char buf[128];
	sprintf(buf,"%s, %ld, %ld\n",sprint128(p),x,y);
	size_t l = strlen(buf);
	size_t n = write(STDOUT_FILENO,buf,l);
	assert(n == l);
}

static void check(int64_t x,int64_t y,int128_t y2) {
	int i,c,c3;
	uint32_t q1,q2,m;
	uint128_t p;
	for (i=nq0;i<nq1;i++)
		if ((c=chi(i,x)) != 1) {
			if (!c) return;
			p = (int128_t)x*x + y2;
			if (p > t[i].thresh || p < B1) return;
			break;
		}
	if (i == nq1) {
		if (oddgcd(x,y) > 1) return;
		p = (int128_t)x*x + y2;
		if (is_cube(p) || p < B1) return;
		for (;i<nq2&&(c=chi(i,x))==1;i++);
		if (i == nq2) { print(p,x,y); return; }
		if (!c) return;
	}
	q1 = (uint32_t)t[i].q;
	while (++i < nq2 && (c=chi(i,x)) == 1);
	if (i == nq2) { print(p,x,y); return; }
	if (!c) return;
	q2 = (uint32_t)t[i].q;
	while (++i < nq2 && (c3=chi(i,x)) != -c)
		if (!c3) return;
	if (i == nq2) { print(p,x,y); return; }
	m = (uint32_t)t[i].q;
	if ((uint64_t)q1*q2*max(3*m,10*q1) > p) print(p,x,y);
}

// assumes res has space for at least one element
static uint32_t crt(uint32_t *res,uint32_t **x,const uint32_t *l,const uint32_t *m,int n) {
	int i;
	uint32_t M,L;
	for (i=0,M=1;i<n;i++) M *= m[i];
	L = 1, *res = 0;
	for (i=0;i<n;i++) {
		uint32_t q=M/m[i],qinv=modinv(q,m[i]),y0;
		for (int j=0;j<l[i];j++) {
			uint32_t y = mulmod(x[i][j],qinv,m[i]) * q;
			if (!j) y0 = y;
			y = submod(j ? y0 : 0,y,M);
			for (int k=0;k<L;k++) res[j*L+k] = submod(res[k],y,M);
		}
		L *= l[i];
	}
	return L;
}

typedef struct {
	uint32_t *x;
	size_t nx;
} crt_list_t;

// number of admissible pairs (x,y) mod q (divided by 3 when q=3)
static inline size_t npairs(uint32_t q) {
	if (q == 2) return 2;
	if (q == 3) return 1; // return 1 when q=3 since x does not depend on y
	return (q % 3 == 1 ? (uint64_t)(q-1)*(q-1) : (uint64_t)q*q-1)/3;
}

// compute all admissible residues x mod M
// and replace each x by A*x mod Mp; for the stride kernel Mp = M and
// A is the inverse of the complementary modulus, while for the
// doubly-focused kernel Mp is the modulus of the half that M belongs
// to and A is its CRT coefficient
static uint32_t get_crt_list(crt_list_t **L,uint32_t M,uint32_t A,uint32_t Mp) {
	int i,n=0;
	uint32_t m[nq0],qindex[nq0],tm=1;
	size_t na=0;
	for (i=0;i<nq0;i++)
		if (!(M % t[i].q)) {
			qindex[n] = i;
			m[n] = t[i].q;
			tm *= t[i].q;
			size_t tn = npairs(t[i].q);
			if (tn > na) na = tn;
			n++;
		}
	assert(tm == M); // ensure M is a squarefree product of small primes

	uint32_t l[n],a[na],*r[n];
	uint32_t Mdiv3 = (M % 3) ? M : M/3;
	crt_list_t *res = (crt_list_t *)malloc(Mdiv3*sizeof(crt_list_t));
	for (uint32_t y=0;y<Mdiv3;y++) {
		res[y].nx = 1;
		for (i=0;i<n;i++) {
			l[i] = 0;
			r[i] = i ? r[i-1]+l[i-1] : a;
			uint32_t ymodq = modq(qindex[i],y);
			for (uint32_t x=0;x<m[i];x++)
				if (t[qindex[i]].chi[ymodq][x] == 1)
					r[i][l[i]++] = x;
			res[y].nx *= l[i];
		}
		if (res[y].nx > 0) {
			res[y].x = (uint32_t *)malloc(res[y].nx*sizeof(uint32_t));
			crt(res[y].x,r,l,m,n);
			for (size_t j=0;j<res[y].nx;j++)
				res[y].x[j] = (uint64_t)A*res[y].x[j] % Mp;
		}
	}
	*L = res;
	return Mdiv3;
}

// LSD radix sort of n uint32 keys, 3 passes of 11 bits;
// b is scratch space of the same size, and the return value
// points at whichever of a,b holds the sorted output
static uint32_t *sort_u32(uint32_t *a,uint32_t *b,size_t n) {
	for (int pass=0;pass<3;pass++) {
		size_t cnt[2048];
		int shift = 11*pass;
		memset(cnt,0,sizeof(cnt));
		for (size_t k=0;k<n;k++) cnt[(a[k]>>shift)&2047]++;
		for (int i=0,s=0;i<2048;i++) {
			size_t c = cnt[i];
			cnt[i] = s, s += c;
		}
		for (size_t k=0;k<n;k++) b[cnt[(a[k]>>shift)&2047]++] = a[k];
		uint32_t *z = a; a = b; b = z;
	}
	return a; // three swaps: sorted data is in the original b
}

static int is_prime(uint32_t x) {
	for (uint32_t p=2;(uint64_t)p*p<=x;p=(p+1)|1)
		if (!(x % p)) return 0;
	return (x > 1);
}

static double get_time(void) {
	struct timeval time;
	gettimeofday(&time,NULL);
	return time.tv_sec + time.tv_usec / 1000000.0;
}

static int128_t B2;
static uint64_t M;                  // the full wheel modulus
static int dfe;                     // 1 for the doubly-focused kernel
static crt_list_t *L1,*L2,*L3,*L4;
static uint32_t ny1,ny2,ny3,ny4;
static uint32_t sc1,sc2;            // stride kernel: scale factors M2, M1
static uint32_t M12,M34;            // doubly-focused kernel: half moduli
static uint32_t *t12,*s12,*t34,*s34; // per-process scratch, grown on demand
static uint64_t *v12,*v34;
static size_t cap12,cap34;

static void grow(size_t n,uint32_t **tb,uint32_t **sb,uint64_t **vb,size_t *cap) {
	if (n <= *cap) return;
	*cap = n + n/4;
	free(*tb); *tb = (uint32_t *)malloc(*cap*sizeof(uint32_t));
	free(*sb); *sb = (uint32_t *)malloc(*cap*sizeof(uint32_t));
	free(*vb); *vb = (uint64_t *)malloc(*cap*sizeof(uint64_t));
	assert(*tb && *sb && *vb);
}

// step through each admissible residue class mod M in [-X,X]
// in increments of M (for M much smaller than X)
static void stride_y(uint64_t y,int64_t X,int128_t Y2) {
	uint32_t y1 = y % ny1;
	uint32_t y2 = y % ny2;
	if (!L1[y1].nx || !L2[y2].nx) return;

	for (int i=nq0;i<nq2;i++)
		chiptr[i] = t[i].chi[modq(i,y)];

	for (size_t i=0;i<L1[y1].nx;i++) {
		uint64_t x1 = (uint64_t)L1[y1].x[i]*sc1;
		for (size_t j=0;j<L2[y2].nx;j++) {
			int64_t x0 = x1 + (uint64_t)L2[y2].x[j]*sc2;
			for (int64_t x=x0;x<=X;x+=M) check(x,y,Y2);
			for (int64_t x=x0-M;x>=-X;x-=M) check(x,y,Y2);
		}
	}
}

// doubly-focused enumeration: each admissible x mod M in [-X,X] appears
// exactly once as v12[i]+v34[j]-{0,M,2M} with the sum in one of the
// windows [0,X], [M-X,M+X], [2M-X,2M)
static void merge_y(uint64_t y,int64_t X,int128_t Y2) {
	uint32_t y1 = y % ny1;
	uint32_t y2 = y % ny2;
	uint32_t y3 = y % ny3;
	uint32_t y4 = y % ny4;
	if (!L1[y1].nx || !L2[y2].nx || !L3[y3].nx || !L4[y4].nx) return;

	for (int i=nq0;i<nq2;i++)
		chiptr[i] = t[i].chi[modq(i,y)];

	grow(L1[y1].nx*L2[y2].nx,&t12,&s12,&v12,&cap12);
	grow(L3[y3].nx*L4[y4].nx,&t34,&s34,&v34,&cap34);

	// combine the pair lists, sort, and scale
	size_t n12=0,n34=0;
	for (size_t i=0;i<L1[y1].nx;i++) {
		uint32_t u1 = L1[y1].x[i];
		for (size_t j=0;j<L2[y2].nx;j++) {
			uint32_t u = u1+L2[y2].x[j];
			t12[n12++] = (u >= M12) ? u-M12 : u;
		}
	}
	for (size_t i=0;i<L3[y3].nx;i++) {
		uint32_t u3 = L3[y3].x[i];
		for (size_t j=0;j<L4[y4].nx;j++) {
			uint32_t u = u3+L4[y4].x[j];
			t34[n34++] = (u >= M34) ? u-M34 : u;
		}
	}
	uint32_t *w12 = sort_u32(t12,s12,n12);
	uint32_t *w34 = sort_u32(t34,s34,n34);
	for (size_t i=0;i<n12;i++) v12[i] = (uint64_t)w12[i]*M34;
	for (size_t j=0;j<n34;j++) v34[j] = (uint64_t)w34[j]*M12;

	for (size_t i=0;i<n12 && v12[i]+v34[0] <= (uint64_t)X;i++)
		for (size_t j=0;j<n34;j++) {
			uint64_t s = v12[i]+v34[j];
			if (s > (uint64_t)X) break;
			check((int64_t)s,y,Y2);
		}
	uint64_t lo2 = M-X, hi2 = M+X, lo3 = 2*M-X;
	for (size_t i=0,j2=n34,j3=n34;i<n12;i++) {
		uint64_t a = v12[i];
		while (j2 > 0 && a+v34[j2-1] >= lo2) j2--;
		for (size_t j=j2;j<n34;j++) {
			uint64_t s = a+v34[j];
			if (s > hi2) break;
			check((int64_t)(s-M),y,Y2);
		}
		while (j3 > 0 && a+v34[j3-1] >= lo3) j3--;
		for (size_t j=j3;j<n34;j++)
			check((int64_t)(a+v34[j]-2*M),y,Y2);
	}
}

int main(int argc,char *argv[]) {
	uint32_t M1,M2,M3,M4;

	// the wheel silently discards candidates whose least nonresidue q_1
	// is a wheel prime, which is justified by Table 1 only for
	// p >= f_0(q_1); hence B1 >= thresholds[nq0-1] is required
	if (argc == 2 && !strcmp(argv[1],"small")) {
		B1 = (uint64_t)2e14;
		B2 = (int128_t)3.26e16; // = f_0(41), where the large range takes over
		M1 = 2*3*5;
		M2 = 7;
		M3 = M4 = 0;
		nq0 = 4;
		dfe = 0;
	} else if (argc == 2 && !strcmp(argv[1],"large")) {
		B1 = (uint64_t)3.26e16; // = f_0(41)
		// the analytic argument (Theorem 7.1 with q2 <= 36.88 f^0.21769
		// and Prop 4.2 with m <= 86 f^{5/16}) covers f >= 1e22
		B2 = (int128_t)1e22;
		// wheel 2*3*5*...*41 split as (2*3*31*37)(11*13*23) x (5*19*41)(7*17*29),
		// balancing the sizes of the two per-y lists
		M1 = 2*3*31*37;
		M2 = 11*13*23;
		M3 = 5*19*41;
		M4 = 7*17*29;
		nq0 = 13;
		dfe = 1;
	} else {
		printf("usage: %s small|large\n",argv[0]);
		return 0;
	}

	double t0 = get_time();
	for (int i=0;i<nq2;i++) {
		if (i) {
			t[i].q = t[i-1].q;
			do t[i].q = (t[i].q+1) | 1; while (!is_prime(t[i].q));
		} else
			t[i].q = 2;
		t[i].negqinv = -1.0/t[i].q;
		t[i].chi = chi_table(t[i].q);
		if (i < nq1) t[i].thresh = (uint128_t)thresholds[i];
	}

	// discarding candidates with q_1 in the wheel must be justified by Table 1
	assert(B1 >= (uint64_t)thresholds[nq0-1]);

	if (dfe) {
		assert((uint64_t)M1*M2 <= 0xFFFFFFFF && (uint64_t)M3*M4 <= 0xFFFFFFFF);
		M12 = M1*M2, M34 = M3*M4;
		M = (uint64_t)M12*M34;
		// CRT coefficients: A1*x1+A2*x2 = x12/M34 mod M12 for x12 = x1 mod M1,
		// x2 mod M2, so that (A1*x1+A2*x2 mod M12)*M34 = x12 mod M12 and
		// 0 mod M34; similarly A3,A4
		uint32_t A1 = modinv(mulmod(M2%M1,M34%M1,M1),M1)*M2;
		uint32_t A2 = modinv(mulmod(M1%M2,M34%M2,M2),M2)*M1;
		uint32_t A3 = modinv(mulmod(M4%M3,M12%M3,M3),M3)*M4;
		uint32_t A4 = modinv(mulmod(M3%M4,M12%M4,M4),M4)*M3;
		ny1 = get_crt_list(&L1,M1,A1,M12);
		ny2 = get_crt_list(&L2,M2,A2,M12);
		ny3 = get_crt_list(&L3,M3,A3,M34);
		ny4 = get_crt_list(&L4,M4,A4,M34);
		// each admissible residue class mod M meets [-X,X] at most once
		assert(2*intsqrt(B2) < M);
	} else {
		M = (uint64_t)M1*M2;
		ny1 = get_crt_list(&L1,M1,modinv(M2%M1,M1),M1);
		ny2 = get_crt_list(&L2,M2,modinv(M1%M2,M2),M2);
		sc1 = M2, sc2 = M1;
	}

	uint32_t nprocs = get_nprocs();
	srandom(42); // ensure all processes generate the same sequence
	setbuf(stdout,NULL);
	double t1 = get_time();
	printf("%s %s: setup took %.3fs, starting %u processes\n",argv[0],argv[1],t1-t0,nprocs);

	for (uint32_t k=0;k<nprocs;k++)
	if (!fork()) {
		// y ranges up to sqrt(B2/243) ~ 6.4e9, so it must be 64-bit
		for (uint64_t b=1;;b+=nprocs) {
			// randomly permute the residue classes
			// to avoid any correlation with nprocs
			uint64_t y = b+((uint64_t)random()+k)%nprocs;
			int128_t Y2 = (int128_t)y*y*243;
			int128_t X2 = B2-Y2;
			if (X2 <= 0) break;
			int64_t X = intsqrt(X2);
			if (dfe) merge_y(y,X,Y2);
			else stride_y(y,X,Y2);
		}
		return 0;
	}

	int status,res=0;
	while (wait(&status) > 0)
		if (WIFSIGNALED(status)) {
			fprintf(stderr,"child process killed by signal %d\n",WTERMSIG(status));
			res = 1;
		} else if (WIFEXITED(status) && WEXITSTATUS(status)) {
			fprintf(stderr,"child process unexpected exit status %d\n",WEXITSTATUS(status));
			res = 1;
		}
	printf("%s %s: sieve took %.3fs, %s exit\n",
		argv[0],argv[1],get_time()-t1,res?"ABNORMAL":"normal");
	return res;
}
