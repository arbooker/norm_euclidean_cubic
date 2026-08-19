# Norm-Euclidean cyclic cubic fields

Source code accompanying

> G. K. Bagger, A. R. Booker, B. Kerr, K. J. McGown, V. Starichkova and
> T. Trudgian, *The determination of norm-Euclidean cyclic cubic fields*.

The paper proves that a cyclic cubic field is norm-Euclidean if and only if
its conductor is one of

    7, 9, 13, 19, 31, 37, 43, 61, 67, 103, 109, 127, 157.

Lezowski and McGown had reduced the problem to the range
`2 * 10^14 < f < 10^50`. The paper closes that range in two pieces: an
analytic argument for `f >= 10^22`, and a computation over
`2 * 10^14 < f < 10^22`.

This repository contains that computation, together with the scripts that
verify or produce the explicit constants appearing in the paper. Nothing here
is intended as a certificate: the point is to let an interested reader re-run
the computation, or re-use parts of it, if they wish.

## The computation

```
make
./cubic small     # 2e14 < f < 3.26e16
./cubic large     # 3.26e16 < f < 1e22
```

`cubic.c` enumerates the integers `f = x^2 + 243y^2` in the given range —
every prime conductor of a cyclic cubic field has this form — and discards
those that Criterion 3.1 of the paper rules out. It forks one process per
core, and prints a line `f, x, y` for each surviving candidate. Both ranges
were run to completion and neither printed anything, which is what the paper
asserts.

Two details are worth knowing before reading the source:

* Candidates are enumerated through a wheel, so an `f` whose least prime
  cubic non-residue `q_1` is one of the wheel primes is silently discarded.
  That is justified by Table 1 of the paper, and only above the corresponding
  threshold `f_0(q_1)`; hence the lower end of each range, and hence the split
  into two ranges with different wheels.
* `f` is *not* tested for primality. A composite `f` cannot be a conductor,
  so passing one through the remaining tests is harmless — and cheaper than
  a primality test would be.

The enumeration kernel is selected per range. The small range uses a stride
kernel; the large range uses Bernstein's doubly-focused enumeration, which is
what makes the extension to `10^22` practicable. On a desktop with an Intel
Core i7-8700 the two runs took about 15 hours and 13 days respectively.

## The verification scripts

Each is standalone and takes no arguments; together they need `mpmath` and
`sympy`.

| Script | What it checks |
| --- | --- |
| `certify_table1.py` | The thresholds `f_0(q_1)` and parameters `lambda` of Table 1: that the verification at the tabulated `(f_0, lambda)` extends to every larger `f`, by sweeping the intervals on which `ceil(lambda f^{1/6})` is constant. Note that the quantity involved is *not* globally monotone — see the comment at the foot of the script. |
| `verify_sections3-5.py` | The constants in Theorem 3.3, Lemma 3.9, Theorem 4.1 and Propositions 4.2 and 5.2, including the disjointness hypotheses `X >= 2` and `2HX <= f`. |
| `verify_theorem71.py` | The constants in Theorem 7.1 (`q_2 <= 36.88 f^0.21769`), Proposition 4.2 and the endgame of Section 8. It also regenerates the LaTeX table in the proof of Theorem 7.1, with every entry rounded outward so that the printed table is itself a valid proof. |
| `optimize_constants.py` | Produces those constants, rather than checking them: a joint optimization over the parameters of Theorem 7.1 and Proposition 4.2, coupled through Case II, and the scan that locates the analytic floor at `f = 10^22`. |

## License

MIT; see `LICENSE`.
