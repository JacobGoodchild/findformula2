# Research log (AI-assisted, Claude Code; directed by Jacob Goodchild)

## Session 1 — 2026-10-03

### Lead 1: exact area of the Kerr shadow  -> RESULTS [1], [2]
* Motivation: EHT measures shadow sizes. Literature check found a Jan-2026 paper
  (arXiv:2601.09655) that *fits* the area with a 15-parameter surrogate, and a Jul-2026
  paper (arXiv:2607.11979) with only the O(a^2) term. So an exact formula is useful.
* Extremal edge-on case is elementary: 16 pi + 15 sqrt(3).
* General spin edge-on: after noticing dxi/dr = -2((r-1)^3+1-a^2)/(a(r-1)^2), sympy
  reduction gives integrand (15r^2-21r+6 + 6(1-a^2)/(r-1))/sqrt(W). Complete elliptic K,E,Pi.
* Small-spin series has rational coefficients (times 27 pi).
* Pitfall hit: importing a module that sets mp.dps reset my precision -> fake
  "non-analytic" series. Fixed by setting dps after import.
* Non-equatorial observers: beta^2 numerator is a SEXTIC in r -> genus-2 (hyperelliptic)
  in general, so no K/E/Pi form expected. BUT for extremal spin it drops to a quartic.
  -> next lead.

### Lead 1b: extremal Kerr at any inclination -> RESULT [3]
* beta^2 polynomial for a=1 factorises: -(r^2-2(1+s)r-c^2)(r^2-2(1-s)r-c^2). Spotted by
  sympy factoring at c^2 = 3/4, 1/2, 1/3.
* Reduction gives A = 4 Int (3r^2 - cos^2)/sqrt(R) + boundary term (for theta > 47.06 deg).
* Near-axis series in sin(theta) has coefficients sqrt2*pi*(rational, power-of-2 denominator).
* Tried PSLQ at 30 deg (singular modulus k = sin 15 deg) to get pure Gamma(1/3) form: failed
  -- a complete Pi(1/2|k) term survives (third-kind part with residue at infinity).

### Literature correction
* Entry [2] (16 pi + 15 sqrt3) turned out to be KNOWN: Cunha-Herdeiro-Radu 2019
  (arXiv:1909.08039). Relabelled. Lesson: search harder for "special case" values.

### Lead 1c: general (a, theta) -> RESULTS [4], [5]
* Edge polynomial N(r) is degree 6 -> genus 2. Tested the "involution" criterion for split
  Jacobian numerically: not satisfied -> no K/E/Pi closed form expected in general.
* Hermite reduction still gives a remarkably compact single integral (coeffs depend on
  u = a^2 cos^2 theta and a^2).
* Small-spin series at all inclinations with rational coefficients; pole-on case checked
  against exact polar formula. Rule of thumb: area shrinks by (1+cos^2)a^2/18.
* Background job hit "polyroots no convergence" at cos^2=0 (degenerate double root at r=0)
  -> sampled cos^2 at odd sixteenths instead.

### Lead 2: escape cones / escape probability (extremal) -> RESULT [6]
* Baines & Visser (2405.08875) said exact solid angles were intractable. Found: their cos(Theta)
  simplifies, (D - 2(1+w)z(1-z)) = (wz-(1-z))^2, but the general-z solid angle is still a
  genus-2 integral (three quadratics, determinant 1 - 3z^2 + z^2 sin^2 theta != 0).
  Special splitting radius r* = M sqrt(3 - sin^2 theta) noted, NOT pursued yet (lead).
* Far-field limit z -> 0 of their solid angle is A(theta) z^2 = result [3] (checked numerically).
* Near-horizon limit z -> 1 is elementary -> P = 7/24 on the equator. Turned out KNOWN.
  Generalised (own derivation) to all latitudes and to extremal Kerr-Newman -> new formula.
* Simulation pitfalls: (1) forgot photons emitted inward that bounce back (gave 1/6!);
  (2) double-precision cancellation at r-1 ~ 1e-7 (fixed with exact algebra for 1+a^2-a*lambda).

### Lead 3: extremal Kerr-Newman shadow -> RESULT [7]
* Elementary for all a. Threshold a = M/2 = same as escape threshold.
* Shoelace check initially failed: 0/0 at r=1 in the general formulas produced a spike;
  starting at r = 1 + 1e-5 fixed it.
