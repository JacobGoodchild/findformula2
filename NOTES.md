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
