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
