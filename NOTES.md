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

### Lead 4: boosted emitters (ISCO) -> RESULT [8]
* The escape wedge of [6], boosted with v = 1/2, reproduces the known 54.65% ISCO value to
  25 digits -- strong confirmation of [6].
* Extremal KN: near-horizon circular orbits have v = 1/(2a) = k; ISCO on the horizon iff a >= 1/sqrt2.
* Closed form found by splitting the integral; threshold value 3/8 + sqrt3/9.

### Dead ends this session
* Horizon-averaged escape probability (uniform over horizon area) for extremal Kerr:
  0.15321959342592516335569958794462...; involves elliptic integrals (sqrt of a quartic), no
  neat closed form found by PSLQ.
* Orientation-averaged extremal Kerr shadow area: 75.01636779376937802536268135...;
  PSLQ with pi, sqrt2, sqrt3, logs found nothing.
* Orientation-averaged small-spin series (from [5]): <A>/pi = 27 - 2a^2 - 11a^4/27
  - 100a^6/567 - 13688a^8/137781 - ...  (just a corollary of [5]).
* Extremal Kerr shadow at 30 deg (singular modulus k = sin 15 deg): Pi(1/2|k) survives.

### Lead 3b: extremal KN at any inclination -> RESULT [9]
* Factorisation p,q = 2 +- 2 a sin(theta) (Kerr case was 2 +- 2 sin(theta)).
* Pitfall: my direct-check code first had an extra 1/sin(theta) (ratio = 1/sin40 exactly).
* Braneworld extension (a > M) checked numerically for [6],[7],[8].
* sympy could not do the ISCO integral symbolically; the closed form in [8] rests on 30-digit
  numerical identification at 4 k values + the known k=1/2 value + ray simulation.

### Lead 5: energy-weighted escape -> RESULT [10]
* ZAMO: 3/32 + sqrt3/16 (Kerr) found by PSLQ immediately; general a via PSLQ at 2 values.
* ISCO: by-parts + same elementary integrals. Check: max of E_inf/E_emit = sqrt3 (known max blueshift).

### Lead 6: near-extremal expansions
* Edge-on: A = 16pi+15sqrt3 + (3pi/sqrt2) sqrt(eps) + eps[(8/sqrt3) ln(27/(2eps)) - 6 sqrt3] (numerical ID).
* Any angle: sqrt(eps) coefficient = pi beta_m^2/(sqrt2 sin theta), beta_m = NHEK-line half-length.
* Pitfall: tanh-sinh nodes landing exactly on a root of N -> ZeroDivision; fixed with
  r = r1 + (r2-r1)(1-cos t)/2 substitution.

### Session 1 summary & plan for the overnight run
* `code/verify_all.py` re-checks the headline numbers (all OK at end of session 1).
* Strongest results: [1]/[4]/[11]/[12] (exact shadow areas -> directly usable by EHT-type
  analyses), [6]/[8]/[10] (exact escape probabilities/energy fractions near extremal horizons,
  turning numerical tables of a 2020 PRD into formulas).
* Leads to pursue next:
  1. Prove the near-extremal coefficients in [1](d)/[3] (3 pi/sqrt2; 8/sqrt3; ln(27/2)) analytically.
  2. Escape probability for extremal Kerr at ANY distance (Baines-Visser open problem): genus 2
     in general; check the special radius r* = M sqrt(3 - sin^2 theta) where the curve splits,
     and/or a series in M/r*.
  3. Third-order terms of the Kerr-Newman series [12]; spin from shadow area inversion formula.
  4. Other emitters near extremal horizons (plunging / free-fall) using the exact escape wedge.
  5. Branch out of shadows: e.g. exact results for photon-ring Lyapunov exponents or
     quasinormal-mode (eikonal) quantities averaged over inclination, or exact lensing
     quantities for extremal Kerr-Newman.

### Lead 7: Baines-Visser solid angle at finite distance -> RESULT [13]
* CORRECTION: earlier note about a "special splitting radius r* = M sqrt(3 - sin^2)" was based on
  a sign lost in the PDF text of their eq. (13). Correct form: sin Phi = (w^2 - 2 + s^2)/(2(1+w)s).
  With it, the determinant is 4s(1 + z^2 cos^2) != 0: no splitting at any z for theta != 90.
* At theta = 90 deg one quadratic becomes a perfect square -> elliptic. Solved; series in pi, sqrt3.
* Simulation confirms BV angles are Carter-frame angles (ZAMO numbers differ at finite r).

## Session 2 (overnight)
### Lead 8: capture of slow particles (dark matter) -> RESULTS [15], [16]
* Marginally bound spherical orbits rational in y = sqrt(r); edge-on area reduces to an EVEN
  quartic -> 7pi + pi sqrt(1-a^2) + 16 sqrt(1+a) E(2a/(1+a)). Extremal: 7pi + 16 sqrt2.
* General incidence: degree-8 polynomial (genus 3). Series: a^2 term direction independent (!).
* Dead end: capture at general speed v (0 < v < 1), edge-on: with r = 1/(t^2 - p^2) everything is
  rational; Q ~ (t+E)^2 * sextic -> genus 2. Only v -> 0 and v = 1 are elliptic.
### Lead 7b: extremal KN equatorial solid angle -> RESULT [14] (K = 4r^2 identity).
* Dead end: slow capture for Kerr-Newman (edge-on): with r = y^2, Q ~ a^2(y^2-q) - (y^3 - 2y^2 + q)^2
  (sextic) -> genus 2 for q != 0. Only Kerr (q = 0) has the y^2 factor that makes it elliptic.
* Literature: Will (2012, arXiv:1208.3931) critical L_c(i) series for E = 1 reproduces [16] exactly
  -> [16] relabelled "probably known/derivable". Analytic perturbation for general E running
  (code/anyv_pert2.py, rational parametrisation r_c = (u^2+3)^2/(4u^2)).
* Dead end (caught by numerical check): tried the "average over orbit tilt" trick of [20] for the
  PHOTON shadow (would give an elementary extremal value 75.1198...). Wrong: for photons the
  impact radius^2 is C + a^2 cos^2(theta_obs), which depends on the viewing direction as well as
  on the tilt, so the reduction fails (true value 75.016367793769378...). For slow particles the
  extra term is a^2 p^2 cos^2 -> 0, which is why [20] is fine (and brute-force verified).
* [19] massive escape from extremal horizon: threshold c/sqrt(1+a^2), elementary for beta <= sqrt3/2.
### Lead 9: Kerr-Sen (string BH) -> [25] shadow edge-on, [26] near-horizon escape, [28] any-angle extremal shadow, [29] slow capture + torque
* Kerr-Sen slow capture is elliptic (even quartic after r = y^2 - b), unlike Kerr-Newman (genus 2).
* Near-horizon escape structure is universal: k = (lambda_H - a)/2 for Kerr, KN, Kerr-Sen (M=1).
### Lead 10: higher orders of [17]
* Analytic perturbation (anyv_pert2.py) got L_1..L_3 quickly; L_4 very slow in sympy (running).
* Order-a^3 torque with exact impact-plane mapping (torque_any.py) too slow in sympy -> stopped.
* Numerical route: polar a^4 coefficient = -4/(r_c^4 (6 - r_c)) (exact rationals at 8 speeds).
  Edge-on / mixed C running (edge_a4.py).
* Dead end: orientation-averaged extremal photon shadow (75.0163677937693780253...) -- PSLQ with
  pi, sqrt2, sqrt3, sqrt6, ln(1+sqrt2), ln(2+sqrt3), ln((1+sqrt2+sqrt6)/(1+sqrt2)) finds nothing.
* a^4 for general direction & speed: edge-on numerical coefficients (r_c = 15/4, 7/2, 10/3) did not give
  clean rationals at ~13-digit precision; not recorded. Polar a^4 = -4/(r_c^4 (6-r_c)) recorded in [17].
* a^4 edge-on coefficient (any speed): with analytic dL/dr and 50 digits the rationals are exact:
  alpha(r_c) = -(2 r^3 + 42 r^2 - 549 r + 1458)/(8 r^2 (6-r)^5)  [6-point interpolation predicts the 7th exactly].
  Running cos^2 = 1/2 to get the full C-dependence (mid_a4b.py).
* Dead end: inner shadow of polar-viewed extremal Kerr: b_in = 2.03468720358978884867..., no closed form found.
* [17](b') torque to a^3: C-dependence from C = 1/2 data fitted (3 params, 7 points), blind check at C = 1/4, 4/5 (6/6 to 1e-22).
  Isotropic corollary; slow limit reproduces [27] series -(2a/3 + a^3/15) -- independent cross-check.
### Lead 11: Kerr-Newman version of [17] -> RESULT [31]
* First sympy attempt with sqrt symbols produced Abs() junk; fix: keep E, L0 as symbols, reduce powers with E^2, L0^2.
* Forgot the a^2 cos^2 theta impact-plane shift at first (q=0 limit did not match [17]); adding it fixed it.
* Light limit gives compact KN shadow area at any charge; matches [12] double series.
### Lead 12: universal f(r) version -> RESULT [32]
* Same perturbation with Delta = r^2 f + a^2, f Taylor-expanded (f''' drops out). Light & massive.
* Shadow a^2 coefficient depends only on f and kappa = 2f - r^2 f'' (Lyapunov exponent); massive on J = ISCO combo.
* Verified on Bardeen / Hayward via exact contour (Tsukamoto form) and exact capture boundary.
* Idea for next: a^4 term universally? (needs order-4 perturbation; Kerr a^4 from [5]/[17] as check.)
