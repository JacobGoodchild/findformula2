# findformula2 — hunting for new exact formulas in general relativity

This is an **AI-assisted** project: the computations, searches and write-ups were
done by Claude Code (Anthropic's AI coding assistant), directed by Jacob Goodchild.

The goal: use high-precision computer experiments to discover *new* exact formulas
about black holes and general relativity, check them carefully, and check
the literature (as far as a web search can) to see whether they are already known.

* `formulas.txt` — verified results only, each with a novelty label
  ("not found" / "probably known" / "known").
* `NOTES.md` — research log: what was tried, dead ends, leads.
* `code/` — the Python scripts (mpmath / sympy / numpy / scipy) used to derive and verify everything.

Units: G = c = 1, and lengths are in units of the black hole mass M unless stated otherwise.

Caveat: "not found" means we could not find it with web searches — not a guarantee that
nobody has written it down before. Please treat the results as candidates for expert checking.

## AI disclosure statement

This project is **AI-assisted**. The research was carried out by Claude (an AI system made by
Anthropic), using the Claude Code software, under the direction of Jacob Goodchild. The AI system:
proposed research directions; wrote and ran all the computer-algebra and numerical code (Python with
SymPy, mpmath, NumPy and SciPy); performed the derivations and verification checks; searched the
literature; and drafted `formulas.txt`, `NOTES.md` and the two papers in `papers/`. Jacob Goodchild
directed the work, set the requirements (verify every result, label novelty honestly, never invent
citations) and is responsible for publishing it. All results were checked against independent
high-precision numerical computations (see `code/verify_all.py` and the `papers/*/verify_*.py`
scripts with their logs). Novelty assessments come from literature searches (see
`papers/LITERATURE_CHECK.md`) and may be incomplete. Please treat the results as candidates for
expert checking.

## Papers

* `papers/paperA/paperA.pdf` — *Universal small-spin formulas for the shadow area, the capture
  cross-section and the spin-down of Kerr-like rotating black holes.*
* `papers/paperB/paperB.pdf` — *Exact results for Kerr black holes: shadow areas and perimeters,
  capture of slow particles, and escape from near extremal horizons.*
* LaTeX sources, verification scripts and their output logs are next to each PDF.

## Highlights (crisp exact results; details and checks in `formulas.txt`)

* **One formula for many black holes** [32]: for any spinning black hole built from a non-spinning
  profile f(r) by the standard Newman–Janis/Azreg-Aïnou recipe (regular, quantum-corrected,
  dark-matter-dressed, …), the small-spin shadow area is
  A = (πr²/f)[1 + a²{(1−cos²θ)[1/(2fr²) − 4/(r²κ)] + cos²θ(2f−1)/(fr²)}], κ = 2f − r²f″,
  evaluated at the photon sphere. Same for particle capture at any speed and the spin-down torque
  −a·sin²θ(1−f)/f. Checked to 11–13 digits on Bardeen and Hayward black holes. A new
  term-by-term trick pushes the shadow formula to a⁴, a⁶ and a⁸ for all such black holes at once
  (reproducing every known Kerr coefficient exactly), and gives the charged (Kerr–Newman) shadow at
  any charge in closed form [31].

* **Shadow size of a spinning black hole, exactly** (EHT-relevant). Edge-on, any spin: one formula
  with complete elliptic integrals [1]; any spin and viewing angle: a single short integral [4];
  rule of thumb: area = 27πM²[1 − (1+cos²θ)a²/18 − …] [5] (the a² and a⁴ terms of this rule were
  already known — Bozza et al. 2006; Kobialko & Gal'tsov 2026; the exact formulas and higher terms were not found).
* **Dark matter (slow particle) capture by a spinning black hole, exactly** [15]:
  σ = [7π + π√(1−a²) + 16√(1+a)·E(2a/(1+a))](M/v)²; at maximal spin (7π + 16√2)(M/v)².
  The matching spin-down torque in closed form [18]; how spin changes capture at any speed and
  direction: factor 1 − a²(r_c + 3cos²θ(4−r_c))/(2r_c²(6−r_c)) [17] (analytically derived here, but it
  turned out to follow from Mach, Momennia & Sarbach, PRD 113, 044068 (2026) — so NOT new).
* **Light escaping from right next to a maximally spinning black hole** [6]:
  P = 1/2 − arcsin(k)/(2π) − k/4 for every latitude and for charged holes (7/24 on the equator is
  known; Yan, Guo & Chen 2021); from the innermost orbit of a charged extremal hole 3/8 + √3/9 at the threshold [8];
  fraction of emitted energy that escapes 3/32 + √3/16 [10].
* **Massive particles** near a maximally spinning hole [19], [24]: nothing launched slower than
  c/√2 escapes from rest near the horizon (exactly 1/8 escape at 0.8c); from the innermost
  orbit the threshold is (3√2−2)/7·c ≈ 0.32c and exactly 1/4 escape at c/2.
* **Shadow circumference** [30]: maximal spin seen side-on has perimeter exactly 18√3·M
  (non-spinning: 6√3π·M); Wei, Liu & Mann (2019) gave 16√3·M, the length without the straight NHEK segment.
* **Open question answered (partly)**: the dark-sky solid angle of a maximally spinning black hole
  for an equatorial observer at any distance [13], [14] (called "intractable" in a 2024 paper).
* Also: charged (Kerr–Newman) and string-theory (Kerr–Sen) versions of several of these
  ([7], [9], [11], [12], [22], [25], [26], [28], [29]).

## Practical formula box: capture of particles by a spinning black hole (G = c = 1)

For particles arriving with speed v at angle θ to the spin axis (C = cos²θ), with
r_c = (4v² − 1 + √(1+8v²))/(2v²) (r_c = 4 for slow particles, 3 for light) and the non-spinning value
σ_S = π r_c² M² / ((r_c − 3)(E² − 1)), E = 1/√(1−v²):

    σ/σ_S = 1 − a²·(r_c + 3C(4 − r_c)) / (2 r_c² (6 − r_c))
              + a⁴·[α8 + β8·C + γ8·C²] / (8 r_c⁴ (6 − r_c)⁵) + O(a⁶)          ([17])
    α8 = −r²(2r³ + 42r² − 549r + 1458),  β8 = 6r(2r⁴ + 2r³ − 239r² + 1302r − 2016),
    γ8 = −10r⁵ − 2r⁴ + 1653r³ − 13266r² + 39744r − 41472      (r = r_c)

(Credit: the a² term and the O(a) spin-down below follow from Mach, Momennia & Sarbach, PRD 113, 044068 (2026);
the a⁴ term was identified numerically here and was only available as fits before — Karydas et al. 2026.)

Exact (all orders in spin) for slow particles side-on ([15]):
σ = [7π + π√(1−a²) + 16√(1+a)·E(2a/(1+a))]·(M/v)², and from any direction as one integral ([23]).
Spin-down: mean swallowed L_z per unit energy = −2a·sin²θ·M/(r_c − 2) + O(a³) ([17]b); exact
side-on slow case in closed form ([18]).

## Results so far (see `formulas.txt` for full details, numbers and checks)

All in units G = c = M = 1. "Not found" = we searched and could not find it published.

| # | What it answers (plain English) | Headline formula | Status |
|---|---|---|---|
| 1 | Size (area) of a spinning black hole's shadow seen side-on, any spin | exact formula with elliptic integrals; small-spin series 27π(1 − a²/18 − a⁴/72 − …); near-max-spin 16π+15√3+(3π/√2)√(1−a)+… | not found (a² term known) |
| 2 | Same, maximum spin | 16π + 15√3 | known (2019) |
| 3 | Max-spin shadow area from any viewing angle | one-line elliptic integral; near-axis series | not found |
| 4 | Shadow area for any spin and any angle | a single compact integral ("master formula") | not found |
| 5 | Rule of thumb: how spin shrinks the shadow, any angle | area ≈ 27π[1 − (1+cos²θ)a²/18 − …] | not found (side-on a² term known) |
| 6 | Chance light escapes from just outside a max-spin (charged) horizon, any latitude | 1/2 − arcsin(k)/2π − k/4 | not found (7/24 special case known) |
| 7 | Side-on shadow of a max-spin *charged* black hole | 16π + 8πa² (a ≤ 1/2), elementary arcsin formula above | not found |
| 8 | Escape chance from the innermost stable orbit, charged case | closed form; 3/8 + √3/9 at the threshold a = 1/√2 | not found (uncharged value known) |
| 9 | Max-spin charged shadow from any angle | one-line elliptic integral | not found |
| 10 | Fraction of *energy* that escapes from a max-spin horizon | 3/32 + √3/16 (lamp at rest); ≈59.56% closed form (orbiting lamp) | not found |
| 11 | Side-on shadow of a spinning *charged* black hole, any spin and charge | one compact integral (elliptic) | not found |
| 12 | Shadow of a spinning charged black hole from any angle; rule of thumb | 1 − Q²/3 − (1+cos²θ)a²/18 − … | not found (charge-only terms known) |
| 13 | How much sky a max-spin black hole blots out for a nearby observer on its equator (open question from 2024) | elementary + one elliptic integral; series in M/r with π, √3 | not found |
| 14 | Same, for charged max-spin black holes | same structure; series with polynomial coefficients | not found |
| 15 | How big a target a spinning black hole is for slow particles like dark matter (side-on) | 7π + π√(1−a²) + 16√(1+a)·E(2a/(1+a)), times (M/v)² | not found |
| 16 | Same, from any direction, and averaged over directions | 16π(M/v)²[1 − a²/16 − 17a⁴/1280 − …]; the a² term doesn't depend on direction | probably known (follows from Will 2012) |
| 17 | How spin changes the capture target for particles of ANY speed and direction | 1 − a²(r_c + 3cos²θ(4−r_c))/(2r_c²(6−r_c)), with r_c from the speed | not found |
| 18 | How fast slow matter spins a black hole down (side-on), exactly | closed form with K and E; swallowed L_z = −aM(1 + 3a²/32 + …) | not found |
| 19 | Can massive particles escape from right next to a max-spin black hole? | only if faster than c/√2; P = (1 − √(1−β²)/β)/2 up to (√3/2)c | not found |
| 20 | Dark-matter capture by a spinning black hole averaged over all directions | one integral; exact elementary value at max spin ≈ 45.28 (M/v)² | not found |
| 21 | Max-spin black hole: capture target for particles of any speed (side-on) | one elliptic family joining 7π+16√2 (slow) to 16π+15√3 (light) | not found |
| 22 | Slow-particle capture by max-spin *charged* black holes (side-on) | one elliptic integral; golden-ratio endpoint πφ⁵ | not found (endpoint probably known) |
| 23 | Exact rule for when a slow star/particle gets swallowed, at any orbital tilt | simple parametric formula replacing published fits | not found in this form |
| 24 | How fast must debris be thrown from the innermost orbit of a max-spin black hole to escape? | faster than (3√2−2)/7·c ≈ 0.32c; exactly 1/4 escape at c/2 | not found |
| 25 | Shadow size of the string-theory (Kerr–Sen) black hole, side-on | compact integral; elementary closed form at extremality | not found |
| 26 | Light escaping from the horizon of an extreme string-theory (Kerr–Sen) black hole | 1/2 − arcsin(k)/2π − k/4, k = (1+b)/2; matches published numbers | not found |
| 27 | How fast a black hole spins down in a bath of slow dark matter | ⟨ℓ_z⟩ = −(2a/3 + a³/15 + …)M; exact at max spin: −0.8193M | not found (exact value) |
| 28 | Max-spin Kerr–Sen shadow from any angle | factorises into two quadratics; one-line elliptic formula | not found |
| 29 | Dark-matter capture by a string-theory (Kerr–Sen) black hole, side-on | 16√(B+a)·E(m) + π(3B + 4 + √(B²−a²)) | not found |
| 30 | Circumference of a spinning black hole's shadow | exact; 18√3·M for maximal spin seen side-on | not found |
| 31 | Capture by a charged, slowly spinning black hole at any speed, direction and charge (and its shadow size) | one-line spin correction; RN ISCO polynomial in the denominator | not found |
| 32 | One formula for the small-spin shadow size, capture cross-section and spin-down of ANY Kerr-like rotating black hole (regular, quantum-corrected, …) from its non-rotating profile f(r) | shadow uses f and f″ at the photon sphere | not found |
