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

## Highlights (crisp exact results; details and checks in `formulas.txt`)

* **Shadow size of a spinning black hole, exactly** (EHT-relevant). Edge-on, any spin: one formula
  with complete elliptic integrals [1]; any spin and viewing angle: a single short integral [4];
  rule of thumb: area = 27πM²[1 − (1+cos²θ)a²/18 − …] [5]. Previously only fits/numerics.
* **Dark matter (slow particle) capture by a spinning black hole, exactly** [15]:
  σ = [7π + π√(1−a²) + 16√(1+a)·E(2a/(1+a))](M/v)²; at maximal spin (7π + 16√2)(M/v)².
  The matching spin-down torque in closed form [18]; how spin changes capture at any speed and
  direction: factor 1 − a²(r_c + 3cos²θ(4−r_c))/(2r_c²(6−r_c)) [17] (analytically derived).
* **Light escaping from right next to a maximally spinning black hole** [6]:
  P = 1/2 − arcsin(k)/(2π) − k/4 for every latitude and for charged holes (7/24 on the equator,
  known); from the innermost orbit of a charged extremal hole 3/8 + √3/9 at the threshold [8];
  fraction of emitted energy that escapes 3/32 + √3/16 [10].
* **Massive particles** near a maximally spinning hole [19], [24]: nothing launched slower than
  c/√2 escapes from rest near the horizon (exactly 1/8 escape at 0.8c); from the innermost
  orbit the threshold is (3√2−2)/7·c ≈ 0.32c and exactly 1/4 escape at c/2.
* **Open question answered (partly)**: the dark-sky solid angle of a maximally spinning black hole
  for an equatorial observer at any distance [13], [14] (called "intractable" in a 2024 paper).
* Also: charged (Kerr–Newman) and string-theory (Kerr–Sen) versions of several of these
  ([7], [9], [11], [12], [22], [25], [26], [28], [29]).

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
