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
| 11 | Side-on shadow of a spinning *charged* black hole, any spin and charge | one compact integral (elliptic) | not found |
| 10 | Fraction of *energy* that escapes from a max-spin horizon | 3/32 + √3/16 (lamp at rest); ≈59.56% closed form (orbiting lamp) | not found |
