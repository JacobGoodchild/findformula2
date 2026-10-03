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
