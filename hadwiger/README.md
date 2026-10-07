# Hadwiger–Nelson attempt (7 Oct 2026, ~2 h session)

Goal: a smaller **Moser-spindle-free 5-chromatic unit-distance graph** than the current record
(1441 vertices, Heule), starting from Haugland's heptagon construction (arXiv:2608.04542, 2131 vertices).

* `build.py` — rebuilds the path graphs T5, T6 from the 84 heptagon unit vectors. Reproduces |T6| = 12856, |T5| = 1042.
* `g1.py` — builds G0 and its 7-core G1. Reproduces |G1| = 740 vertices, 3985 edges (as in the paper).
* `variants.py` — stricter selection thresholds: 400 and 364 vertices. **Both turn out 4-colourable with A, B
  equal (SAT, seconds), so they do not work.**
* `cegar.py` — counterexample-guided growth from the 400-vertex variant back towards G1, adding vertices that
  block the current 4-colouring. It grew 400 -> 426 -> 437 -> ... -> 563 vertices and was STILL 4-colourable
  (with A, B the same colour) when stopped (log `cegar_all0.log`).
* `finalsize.py` — computes the size of the final 5-chromatic graph (two rotated copies, then the spindle step)
  from a given G1; reproduces Haugland's 1066 and 2131 exactly.

## Outcome (negative, but informative)

* To beat the 1441-vertex spindle-free record with this construction, the final graph 2|G2| - 1 must be < 1441,
  i.e. |G2| <= 720. Since |G2| >= 2|G1| - 414 (414 = overlap of the two copies for the full G1), any working
  G1 must have at most 567 vertices.
* The counterexample-guided search showed that a 563-vertex subgraph (grown greedily from the stricter
  400-vertex variant) still admits a 4-colouring with A, B equal. So this greedy route cannot reach the record;
  a record via the heptagon framework would need either a different (non-greedy) subset or a different way of
  combining copies (e.g. a cheaper final step than the two-copy + spindle assembly).
* Certifying non-colourability is the bottleneck: CaDiCaL did not prove the 740-vertex G1 non-4-colourable
  (with A, B equal) within ~70 minutes on this 4-core machine, with or without at-most-one clauses.

Nothing here is a verified new record.
* `sat4.py`, `sat_noamo.py`, `satfile.py`, `minimize.py` — SAT checks (CaDiCaL via python-sat) and an
  incremental core-guided minimiser.

Key obstacle: proving NON-4-colourability of the 740-vertex G1 with CaDiCaL did not finish within ~50 minutes
on this machine, so no smaller graph could be certified in the session. Nothing here is a verified new record.
