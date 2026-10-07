# Hadwiger–Nelson attempt (7 Oct 2026, ~2 h session)

Goal: a smaller **Moser-spindle-free 5-chromatic unit-distance graph** than the current record
(1441 vertices, Heule), starting from Haugland's heptagon construction (arXiv:2608.04542, 2131 vertices).

* `build.py` — rebuilds the path graphs T5, T6 from the 84 heptagon unit vectors. Reproduces |T6| = 12856, |T5| = 1042.
* `g1.py` — builds G0 and its 7-core G1. Reproduces |G1| = 740 vertices, 3985 edges (as in the paper).
* `variants.py` — stricter selection thresholds: 400 and 364 vertices. **Both turn out 4-colourable with A, B
  equal (SAT, seconds), so they do not work.**
* `cegar.py` — counterexample-guided growth from the 400-vertex variant back towards G1, adding vertices that
  block the current colouring. Status at end of session: see log `cegar_all0.log`.
* `sat4.py`, `sat_noamo.py`, `satfile.py`, `minimize.py` — SAT checks (CaDiCaL via python-sat) and an
  incremental core-guided minimiser.

Key obstacle: proving NON-4-colourability of the 740-vertex G1 with CaDiCaL did not finish within ~50 minutes
on this machine, so no smaller graph could be certified in the session. Nothing here is a verified new record.
