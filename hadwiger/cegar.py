# Counterexample-guided growth: start from a small subset S of the 740-vertex G1 and add vertices
# that kill the current 4-colouring (with colour(A)=colour(B)) until no colouring exists.
import pickle, time, sys, random
import numpy as np
from pysat.solvers import Solver
Q, E, a, b = pickle.load(open('G1.pkl', 'rb'))
nv = len(Q)
seedfile = sys.argv[1]; rule = sys.argv[2] if len(sys.argv) > 2 else 'maxdeg'; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
Qs, Es, as_, bs = pickle.load(open(seedfile, 'rb'))
key = lambda p: (round(p[0]*1e6), round(p[1]*1e6))
pos = {key(p): i for i, p in enumerate(Q)}
S = set(pos[key(p)] for p in Qs)
adj = [set() for _ in range(nv)]
for i, j in E: adj[i].add(j); adj[j].add(i)
var = lambda v, c: 4*v + c + 1
sel = lambda v: 4*nv + v + 1
s = Solver(name='cadical153')
for v in range(nv):
    s.add_clause([-sel(v)] + [var(v, c) for c in range(4)])
    for c in range(4):
        for d in range(c+1, 4): s.add_clause([-var(v, c), -var(v, d)])
for i, j in E:
    for c in range(4): s.add_clause([-var(i, c), -var(j, c)])
s.add_clause([var(a, 0)]); s.add_clause([var(b, 0)])
tri = None
for x in adj[a]:
    for y in adj[a] & adj[x]: tri = (x, y); break
    if tri: break
s.add_clause([var(tri[0], 1)]); s.add_clause([var(tri[1], 2)])
S |= {a, b, tri[0], tri[1]}
rng = random.Random(seed); t0 = time.time(); it = 0
while True:
    it += 1
    if not s.solve(assumptions=[sel(v) for v in S]):
        core = s.get_core(); C = set(l - 4*nv - 1 for l in core if l > 4*nv) | {a, b, tri[0], tri[1]}
        print('UNSAT at |S|=%d (core %d) after %d iterations, %.0fs' % (len(S), len(C), it, time.time()-t0), flush=True)
        S = C; break
    model = s.get_model(); col = {}
    for v in S:
        for c in range(4):
            if model[var(v, c)-1] > 0: col[v] = c
    # candidates: vertices outside S whose S-neighbours use all 4 colours (colouring cannot extend)
    cands = []
    for v in range(nv):
        if v in S: continue
        cs = set(col[u] for u in adj[v] if u in S)
        if len(cs) == 4: cands.append((len(adj[v] & S), rng.random(), v))
    if not cands:
        # colouring extends to all of G1?? then G1 itself colourable - impossible; add max-degree vertex
        rest = [(len(adj[v] & S), rng.random(), v) for v in range(nv) if v not in S]
        cands = rest
    cands.sort(reverse=True)
    k = 1 if rule == 'one' else (len(cands) if rule == 'all' else max(1, len(cands)//8))
    for _, _, v in cands[:k]: S.add(v)
    print('it %d |S|=%d cands %d %.0fs' % (it, len(S), len(cands), time.time()-t0), flush=True)
keep = sorted(S); m = {v: i for i, v in enumerate(keep)}
E2 = [(m[i], m[j]) for i, j in E if i in S and j in S]
pickle.dump((Q[keep], E2, m[a], m[b]), open('cegar_%s_%d.pkl' % (rule, seed), 'wb'))
print('FINAL', len(keep), 'vertices', len(E2), 'edges', flush=True)
