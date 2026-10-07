# Core-guided + greedy vertex minimisation of a graph that must NOT be 4-colourable with colour(a)=colour(b).
# Each vertex v has a selector literal s_v; its clauses are active only when s_v is assumed true.
import pickle, time, random, sys
from pysat.solvers import Solver
def run(Q, E, a, b, solver='cadical153', seed=0, log=print, deadline=None):
    nv = len(Q)
    var = lambda v, c: 4*v + c + 1
    sel = lambda v: 4*nv + v + 1
    s = Solver(name=solver)
    adj = [set() for _ in range(nv)]
    for i, j in E: adj[i].add(j); adj[j].add(i)
    for v in range(nv):
        s.add_clause([-sel(v)] + [var(v, c) for c in range(4)])
        for c in range(4):
            for d in range(c+1, 4): s.add_clause([-var(v, c), -var(v, d)])
    for i, j in E:
        for c in range(4): s.add_clause([-var(i, c), -var(j, c)])
    # symmetry breaking (kept vertices only): A,B colour 0 and a triangle at A
    s.add_clause([var(a, 0)]); s.add_clause([var(b, 0)])
    tri = None
    for x in adj[a]:
        for y in adj[a] & adj[x]: tri = (x, y); break
        if tri: break
    fixed = {a, b} | set(tri)
    s.add_clause([var(tri[0], 1)]); s.add_clause([var(tri[1], 2)])
    alive = set(range(nv))
    def unsat(active):
        return not s.solve(assumptions=[sel(v) for v in active])
    t0 = time.time()
    assert unsat(alive), 'graph is colourable!'
    core = set(-1 for _ in [])
    c = s.get_core()
    if c: alive = set(l - 4*nv - 1 for l in c if l > 4*nv) | fixed
    log('initial UNSAT %.1fs, core size %d' % (time.time()-t0, len(alive)))
    rng = random.Random(seed); order = sorted(alive - fixed, key=lambda v: rng.random())
    for v in order:
        if deadline and time.time() > deadline: break
        if v not in alive: continue
        trial = alive - {v}
        if unsat(trial):
            c = s.get_core()
            alive = (set(l - 4*nv - 1 for l in c if l > 4*nv) | fixed) if c else trial
            log('removed %d -> %d vertices (%.0fs)' % (v, len(alive), time.time()-t0))
    s.delete(); return alive
if __name__ == "__main__":
    Q, E, a, b = pickle.load(open(sys.argv[1] if len(sys.argv) > 1 else 'G1.pkl', 'rb'))
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hours = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    alive = run(Q, E, a, b, seed=seed, log=lambda m: print(m, flush=True), deadline=time.time()+hours*3600)
    keep = sorted(alive); m = {v: i for i, v in enumerate(keep)}
    E2 = [(m[i], m[j]) for i, j in E if i in alive and j in alive]
    pickle.dump((Q[keep], E2, m[a], m[b]), open('G1min_seed%d.pkl' % seed, 'wb'))
    print('FINAL', len(keep), 'vertices', len(E2), 'edges', flush=True)
