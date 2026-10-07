import numpy as np, pickle, math, time
from pysat.solvers import Solver
from build import key, A, B
V6, dA, dB = pickle.load(open('T6.pkl', 'rb'))
keys = list(V6); P = np.array([V6[k] for k in keys]); idx = {k: i for i, k in enumerate(keys)}
T5 = np.array([dA[k] + dB[k] <= 5 for k in keys])
def unit_edges(Pts, tol=1e-7):
    # all pairs at distance 1 (grid hashing)
    from collections import defaultdict
    cell = defaultdict(list)
    for i, (x, y) in enumerate(Pts): cell[(math.floor(x), math.floor(y))].append(i)
    E = []
    for i, (x, y) in enumerate(Pts):
        cx, cy = math.floor(x), math.floor(y)
        for dx in (-1, 0, 1, 2) if False else (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in cell.get((cx+dx, cy+dy), []):
                    if j > i and abs(math.hypot(Pts[j][0]-x, Pts[j][1]-y) - 1) < tol: E.append((i, j))
    return E
def kcore(nv, E, k, keep=()):
    adj = [set() for _ in range(nv)]
    for i, j in E: adj[i].add(j); adj[j].add(i)
    alive = set(range(nv)); changed = True
    while changed:
        changed = False
        for v in list(alive):
            if v in keep: continue
            if len(adj[v] & alive) < k: alive.discard(v); changed = True
    return alive
def four_col_AB_same(nv, E, a, b, extra_same=None, solver='cadical153'):
    # SAT: 4-colouring with colour(a)==colour(b). Returns True if satisfiable.
    var = lambda v, c: 4*v + c + 1
    s = Solver(name=solver)
    for v in range(nv):
        s.add_clause([var(v, c) for c in range(4)])
        for c in range(4):
            for d in range(c+1, 4): s.add_clause([-var(v, c), -var(v, d)])
    for i, j in E:
        for c in range(4): s.add_clause([-var(i, c), -var(j, c)])
    s.add_clause([var(a, 0)]); s.add_clause([var(b, 0)])   # wlog both colour 0
    r = s.solve(); s.delete(); return r
if __name__ == "__main__":
    t0 = time.time()
    Eall_T5 = None
    # G0: T5 vertices plus T6 vertices adjacent to >= 7 vertices of T5
    t5idx = np.where(T5)[0]
    from collections import defaultdict
    cell = defaultdict(list)
    for i in t5idx: cell[(math.floor(P[i][0]), math.floor(P[i][1]))].append(i)
    cnt = np.zeros(len(keys), int)
    for i, (x, y) in enumerate(P):
        cx, cy = math.floor(x), math.floor(y)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in cell.get((cx+dx, cy+dy), []):
                    if abs(math.hypot(P[j][0]-x, P[j][1]-y) - 1) < 1e-7: cnt[i] += 1
    G0 = [i for i in range(len(keys)) if T5[i] or cnt[i] >= 7]
    Q = P[G0]; E0 = unit_edges(Q)
    a = G0.index(idx[key(A)]); b = G0.index(idx[key(B)])
    core = kcore(len(G0), E0, 7)
    Vc = sorted(core); m = {v: i for i, v in enumerate(Vc)}
    E1 = [(m[i], m[j]) for i, j in E0 if i in core and j in core]
    print('|G0| =', len(G0), ' |G1| =', len(Vc), ' edges =', len(E1), ' A,B in core:', a in core, b in core)
    pickle.dump((Q[Vc], E1, m[a], m[b]), open("G1.pkl", "wb")); raise SystemExit
    print('4-colourable with A,B same colour?', sat, ' (expect False)  %.1fs' % (time.time()-t0))
    pickle.dump((Q[Vc], E1, m[a], m[b]), open('G1.pkl', 'wb'))
