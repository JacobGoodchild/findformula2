import pickle, time, sys
from pysat.solvers import Solver
f = sys.argv[1]; solver = sys.argv[2] if len(sys.argv) > 2 else 'cadical153'
Q, E, a, b = pickle.load(open(f, 'rb')); nv = len(Q)
adj = [set() for _ in range(nv)]
for i, j in E: adj[i].add(j); adj[j].add(i)
var = lambda v, c: 4*v + c + 1
s = Solver(name=solver)
for v in range(nv): s.add_clause([var(v, c) for c in range(4)])
for i, j in E:
    for c in range(4): s.add_clause([-var(i, c), -var(j, c)])
for c in range(1, 4): s.add_clause([-var(a, c)]); s.add_clause([-var(b, c)])
tri = None
for x in adj[a]:
    for y in adj[a] & adj[x]: tri = (x, y); break
    if tri: break
s.add_clause([var(tri[0], 1)]); s.add_clause([var(tri[1], 2)])
t0 = time.time(); r = s.solve()
print(f, solver, 'noAMO |V|=%d' % nv, 'colourable with A=B:', r, '%.1fs' % (time.time()-t0), flush=True)
