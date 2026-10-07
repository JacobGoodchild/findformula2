# 4-colourability test with symmetry breaking: A=B colour 0; a triangle A,x,y gets colours 0,1,2.
import pickle, time, sys
from pysat.solvers import Solver
def adjlist(nv, E):
    adj = [set() for _ in range(nv)]
    for i, j in E: adj[i].add(j); adj[j].add(i)
    return adj
def find_triangle(adj, a):
    for x in adj[a]:
        for y in adj[a] & adj[x]:
            return x, y
    return None
def colourable(nv, E, a, b, solver='cadical153', time_budget=None):
    var = lambda v, c: 4*v + c + 1
    s = Solver(name=solver)
    for v in range(nv):
        s.add_clause([var(v, c) for c in range(4)])
        for c in range(4):
            for d in range(c+1, 4): s.add_clause([-var(v, c), -var(v, d)])
    for i, j in E:
        for c in range(4): s.add_clause([-var(i, c), -var(j, c)])
    adj = adjlist(nv, E)
    s.add_clause([var(a, 0)]); s.add_clause([var(b, 0)])
    t = find_triangle(adj, a)
    if t: s.add_clause([var(t[0], 1)]); s.add_clause([var(t[1], 2)])
    else:
        x = next(iter(adj[a])); s.add_clause([var(x, 1)])
    r = s.solve(); s.delete(); return r
if __name__ == "__main__":
    Q, E, a, b = pickle.load(open('G1.pkl', 'rb'))
    for name in sys.argv[1:] or ['cadical153']:
        t0 = time.time(); r = colourable(len(Q), E, a, b, solver=name)
        print(name, 'satisfiable:', r, '%.1fs' % (time.time()-t0), flush=True)
