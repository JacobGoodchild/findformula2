import pickle, time, sys
from sat4 import colourable
f = sys.argv[1]; solver = sys.argv[2] if len(sys.argv) > 2 else 'cadical153'
Q, E, a, b = pickle.load(open(f, 'rb'))
t0 = time.time(); r = colourable(len(Q), E, a, b, solver=solver)
print(f, solver, '|V|=%d' % len(Q), 'colourable with A=B:', r, '%.1fs' % (time.time()-t0), flush=True)
