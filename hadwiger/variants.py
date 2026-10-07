import numpy as np, pickle, math, sys
from collections import defaultdict
from build import key, A, B
from g1 import unit_edges, kcore
V6, dA, dB = pickle.load(open('T6.pkl', 'rb'))
keys = list(V6); P = np.array([V6[k] for k in keys]); idx = {k: i for i, k in enumerate(keys)}
T5 = np.array([dA[k] + dB[k] <= 5 for k in keys])
t5idx = np.where(T5)[0]
cell = defaultdict(list)
for i in t5idx: cell[(math.floor(P[i][0]), math.floor(P[i][1]))].append(i)
cnt = np.zeros(len(keys), int)
for i, (x, y) in enumerate(P):
    cx, cy = math.floor(x), math.floor(y)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cell.get((cx+dx, cy+dy), []):
                if abs(math.hypot(P[j][0]-x, P[j][1]-y) - 1) < 1e-7: cnt[i] += 1
ia, ib = idx[key(A)], idx[key(B)]
out = {}
for thr in [6, 7, 8, 9, 10]:
    G0 = [i for i in range(len(keys)) if T5[i] or cnt[i] >= thr]
    Q = P[G0]; E0 = unit_edges(Q); a = G0.index(ia); b = G0.index(ib)
    for k in [7, 8, 9]:
        core = kcore(len(G0), E0, k)
        if a not in core or b not in core: print('thr', thr, 'k', k, ': A/B dropped'); continue
        Vc = sorted(core); m = {v: i for i, v in enumerate(Vc)}
        E1 = [(m[i], m[j]) for i, j in E0 if i in core and j in core]
        print('thr %2d k %d : |G0| %5d  |core| %4d  edges %5d' % (thr, k, len(G0), len(Vc), len(E1)), flush=True)
        pickle.dump((Q[Vc], E1, m[a], m[b]), open('var_t%d_k%d.pkl' % (thr, k), 'wb'))
