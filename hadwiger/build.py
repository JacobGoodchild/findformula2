# Rebuild Haugland's heptagon-based graph G1 (arXiv:2608.04542, Section 3).
import numpy as np, math, pickle, sys
from collections import defaultdict
S = 10**7
def key(p): return (round(p[0]*S), round(p[1]*S))
phi = math.pi - math.atan((4*math.cos(math.pi/7) - 1)/math.sqrt(3))
theta = phi/(math.pi/21) - 14
U = []
for j in range(42):
    U.append((math.cos(math.pi*j/21), math.sin(math.pi*j/21)))
    U.append((math.cos(math.pi*(j+theta)/21), math.sin(math.pi*(j+theta)/21)))
U = np.array(U)
A = np.array([0.0, 0.0]); B = np.array([0.0, math.sqrt(3)])
def ball(center, depth):
    d = {key(center): (0, center)}; frontier = [center]
    for t in range(1, depth+1):
        new = []
        for p in frontier:
            for u in U:
                q = p + u; k = key(q)
                if k not in d: d[k] = (t, q); new.append(q)
        frontier = new
    return d
def dist_le(p, k, ballC3, ball0):   # is graph distance from centre C to p <= k ? (k <= 6) via ballC3 + ball0[r]
    # ballC3: dict of points with distance<=3 from C; ball0: list of (dist, vec) for origin ball up to 3
    for (r, x) in ball0:
        if r > k - 3: continue
        e = ballC3.get(key(p - x))
        if e is not None and e[0] + r <= k: return True
    return False
def build(n):
    bA = ball(A, 3); bB = ball(B, 3); b0 = ball(np.zeros(2), 3)
    b0list = sorted(((v[0], v[1]) for v in b0.values()), key=lambda t: t[0])
    def exact_dist(p, bC, maxd):
        k = key(p)
        if k in bC: return bC[k][0]
        for dd in range(4, maxd+1):
            if dist_le(p, dd, bC, b0list): return dd
        return 99
    V = {}
    for (bS, bO) in [(bA, bB), (bB, bA)]:
        for k, (a, p) in bS.items():
            if a > n: continue
            b = exact_dist(p, bO, n - a)
            if a + b <= n: V[k] = p
    dA = {k: exact_dist(p, bA, 6) for k, p in V.items()}
    dB = {k: exact_dist(p, bB, 6) for k, p in V.items()}
    return V, dA, dB
if __name__ == "__main__":
    import time; t0 = time.time()
    V6, dA, dB = build(6)
    T5 = [k for k in V6 if dA[k] + dB[k] <= 5]
    print('|T6| =', len(V6), ' |T5| =', len(T5), ' time %.1fs' % (time.time()-t0))
    pickle.dump((V6, dA, dB), open('T6.pkl', 'wb'))
