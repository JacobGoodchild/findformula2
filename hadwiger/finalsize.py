import pickle, math, sys, numpy as np
def rot(P, ang, c):
    c = np.array(c); R = np.array([[math.cos(ang), -math.sin(ang)], [math.sin(ang), math.cos(ang)]])
    return (P - c) @ R.T + c
key = lambda p: (round(p[0]*1e6), round(p[1]*1e6))
def final_size(Q):
    s3 = math.sqrt(3)
    for s1, s2 in [(-1, 1), (1, -1)]:
        V1 = rot(Q, s1*math.pi/3, (-0.5, s3/2)); V2 = rot(Q, s2*math.pi/3, (0.5, s3/2))
        U = {key(p) for p in V1} | {key(p) for p in V2}
        pts = {(-1000000, 0), (1000000, 0), (0, 0)}
        if pts <= U: return len(U), 2*len(U) - 1
    return None
if __name__ == "__main__":
    Q, E, a, b = pickle.load(open(sys.argv[1], 'rb'))
    print(sys.argv[1], '|G1| =', len(Q), ' (|G2|, |G3|) =', final_size(Q))
