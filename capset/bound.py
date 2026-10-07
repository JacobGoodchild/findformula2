# Cap set capacity lower bound from a constant-weight admissible set A(n,w) of size s,
# via the extended product construction (Tyrrell 2023; FunSearch 2023; X-evolve 2025):
# cap set in F_3^{6 m n} of size s * b0^(n-w) * b1^w, with b1 = a1^m, b0 = a0*m*a1^(m-1),
# (a0,a1) = (12,112) from 6-dimensional cap sets, m from I(m,m-1).
from math import comb, log, exp
def C(s, n, w, m):
    a0, a1 = 12, 112
    b1 = a1**m; b0 = a0*m*a1**(m-1)
    return exp((log(s) + (n-w)*log(b0) + w*log(b1))/(6*m*n))
def best(s, n, w):
    return max((C(s, n, w, m), m) for m in range(2, 60))
if __name__ == "__main__":
    for s, n, w, label in [(3003, 15, 10, 'I(15,10) full'), (43650, 21, 15, 'A(21,15) X-evolve'),
                           (237984, 24, 17, 'A(24,17) FunSearch'), (comb(10, 5), 10, 5, 'I(10,5) Edel'),
                           (comb(11, 7), 11, 7, 'I(11,7) Tyrrell')]:
        c, m = best(s, n, w); print('%-22s s=%7d  C=%.6f (m=%d)' % (label, s, c, m))
