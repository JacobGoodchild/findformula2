from math import log, exp
def Cgen(s, n, w, m, d, a0, a1):
    b1 = a1**m; b0 = a0*m*a1**(m-1)
    return exp((log(s) + (n-w)*log(b0) + w*log(b1))/(d*m*n))
def best(s, n, w, d, a0, a1):
    return max((Cgen(s, n, w, m, d, a0, a1), m) for m in range(2, 80))
ING = [(1270863, 27, 19), (237984, 24, 17), (43650, 21, 15), (3003, 15, 10), (252, 10, 5)]
REC = 2.2203085
if __name__ == "__main__":
    print('check d=6 (12,112):', best(1270863, 27, 19, 6, 12, 112))
    for d, a1s in [(6, [112, 110, 108, 104, 100, 96]), (7, [236, 230, 224, 216, 200]), (8, [512, 496, 480, 448])]:
        for a1 in a1s:
            # minimal a0 to beat the record with the best ingredient
            for a0 in range(1, 3*a1):
                c = max(best(s, n, w, d, a0, a1)[0] for s, n, w in ING)
                if c > REC: print('d=%d a1=%d: need a0 >= %d  (C=%.6f)' % (d, a1, a0, c)); break
            else: print('d=%d a1=%d: impossible' % (d, a1))
