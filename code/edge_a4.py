from mpmath import mp, mpf, sqrt, matrix, lu_solve, nstr, pi, acos
from fractions import Fraction
import anyv_area as AA
mp.dps=40
def coeffs(E,th):
    hs=[mpf('0.004')*k for k in range(1,7)]
    ys=[AA.area(E,h,th) for h in hs]
    M=matrix([[h**(2*j) for j in range(6)] for h in hs]); c=lu_solve(M,matrix(ys))
    return [c[j]/c[0] for j in range(3)]
for rc in [Fraction(15,4),Fraction(7,2),Fraction(10,3)]:
    r=mpf(rc.numerator)/rc.denominator; E=sqrt((r-2)**2/(r*(r-3)))
    for C in [mpf(0),mpf(1)/2]:
        c=coeffs(E,acos(sqrt(C)))
        print(rc, nstr(C,2), nstr(c[1],18), nstr(c[2],16), Fraction(str(nstr(c[2],14))).limit_denominator(10**7), flush=True)
