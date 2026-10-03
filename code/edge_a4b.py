from mpmath import mp, mpf, sqrt, matrix, lu_solve, nstr, pi, acos
from fractions import Fraction
import anyv_area2 as AA
mp.dps=50
def coeffs(E,th):
    hs=[mpf('0.003')*k for k in range(1,8)]
    ys=[AA.area(E,h,th) for h in hs]
    M=matrix([[h**(2*j) for j in range(7)] for h in hs]); c=lu_solve(M,matrix(ys))
    return [c[j]/c[0] for j in range(3)]
for rc in [Fraction(15,4),Fraction(7,2),Fraction(10,3),Fraction(13,4),Fraction(18,5),Fraction(16,5)]:
    r=mpf(rc.numerator)/rc.denominator; E=sqrt((r-2)**2/(r*(r-3)))
    c=coeffs(E,pi/2)
    print(rc, nstr(c[1],22), nstr(c[2],22), Fraction(str(nstr(c[2],20))).limit_denominator(10**9), flush=True)
