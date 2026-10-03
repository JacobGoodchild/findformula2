from mpmath import mp, mpf, sqrt, matrix, lu_solve, nstr, pi
from fractions import Fraction
from anyv_pole import polar_area
mp.dps=60
def coeffs(E):
    hs=[mpf('0.002')*k for k in range(1,8)]
    ys=[polar_area(E,h) for h in hs]
    M=matrix([[h**(2*j) for j in range(7)] for h in hs]); c=lu_solve(M,matrix(ys))
    return [c[j]/c[0] for j in range(4)]
for rc in [Fraction(15,4),Fraction(7,2),Fraction(10,3),Fraction(13,4),Fraction(31,10),Fraction(18,5),Fraction(37,10),Fraction(16,5)]:
    r=mpf(rc.numerator)/rc.denominator; E=sqrt((r-2)**2/(r*(r-3)))
    c=coeffs(E)
    print(rc, nstr(c[1],25), Fraction(str(nstr(c[1],30))).limit_denominator(10**7), '|', nstr(c[2],22), Fraction(str(nstr(c[2],22))).limit_denominator(10**9), flush=True)
