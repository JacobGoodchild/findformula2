from mpmath import mp, mpf, pi, sqrt, nstr, matrix, lu_solve
from fractions import Fraction
import anyv_area as AA
from anyv_pole import polar_area
mp.dps=60
def coefs(fun):
    hs=[mpf('0.003')*k for k in range(1,8)]
    ys=[fun(h) for h in hs]
    M=matrix([[h**(2*j) for j in range(7)] for h in hs]); c=lu_solve(M,matrix(ys))
    return [c[j]/c[0] for j in range(4)]
for rc in [Fraction(15,4),Fraction(7,2),Fraction(10,3),Fraction(13,4),Fraction(31,10),Fraction(18,5),Fraction(16,5)]:
    r=mpf(rc.numerator)/rc.denominator; E=sqrt((r-2)**2/(r*(r-3)))
    ce=coefs(lambda h: AA.area(E,h,pi/2)); cp=coefs(lambda h: polar_area(E,h))
    print(rc, 'edge a2,a4:', nstr(ce[1],20), nstr(ce[2],16), ' pole:', nstr(cp[1],20), nstr(cp[2],16), flush=True)
