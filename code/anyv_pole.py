from mpmath import mp, mpf, sqrt, findroot, matrix, lu_solve, nstr, pi
from fractions import Fraction
import anyv_area as AA
mp.dps=50
def polar_area(E,a):
    E=mpf(E); a=mpf(a); p2=E*E-1
    E2=E*E; disc=(3*E2-4)**2-16*(1-E2)
    rc=[x for x in [(-(3*E2-4)+sg*sqrt(disc))/(2*(1-E2)) for sg in (1,-1)] if 3<x.real<=4][0]
    r0=findroot(lambda r: AA.Lf(r,a,E), rc)
    return pi*(AA.Qf(r0,a,E)+a*a*p2)/p2
def fpole(E):
    hs=[mpf('0.002')*k for k in (1,2,3,4,5)]
    ys=[polar_area(E,h) for h in hs]
    M=matrix([[1,h**2,h**4,h**6,h**8] for h in hs]); c=lu_solve(M,matrix(ys))
    return -c[1]/c[0]
for rc in [Fraction(15,4),Fraction(7,2),Fraction(10,3),Fraction(13,4),Fraction(31,10)]:
    r=mpf(rc.numerator)/rc.denominator; E=sqrt((r-2)**2/(r*(r-3)))
    f=fpole(E); fe=1/(2*r*(6-r))
    print(rc, nstr(f,20), Fraction(str(nstr(f,18))).limit_denominator(100000), 'minus edge:', Fraction(str(nstr(f-fe,18))).limit_denominator(100000))
