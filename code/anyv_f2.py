from mpmath import mp, mpf, pi, sqrt, nstr, matrix, lu_solve, identify
from fractions import Fraction
import anyv_area as AA
mp.dps=40
def fcoef(E,th):
    hs=[mpf('0.004')*k for k in (1,2,3,4)]
    ys=[AA.area(E,h,th) for h in hs]
    M=matrix([[1,h**2,h**4,h**6] for h in hs]); c=lu_solve(M,matrix(ys))
    return -c[1]/c[0]
if __name__=="__main__":
 for rc in [Fraction(15,4),Fraction(7,2),Fraction(10,3),Fraction(13,4)]:
  r=mpf(rc.numerator)/rc.denominator; E=sqrt((r-2)**2/(r*(r-3)))
  fe=fcoef(E,pi/2); fp=fcoef(E,pi/180*mpf('0.5'))
  print(rc, nstr(fe,14), Fraction(str(nstr(fe,14))).limit_denominator(5000), nstr(fp,14), Fraction(str(nstr(fp,14))).limit_denominator(5000), flush=True)
