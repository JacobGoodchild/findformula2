# mean captured L_z per unit energy, any speed, incidence theta: odd series in a
import sympy as sp
from mpmath import mp, mpf, sqrt, quad, pi, cos, sin, findroot, matrix, lu_solve, nstr, acos
from fractions import Fraction
import anyv_area2 as AA
mp.dps=50
def mean_Lz_over_E(Ev,av,th):
    Ev=mpf(Ev); av=mpf(av); p2=Ev*Ev-1; c2=cos(th)**2; s2=sin(th)**2
    g=lambda rr: AA.Qf(rr,av,Ev)+av*av*p2*c2-AA.Lf(rr,av,Ev)**2*c2/s2
    E2=Ev*Ev; disc=(3*E2-4)**2-16*(1-E2)
    rc=[x for x in [(-(3*E2-4)+sg*sqrt(disc))/(2*(1-E2)) for sg in (1,-1)] if 3<x.real<=4][0]
    w=6*av; rs=[rc-w+2*w*k/2000 for k in range(2001)]
    idx=[i for i,x in enumerate(rs) if g(x).real>0]
    lo=findroot(g,(rs[idx[0]-1],rs[idx[0]]),solver='anderson'); hi=findroot(g,(rs[idx[-1]],rs[idx[-1]+1]),solver='anderson')
    def integ(wt):
        f=lambda t:(lambda x: wt(x)*sqrt(max(g(x).real,0))*abs(AA.dLf(x,av,Ev))*(hi-lo)*sin(t)/2)(lo+(hi-lo)*(1-cos(t))/2)
        return quad(f,[0,pi/2,pi])
    return integ(lambda x: AA.Lf(x,av,Ev))/integ(lambda x:1)/Ev
def odd_coeffs(E,th):
    hs=[mpf('0.003')*k for k in range(1,7)]
    ys=[mean_Lz_over_E(E,h,th)/h for h in hs]
    M=matrix([[h**(2*j) for j in range(6)] for h in hs]); c=lu_solve(M,matrix(ys))
    return c[0],c[1]
if __name__=="__main__":
    for rc in [Fraction(15,4),Fraction(7,2),Fraction(10,3),Fraction(13,4),Fraction(18,5),Fraction(16,5),Fraction(31,10)]:
        r=mpf(rc.numerator)/rc.denominator; E=sqrt((r-2)**2/(r*(r-3)))
        c1,c3=odd_coeffs(E,acos(sqrt(mpf(1)/2)))
        print(rc, nstr(c1,20), nstr(c3,20), Fraction(str(nstr(c3,18))).limit_denominator(10**9), flush=True)
