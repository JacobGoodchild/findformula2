from mpmath import mp, mpf, pi, acos, sqrt, matrix, lu_solve, nstr
from fractions import Fraction
import perimeter_incl as PI
mp.dps=60
def rat(v):
    f=Fraction(str(nstr(v,45))).limit_denominator(10**8); return f, abs(v-mpf(f.numerator)/f.denominator)
def coeffs(C,N=10,h=mpf(1)/10**4):
    th=acos(sqrt(mpf(C)))
    xs=[h*(k+1) for k in range(N)]
    ys=[PI.P_int(sqrt(x),th)/(2*pi*sqrt(27))-1 for x in xs]
    M=matrix([[(x/h)**(j+1) for j in range(N)] for x in xs]); c=lu_solve(M,matrix(ys))
    return [c[j]/h**(j+1) for j in range(N)]
Cs=[mpf(k)/8+mpf(1)/16 for k in range(5)]
tab=[coeffs(C) for C in Cs]
for j in range(2):
    ys=[t[j] for t in tab]; deg=j+1
    M=matrix([[C**i for i in range(deg+1)] for C in Cs[:deg+1]]); p=lu_solve(M,matrix(ys[:deg+1]))
    res=max(abs(sum(p[i]*C**i for i in range(deg+1))-y) for C,y in zip(Cs,ys))
    print('a^%d'%(2*j+2),[str(rat(p[i])[0]) for i in range(deg+1)],'resid',nstr(res,3),'raterr',nstr(max(rat(p[i])[1] for i in range(deg+1)),3), flush=True)
