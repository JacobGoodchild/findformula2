from mpmath import *
import mpmath
import master
from series_incl import coeffs
mp.dps=80
from fractions import Fraction
def rat(v,maxden=10**9):
    f=Fraction(str(nstr(v,40))).limit_denominator(maxden); return f, abs(v-mpf(f.numerator)/f.denominator)
C2s=[mpf(k)/8+mpf(1)/16 for k in range(0,8)]
table=[coeffs(acos(sqrt(C2)),N=16,h=mpf(1)/10**4) for C2 in C2s]
for j in range(4):
    ys=[t[j] for t in table]
    deg=j+2
    M=mpmath.matrix([[C2**i for i in range(deg+1)] for C2 in C2s[:deg+1]])
    p=mpmath.lu_solve(M,mpmath.matrix(ys[:deg+1]))
    resid=max(abs(sum(p[i]*C2**i for i in range(deg+1))-y) for C2,y in zip(C2s,ys))
    print('a^%d:'%(2*j+2), [rat(p[i])[0] for i in range(deg+1)], 'resid',nstr(resid,3), 'raterr',nstr(max(rat(p[i])[1] for i in range(deg+1)),3))
