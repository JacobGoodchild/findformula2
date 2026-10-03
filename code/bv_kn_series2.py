from mpmath import mp, mpf, pi, matrix, lu_solve, nstr
from fractions import Fraction
from bv_kn_series import coeffs
mp.dps=50
As=[mpf(k)/20 for k in range(1,10)]   # 0.05..0.45
tab=[coeffs(a,N=12) for a in As]
def rat(v): 
    f=Fraction(str(nstr(v,40))).limit_denominator(10**6); return f, abs(v-mpf(f.numerator)/f.denominator)
for j in range(5):
    ys=[t[j]/pi for t in tab]
    deg=min(j+1,len(As)-1)+1
    M=matrix([[a**(2*i) for i in range(deg+1)] for a in As[:deg+1]]); p=lu_solve(M,matrix(ys[:deg+1]))
    res=max(abs(sum(p[i]*a**(2*i) for i in range(deg+1))-y) for a,y in zip(As,ys))
    print(j,[str(rat(p[i])[0]) for i in range(deg+1)],'resid',nstr(res,3))
