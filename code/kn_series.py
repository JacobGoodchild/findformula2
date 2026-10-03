from mpmath import mp, mpf, pi, acos, sqrt, matrix, lu_solve, nstr
from fractions import Fraction
import kn_master_check as K
mp.dps=60
def rat(v):
    f=Fraction(str(nstr(v,45))).limit_denominator(10**7); return f, abs(v-mpf(f.numerator)/f.denominator)
# fit A/(27pi)-1 = sum_{i+j<=3, (i,j)!=(0,0)} c_ij x^i q^j, x=a^2, at fixed C=cos^2
def fit(C):
    th=acos(sqrt(mpf(C)))
    h=mpf(1)/10**4
    pts=[(h*i,h*j) for i in range(0,5) for j in range(0,5) if 0<i+j<=4]
    mons=[(i,j) for i in range(0,5) for j in range(0,5) if 0<i+j<=4]
    rows=[];ys=[]
    for (x,q) in pts:
        if x==0: x=h/10**3   # avoid a=0
        rows.append([(x/h)**i*(q/h)**j for (i,j) in mons]); ys.append(K.A_kn(sqrt(x),th,q)[0]/(27*pi)-1)
    c=lu_solve(matrix(rows),matrix(ys))
    return {m:c[k]/h**(m[0]+m[1]) for k,m in enumerate(mons)}
for C in ['0.25','0.75']:
    d=fit(C)
    print(C,{m:(rat(v)[0],nstr(rat(v)[1],2)) for m,v in d.items() if m[0]+m[1]<=2})
