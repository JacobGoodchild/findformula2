from mpmath import *
import mpmath
from ext_closed import A_closed
mp.dps=200
A0=pi*(12+8*sqrt(2))
N=24; h=mpf(1)/10**4
ss=[h*(k+1) for k in range(N)]
ys=[(A_closed(asin(s))-A0) for s in ss]
M=mpmath.matrix([[ (s/h)**(j+1) for j in range(N)] for s in ss]); c=mpmath.lu_solve(M,mpmath.matrix(ys)); c=[c[j]/h**(j+1) for j in range(N)]
for j in range(6):
    print(j+1, nstr(c[j],30), (identify(c[j],['pi','sqrt(2)*pi']) if abs(c[j])>1e-30 else 0))
for j in [1,3,5,7]:
    print(j+1, pslq([c[j],pi,sqrt(2)*pi],maxcoeff=10**8,maxsteps=10**6, tol=mpf(10)**-20))
