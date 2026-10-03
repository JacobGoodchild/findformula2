from mpmath import mp, mpf, pi, taylor, identify, sqrt, pslq
from area_closed import A_closed
mp.dps=60
# A as function of x=a^2; fit series by sampling small x via Richardson/polyfit
import mpmath
f=lambda x: A_closed(sqrt(x))/(27*pi)
xs=[mpf(1)/10**4*(k+1) for k in range(14)]
ys=[f(x) for x in xs]
# polynomial interpolation
N=len(xs)
M=mpmath.matrix([[x**j for j in range(N)] for x in xs]); c=mpmath.lu_solve(M,mpmath.matrix(ys))
for j in range(6): print(j, mpmath.nstr(c[j],25), identify(c[j]), identify(c[j],['pi','sqrt(3)']))
