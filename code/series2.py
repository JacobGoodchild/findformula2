from mpmath import mp, mpf, pi, identify, sqrt
import mpmath

from area_closed import A_closed
mp.dps=250
f=lambda x: A_closed(sqrt(x))/(27*pi)
N=30
xs=[mpf(1)/10**6*(k+1) for k in range(N)]
ys=[f(x) for x in xs]
h=mpf(1)/10**6
M=mpmath.matrix([[(x/h)**j for j in range(N)] for x in xs]); c=mpmath.lu_solve(M,mpmath.matrix(ys)); c=[c[j]/h**j for j in range(N)]
for j in range(7):
    v=c[j]; print(j, mpmath.nstr(v,40))
    for den in [2**p*3**q for p in range(12) for q in range(14)]:
        t=v*den
        if abs(t-mpmath.nint(t))<mpf(10)**-25: print('   =',mpmath.nint(t),'/',den); break
