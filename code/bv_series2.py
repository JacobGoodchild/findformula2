from mpmath import mp, mpf, pi, sqrt, matrix, lu_solve, nstr, pslq, asin, atan
from bv_equator import Omega_eq
mp.dps=60
A=16*pi+15*sqrt(3)
h=mpf(1)/10**5; N=14
zs=[h*(k+1) for k in range(N)]
ys=[Omega_eq(z)[0]/z**2 for z in zs]
M=matrix([[(z/h)**j for j in range(N)] for z in zs]); c=lu_solve(M,matrix(ys)); c=[c[j]/h**j for j in range(N)]
for j in range(6):
    print(j, nstr(c[j],30), pslq([c[j],pi,sqrt(3)],tol=mpf(10)**-22,maxcoeff=10**7,maxsteps=10**6))
# check I2 closed form
z=mpf('0.3'); k=z/(1+z)
from mpmath import quad, sin
I2n=2*quad(lambda p: sqrt(sin(p)**2-k*k)/sin(p),[pi/6,pi/2])
I2c=2*asin(sqrt(3)/(2*sqrt(1-k*k)))-2*k*atan(sqrt(3)*k/sqrt(1-4*k*k))
print('I2', I2n, I2c)
