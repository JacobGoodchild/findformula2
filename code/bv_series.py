from mpmath import mp, mpf, pi, sqrt, matrix, lu_solve, nstr, pslq, quad, sin, asin
import escape
mp.dps=50
th=pi/2
A=16*pi+15*sqrt(3)
def Om(z):
    return escape.Omega(z,th)
h=mpf(1)/10**4; N=10
zs=[h*(k+1) for k in range(N)]
ys=[Om(z)/z**2 for z in zs]
M=matrix([[(z/h)**j for j in range(N)] for z in zs]); c=lu_solve(M,matrix(ys)); c=[c[j]/h**j for j in range(N)]
for j in range(5):
    print(j, nstr(c[j],25), pslq([c[j],pi,sqrt(3)],maxcoeff=10**5,maxsteps=10**6), nstr(c[j]/A,15))
print('--- loose tol')
for j in range(4):
    print(j, pslq([c[j],pi,sqrt(3)],tol=mpf(10)**-15,maxcoeff=10**6,maxsteps=10**6), pslq([c[j],pi,sqrt(3),1,sqrt(3)*pi],tol=mpf(10)**-15,maxcoeff=10**5,maxsteps=10**6))
