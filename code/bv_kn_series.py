from mpmath import mp, mpf, pi, sqrt, matrix, lu_solve, nstr, pslq, asin
from bv_kn import Omega_kn
mp.dps=50
def coeffs(a,N=12,h=mpf(1)/10**5):
    zs=[h*(k+1) for k in range(N)]
    ys=[Omega_kn(z,a)/z**2 for z in zs]
    M=matrix([[(z/h)**j for j in range(N)] for z in zs]); c=lu_solve(M,matrix(ys)); return [c[j]/h**j for j in range(N)]
if __name__=="__main__":
    for a in [mpf(3)/10, mpf(2)/5]:
        c=coeffs(a)
        for j in range(5):
            print(nstr(a,3), j, nstr(c[j],25), nstr(c[j]/pi,20), pslq([c[j],pi],tol=mpf(10)**-20,maxcoeff=10**8,maxsteps=10**6))
