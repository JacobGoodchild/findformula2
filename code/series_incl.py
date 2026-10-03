from mpmath import *
import mpmath
import master
mp.dps=60
def coeffs(th,N=10,h=mpf(1)/10**3):
    xs=[h*(k+1) for k in range(N)]
    ys=[master.A_master(sqrt(x),th)[0].real/(27*pi)-1 for x in xs]
    M=mpmath.matrix([[(x/h)**(j+1) for j in range(N)] for x in xs]); c=mpmath.lu_solve(M,mpmath.matrix(ys))
    return [c[j]/h**(j+1) for j in range(N)]
if __name__=="__main__":
 for C2 in ['0','0.25','0.5','0.75']:
  th=acos(sqrt(mpf(C2)))
  c=coeffs(th)
  print(C2, [nstr(v,20) for v in c[:3]])
