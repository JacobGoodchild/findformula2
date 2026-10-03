from mpmath import mp, mpf, sqrt, pi, cos, sin, quad, polyroots, nstr
import ext_closed, master
mp.dps=50
def A_master_hp(a,th):
    c=cos(th); u=a*a*c*c
    co=master.Npoly(a,c)
    N=lambda r: sum(co[i]*r**(6-i) for i in range(7))
    rts=sorted([x.real for x in polyroots(co,maxsteps=500,extraprec=500) if abs(x.imag)<mpf(10)**-30 and x.real>1])
    r1,r2=rts[0],rts[-1]
    P=lambda r: (2*(u*u-5*a*a*u+6*u+3*a*a-3) -2*(2*u*u+a*a*u-3*u-3*a*a+3)*r +2*(u*u-6*u-3)*r**2
                 +4*(u+3)*r**3 +2*(u-3)*r**4 -2*(a*a-1)*(u*u+6*u-3)/(r-1))/(u-1)
    from mpmath import cos as mcos, sin as msin
    def f(t):
        r=r1+(r2-r1)*(1-mcos(t))/2; n=N(r)
        if n<=0: return mpf(0)
        return P(r)/sqrt(n)*(r2-r1)*msin(t)/2
    return quad(f,[0,pi/1000,pi/2,pi]).real
if __name__=="__main__":
 for deg in [30,60,75,90]:
  th=mpf(deg)*pi/180
  A1=ext_closed.A_closed(th) if deg!=90 else 16*pi+15*sqrt(3)
  out=[]
  for e in [mpf(10)**-8,mpf(10)**-10,mpf(10)**-12]:
      out.append(nstr((A_master_hp(1-e,th)-A1)/sqrt(e),12))
  print(deg, out)
