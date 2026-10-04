# Independent check of the universal a^2 shadow formula: exact critical curve of Kerr-like metric
# with Delta = r^2 f(r) + a^2 (Tsukamoto 2018 form), area = 2 Int beta |d alpha|, small a extrapolation.
from mpmath import mp, mpf, sqrt, quad, pi, findroot, diff, sin, cos, matrix, lu_solve, nstr
mp.dps=40
def make(fun):
    def area(a,th):
        D=lambda r: r*r*fun(r)+a*a
        Dp=lambda r: diff(D,r)
        xi=lambda r:((r*r+a*a)*Dp(r)-4*r*D(r))/(a*Dp(r))
        eta=lambda r: r*r*(16*a*a*D(r)-(4*D(r)-r*Dp(r))**2)/(a*a*Dp(r)**2)
        s,c=sin(th),cos(th)
        b2=lambda r: eta(r)+a*a*c*c-xi(r)**2*c*c/(s*s)
        rph=findroot(lambda r: 2*fun(r)-r*diff(fun,r), 3)
        # roots of b2 near rph
        r1=findroot(b2, rph-0.9*a*2); r2=findroot(b2, rph+0.9*a*2)
        if r1>r2: r1,r2=r2,r1
        g=lambda t:(lambda r: sqrt(max(b2(r),0))*abs(diff(xi,r))/s*(r2-r1)*sin(t)/2)(r1+(r2-r1)*(1-cos(t))/2)
        return 2*quad(g,[0,pi/2,pi]), rph
    return area
def Kform(fun,r,C):
    f0=fun(r); f2=diff(fun,r,2)
    return (C*(8*f0**2-4*f0*f2*r**2+2*f0+3*f2*r**2)-6*f0-f2*r**2)/(2*f0*r**2*(2*f0-f2*r**2))
tests={'Bardeen g=0.4':lambda r: 1-2*r**2/(r**2+mpf('0.16'))**1.5,
       'Hayward l=0.5':lambda r: 1-2*r**2/(r**3+2*mpf('0.25')),
       'Kerr':lambda r: 1-2/r}
if __name__=="__main__":
  for name,fun in tests.items():
      ar=make(fun)
      for thdeg in [90,50]:
          th=thdeg*pi/180; C=cos(th)**2
          hs=[mpf('0.01'),mpf('0.02'),mpf('0.03')]; ys=[]
          for h in hs:
              A,rph=ar(h,th); A0=pi*rph**2/fun(rph); ys.append((A/A0-1)/h**2)
          K=lu_solve(matrix([[1,h**2,h**4] for h in hs]),matrix(ys))[0]
          print(name,thdeg,'r_ph=',nstr(rph,10),'numeric',nstr(K,12),'formula',nstr(Kform(fun,rph,C),12),flush=True)
