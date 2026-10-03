from mpmath import mp, mpf, sqrt, quad, pi, cos, sin, polyroots
import numpy as np
import master
mp.dps=30
def P_int(a,th):
    a=mpf(a); c=cos(th); u=a*a*c*c
    co=master.Npoly(a,c); N=lambda r: sum(co[i]*r**(6-i) for i in range(7))
    rts=sorted([x.real for x in polyroots(co,maxsteps=800,extraprec=400) if abs(x.imag)<mpf(10)**-18 and x.real>1])
    r1,r2=rts[0],rts[-1]
    g=lambda r: 4*((r-1)**3+1-a*a)*sqrt(r*(r*r-u))/((r-1)**2*sqrt(N(r)))
    f=lambda t: g(r1+(r2-r1)*(1-cos(t))/2)*(r2-r1)*sin(t)/2
    return 2*quad(f,[0,pi/2,pi]).real
def P_poly(a,th,N=2000000):
    _,(r1,r2)=master.A_master(mpf(a),mpf(th)); r1,r2=float(r1),float(r2); c=np.cos(th); s=np.sin(th)
    t=np.linspace(0,np.pi,N); r=r1+(r2-r1)*(1-np.cos(t))/2
    xi=(r**2*(3-r)-a*a*(r+1))/(a*(r-1)); eta=r**3*(4*a*a-r*(r-3)**2)/(a*a*(r-1)**2)
    X=-xi/s; Y=np.sqrt(np.clip(eta+a*a*c*c-xi**2*c*c/(s*s),0,None))
    return 2*np.sum(np.hypot(np.diff(X),np.diff(Y)))
for a,deg in [(0.5,60),(0.9,17),(0.99,45)]:
    th=deg*np.pi/180
    print(a,deg, P_int(a,mpf(deg)*pi/180), P_poly(a,th))
