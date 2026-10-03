# Master formula: Kerr shadow area for any spin a (0<a<1) and inclination th, observer at infinity
from mpmath import mp, mpf, sqrt, quad, cos, sin, pi, polyroots
import numpy as np
mp.dps=30
def Npoly(a,c):
    u=a*a*c*c
    # N = -r^6+6r^5-(9+2u)r^4+4a^2 r^3+6u r^2-4a^2 u r - u^2 (r-1)^2
    return [-1,6,-(9+2*u),4*a*a,6*u-u*u,-4*a*a*u+2*u*u,-u*u]
def A_master(a,th):
    c=cos(th); u=a*a*c*c
    co=Npoly(a,c)
    N=lambda r: sum(co[i]*r**(6-i) for i in range(7))
    rts=sorted([x.real for x in polyroots(co,maxsteps=2000,extraprec=300) if abs(x.imag)<mpf(10)**-15 and x.real>1])
    r1,r2=rts[0],rts[-1]
    P=lambda r: (2*(u*u-5*a*a*u+6*u+3*a*a-3) -2*(2*u*u+a*a*u-3*u-3*a*a+3)*r +2*(u*u-6*u-3)*r**2
                 +4*(u+3)*r**3 +2*(u-3)*r**4 -2*(a*a-1)*(u*u+6*u-3)/(r-1))/(u-1)
    return quad(lambda r:P(r)/sqrt(N(r)),[r1,(r1+r2)/2,r2]), (r1,r2)
def A_direct(a,th):
    c=cos(th); s=sin(th)
    xi=lambda r:(r**2*(3-r)-a*a*(r+1))/(a*(r-1)); eta=lambda r:r**3*(4*a*a-r*(r-3)**2)/(a*a*(r-1)**2)
    b2=lambda r: eta(r)+a*a*c*c-xi(r)**2*c*c/(s*s)
    dxi=lambda r: -2*((r-1)**3+1-a*a)/(a*(r-1)**2)
    _,(r1,r2)=A_master(a,th)
    return 2/s*quad(lambda r: sqrt(b2(r))*abs(dxi(r)),[r1,(r1+r2)/2,r2])
def shoelace(a,th,n=400000):
    c=np.cos(th); s=np.sin(th); a=float(a)
    _,(r1,r2)=A_master(mpf(a),mpf(th)); r1,r2=float(r1),float(r2)
    t=np.linspace(0,np.pi,n); r=r1+(r2-r1)*(1-np.cos(t))/2
    xi=(r**2*(3-r)-a*a*(r+1))/(a*(r-1)); eta=r**3*(4*a*a-r*(r-3)**2)/(a*a*(r-1)**2)
    al=-xi/s; be=np.sqrt(np.clip(eta+a*a*c*c-xi**2*c*c/(s*s),0,None))
    X=np.concatenate([al,al[::-1]]); Y=np.concatenate([be,-be[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
if __name__=="__main__":
    for a,deg in [(0.5,17),(0.9,17),(0.5,60),(0.94,30),(0.99,85),(0.3,89),(0.998,45)]:
        th=mpf(deg)*pi/180; a=mpf(a)
        print(a,deg, A_master(a,th)[0], A_direct(a,th), shoelace(a,float(th)))
