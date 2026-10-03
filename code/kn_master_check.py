from mpmath import mp, mpf, sqrt, quad, cos, sin, pi, polyroots
import numpy as np
mp.dps=30
def Ncoef(a,c,q):
    u=a*a*c*c
    return [-1,6,-(2*u+4*q+9),4*a*a+12*q,-u*u+6*u-4*a*a*q-4*q*q,2*u*u-4*a*a*u-4*u*q,-u*u]
def A_kn(a,th,q):
    a=mpf(a); q=mpf(q); c=cos(th); u=a*a*c*c; D=u+q-1
    co=Ncoef(a,c,q); N=lambda r: sum(co[i]*r**(6-i) for i in range(7))
    rts=sorted([x.real for x in polyroots(co,maxsteps=800,extraprec=400) if abs(x.imag)<mpf(10)**-18 and x.real>1])
    r1,r2=rts[0],rts[-1]
    P=lambda r:(-2*(-u*u+5*a*a*u+5*u*q-6*u+4*a*a*q-3*a*a+4*q*q-7*q+3)
                -2*(2*u*u+a*a*u+u*q-3*u+2*a*a*q-3*a*a+2*q*q-7*q+3)*r
                -2*(-u*u+6*u+2*q+3)*r**2 -4*(-u+q-3)*r**3 +2*(u+2*q-3)*r**4
                -2*(a*a+q-1)*(u*u+6*u+4*q-3)/(r-1))/D
    f=lambda t:(lambda r: P(r)/sqrt(N(r))*(r2-r1)*sin(t)/2)(r1+(r2-r1)*(1-cos(t))/2)
    return quad(f,[0,pi/2,pi]).real,(r1,r2)
def shoelace(a,th,q,n=600000):
    _,(r1,r2)=A_kn(a,th,q); r1,r2=float(r1),float(r2); c=np.cos(th); s=np.sin(th)
    t=np.linspace(0,np.pi,n); r=r1+(r2-r1)*(1-np.cos(t))/2
    xi=-(r**3-3*r**2+a*a*r+a*a+2*q*r)/(a*(r-1)); eta=r**2*(4*a*a*(r-q)-(r*r-3*r+2*q)**2)/(a*a*(r-1)**2)
    al=-xi/s; be=np.sqrt(np.clip(eta+a*a*c*c-xi**2*c*c/(s*s),0,None))
    X=np.concatenate([al,al[::-1]]); Y=np.concatenate([be,-be[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
if __name__=="__main__":
    for a,deg,q in [(0.5,17,0.3),(0.7,45,0.2),(0.3,70,0.6),(0.9,30,0.1),(0.5,60,0.0)]:
        th=float(deg)*np.pi/180
        print(a,deg,q, A_kn(a,mpf(th),q)[0], shoelace(a,th,q))
