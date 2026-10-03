import numpy as np
from mpmath import mp, mpf, sqrt, quad, pi, cos, sin
mp.dps=30
def A_formula(a,th):
    a=mpf(a); s=sin(th); c2=cos(th)**2
    N4=lambda r: -((r-a*(1+s))**2-2*a*(1+s))*((r-a*(1-s))**2-2*a*(1-s))
    rD=a*(1+s)+sqrt(2*a*(1+s)); rC=a*(1-s)+sqrt(2*a*(1-s))
    lo=max(rC,a)
    f=lambda t:(lambda r:(3*r*r+2*(1-a)*r-a*a*c2)/sqrt(N4(r))*(rD-lo)*sin(t)/2)(lo+(rD-lo)*(1-cos(t))/2)
    I=4*quad(f,[0,pi/2,pi]).real
    B=(a*s*s+2)*sqrt(N4(a))/(a*s**2) if rC<a else 0
    return I+B
def shoelace(a,th,N=1500000):
    b=1-a; s=np.sin(th); c=np.cos(th)
    rD=a*(1+s)+np.sqrt(2*a*(1+s)); rC=a*(1-s)+np.sqrt(2*a*(1-s)); lo=max(rC,a+1e-6)
    t=np.linspace(0,np.pi,N); r=lo+(rD-lo)*(1-np.cos(t))/2
    xi=(b*b+2*b*r+r*r-2*r-1)/(b-1); eta=r*r*(4*(1-b)-(r+2*b-2)**2)/a**2
    be=np.sqrt(np.clip(eta+a*a*c*c-xi**2*c*c/(s*s),0,None)); al=-xi/s
    X=np.concatenate([al,al[::-1]]); Y=np.concatenate([be,-be[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
for a,deg in [(0.5,60),(0.5,85),(0.8,30),(0.3,70),(1.0,60)]:
    th=mpf(deg)*pi/180
    print(a,deg, A_formula(a,th), shoelace(a,float(th)))
