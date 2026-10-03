# Slow-particle (E=1) capture area for edge-on incidence, Kerr: A = (2/a^2) Int (1+u)^2 sqrt(P) (3u^4+4u^3+1-a^2)/u^3 du
from mpmath import mp, mpf, sqrt, quad, pi
mp.dps=30
def A_mb_edge(a):
    a=mpf(a); P=lambda u: a*a-(u*u-1)**2
    lo,hi=sqrt(1-a),sqrt(1+a)
    from mpmath import cos,sin
    f=lambda t:(lambda u:(1+u)**2*sqrt(max(P(u),0))*abs(3*u**4+4*u**3+1-a*a)/u**3*(hi-lo)*sin(t)/2)(lo+(hi-lo)*(1-cos(t))/2)
    return 2/a**2*quad(f,[0,pi/2,pi])
def A_mb_shoelace(a,N=400000):
    import numpy as np
    y=np.linspace(1+np.sqrt(1-a),1+np.sqrt(1+a),N)
    L=-(y**4-2*y**3+a*a)/(a*(y-1)); Q=y**4*(a*a-(y*y-2*y)**2)/(a*a*(y-1)**2)
    X=np.concatenate([-L,(-L)[::-1]]); Y=np.concatenate([np.sqrt(np.clip(Q,0,None)),-np.sqrt(np.clip(Q,0,None))[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
for a in ['0.001','0.5','0.9','0.999']:
    print(a, A_mb_edge(a), A_mb_shoelace(float(a)), 16*pi)
