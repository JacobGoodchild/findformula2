from mpmath import mp, mpf, sqrt, ellipk, ellipe, ellippi, polyroots, pi, quad
import numpy as np
mp.dps=30
def P_closed(a):
    a=mpf(a)
    r0,rp,rr=sorted([x.real for x in polyroots([1,-6,9,-4*a*a],maxsteps=300,extraprec=300)])
    m=(rr-rp)/(rr-r0); n=(rr-rp)/(rr-1)
    K,E,Pi_=ellipk(m),ellipe(m),ellippi(n,m)
    V=(n*E+(m-n)*K+(2*n*m+2*n-n*n-3*m)*Pi_)/(2*(n-1)*(m-n))
    return 16/sqrt(rr-r0)*(r0*K+(rr-r0)*E-K+(1-a*a)*V/(rr-1)**2)
def P_poly(a,N=2000000):
    c=np.sort(np.roots([1,-6,9,-4*a*a]).real); rp,rr=c[1],c[2]
    t=np.linspace(0,np.pi,N); r=rp+(rr-rp)*(1-np.cos(t))/2
    xi=(r**2*(3-r)-a*a*(r+1))/(a*(r-1)); eta=r**3*(4*a*a-r*(r-3)**2)/(a*a*(r-1)**2)
    X=-xi; Y=np.sqrt(np.clip(eta,0,None))
    return 2*np.sum(np.hypot(np.diff(X),np.diff(Y)))
for a in [0.001,0.3,0.7,0.9,0.999]:
    print(a, P_closed(a), P_poly(a))
print('Schw', 2*pi*sqrt(27), ' extremal 18sqrt3', 18*sqrt(3))
