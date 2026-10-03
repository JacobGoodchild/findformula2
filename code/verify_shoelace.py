# Independent check: shoelace area of the critical curve polygon (equatorial observer)
import numpy as np
from mpmath import mp, mpf
from area_closed import A_closed
def shoelace(a, N=400000):
    c=np.roots([1,-6,9,-4*a*a]).real; c.sort(); r1,r2=c[1],c[2]
    t=np.linspace(0,np.pi,N); r=r1+(r2-r1)*(1-np.cos(t))/2
    xi=(r**2*(3-r)-a*a*(r+1))/(a*(r-1)); eta=r**3*(4*a*a-r*(r-3)**2)/(a*a*(r-1)**2)
    al=-xi; be=np.sqrt(np.clip(eta,0,None))
    X=np.concatenate([al,al[::-1]]); Y=np.concatenate([be,-be[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
for a in [0.2,0.5,0.8,0.95,0.999]:
    print(a, shoelace(a), A_closed(mpf(a)))
