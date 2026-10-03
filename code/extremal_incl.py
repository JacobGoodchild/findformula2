# Extremal Kerr (a=1) shadow area for inclination th (observer at infinity)
from mpmath import mp, mpf, sqrt, quad, cos, sin, pi, polyroots, acos, identify, pslq, asin
import numpy as np
mp.dps=30
def Q(r,c2):  # beta^2 * sin^0 : eta + c^2 - xi^2 cot^2
    s2=1-c2; return r**3*(4-r)+c2-(c2/s2)*(r*r-2*r-1)**2
def Apoly(c2):
    s2=1-c2
    # coefficients of Q in r (degree 4)
    k=c2/s2
    # r^3(4-r) = -r^4+4r^3 ; (r^2-2r-1)^2 = r^4-4r^3+2r^2+4r+1
    return [-(1+k), 4*(1+k), -2*k, -4*k, c2-k]
def area(th):
    c2=cos(th)**2; s=sin(th)
    rts=sorted([x.real for x in polyroots(Apoly(c2),maxsteps=200,extraprec=100) if abs(x.imag)<mpf(10)**-20])
    hi=max(rts); lo=[x for x in rts if x<hi]
    lo=max(lo) if lo else mpf(1)
    lo=max(lo,mpf(1))
    return 4/s*quad(lambda r:(r-1)*sqrt(max(Q(r,c2),0)),[lo,hi]), lo, hi, rts
def shoelace(th,N=200000):
    th=float(th); c2=np.cos(th)**2; s=np.sin(th)
    r=np.linspace(1,4,N); b2=r**3*(4-r)+c2-c2/(1-c2)*(r*r-2*r-1)**2
    m=b2>=0; r=r[m]; al=(r*r-2*r-1)/s; be=np.sqrt(b2[m])
    X=np.concatenate([al,al[::-1]]); Y=np.concatenate([be,-be[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
if __name__=="__main__":
    for deg in [1,10,30,45,47,50,60,75,89.999]:
        th=mpf(deg)*pi/180; A,lo,hi,rts=area(th)
        print(deg, A, shoelace(th), lo, hi)
