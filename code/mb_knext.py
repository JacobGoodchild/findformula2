# Slow-particle edge-on capture area, extremal Kerr-Newman (a^2+Q^2=1)
import numpy as np
from mpmath import mp, mpf, sqrt, quad, pi, polyroots, cos, sin
mp.dps=30
def A_formula(a):
    a=mpf(a)
    B=lambda y: a*a*(2*y+1)-(y*y-y-1)**2
    co=[-1,2,1,2*a*a-2,a*a-1]   # B expanded: -(y^4-2y^3-y^2+2y+1)+2a^2 y+a^2
    rts=sorted([x.real for x in polyroots(co,maxsteps=500,extraprec=300) if abs(x.imag)<1e-20])
    hi=rts[-1]; lo=max(mpf(1),[x for x in rts if x<hi][-1])
    f=lambda t:(lambda y: y*(3*y+1)*(y-1)*sqrt(max(B(y),0))*(hi-lo)*sin(t)/2)(lo+(hi-lo)*(1-cos(t))/2)
    return 2/a**2*quad(f,[0,pi/2,pi]), lo, hi
def A_sim(a,n=3000):
    q=1-a*a; rr=1+np.geomspace(1e-10,80,6000)
    def cap(L,Q):
        R=(rr**2+a*a-a*L)**2-(rr-1)**2*(rr**2+(L-a)**2+Q)
        return np.all(R>0)
    Ls=np.linspace(-8,8,n+1); Ls=(Ls[:-1]+Ls[1:])/2; dL=Ls[1]-Ls[0]; area=0
    for L in Ls:
        if not cap(L,0.0): continue
        lo,hi=0.0,40.0
        for _ in range(50):
            m=(lo+hi)/2
            if cap(L,m): lo=m
            else: hi=m
        area+=2*np.sqrt(lo)*dL
    return area
if __name__=="__main__":
    for a in [1.0,0.8,0.5,0.3]:
        A,lo,hi=A_formula(a); print(a, A, float(lo), float(hi), A_sim(a), flush=True)
    print('7pi+16sqrt2',7*pi+16*sqrt(2))
