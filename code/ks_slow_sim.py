import numpy as np
from mpmath import mpf, sqrt, pi, ellipe
def A_formula(a,b):
    B=1-b; return float(16*sqrt(B+a)*ellipe(2*a/(B+a))+pi*(3*B+4+sqrt(B*B-a*a)))
def A_sim(a,b,n=3000):
    rp=(1-b)+np.sqrt((1-b)**2-a*a)
    rr=rp+np.geomspace(1e-10,80,6000); rho=rr*(rr+2*b); D=rho-2*rr+a*a
    def cap(L,Q):
        R=(rho+a*a-a*L)**2-D*(rho+(L-a)**2+Q)
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
for a,b in [(0.5,0.2),(0.3,0.5),(0.6,0.4)]:
    print(a,b, A_sim(a,b), A_formula(a,b), flush=True)
