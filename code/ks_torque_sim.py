import numpy as np
from mpmath import mpf, sqrt, pi, ellipe, ellipk
def M1_formula(a,b):
    a=mpf(a); B=1-mpf(b); s=sqrt(B+a); m=2*a/(B+a); K=ellipk(m); E=ellipe(m)
    return float((-32*B**2*E*s+32*B**2*K*s-30*pi*B**2-160*B*E*s-32*B*K*a*s+160*B*K*s+30*pi*B*sqrt(B*B-a*a)-96*E*a*a*s-160*K*a*s-105*pi*a*a)/(15*a))
def M1_sim(a,b,n=3000):
    rp=(1-b)+np.sqrt((1-b)**2-a*a); rr=rp+np.geomspace(1e-10,80,6000); rho=rr*(rr+2*b); D=rho-2*rr+a*a
    def cap(L,Q):
        return np.all((rho+a*a-a*L)**2-D*(rho+(L-a)**2+Q)>0)
    Ls=np.linspace(-8,8,n+1); Ls=(Ls[:-1]+Ls[1:])/2; dL=Ls[1]-Ls[0]; m_=0
    for L in Ls:
        if not cap(L,0.0): continue
        lo,hi=0.0,40.0
        for _ in range(48):
            q=(lo+hi)/2
            if cap(L,q): lo=q
            else: hi=q
        m_+=L*2*np.sqrt(lo)*dL
    return m_
for a,b in [(0.5,0.2),(0.3,0.5)]:
    print(a,b, M1_sim(a,b), M1_formula(a,b), flush=True)
