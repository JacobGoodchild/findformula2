import numpy as np
from mb_sim2 import captured
from mpmath import mpf, sqrt, pi, ellipk, ellipe
def M1_formula(a):
    a=mpf(a); m=2*a/(1+a); s=sqrt(1+a)
    return float(-(32*s*((a*a+2)*ellipe(m)-2*(1-a)*ellipk(m))+35*pi*a*a+10*pi*(1-sqrt(1-a*a)))/(5*a))
def moment_sim(a,n=3000):
    rp=1+np.sqrt(1-a*a); rr=rp+np.geomspace(1e-10,80,6000)
    Ls=np.linspace(-8,8,n+1); Ls=(Ls[:-1]+Ls[1:])/2; dL=Ls[1]-Ls[0]
    m=0.0
    for L in Ls:
        if not captured(a,L,0.0,rr): continue
        lo,hi=0.0,40.0
        for _ in range(45):
            q=(lo+hi)/2
            if captured(a,L,q,rr): lo=q
            else: hi=q
        m+=L*2*np.sqrt(lo)*dL
    return m
for a in [0.5,0.9]:
    print(a, moment_sim(a), M1_formula(a), flush=True)
