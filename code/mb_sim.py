# Brute-force capture region for E=1 (slow) particles, edge-on: grid over (L, Q>=0), capture iff R(r)>0 for all r>r_+
import numpy as np
def capture_area(a, n=1500):
    rp=1+np.sqrt(1-a*a)
    rr=rp+np.geomspace(1e-9,60,3000)
    Ls=np.linspace(-8,8,n); Qs=np.linspace(0,30,n//2)
    dL=Ls[1]-Ls[0]; dQ=Qs[1]-Qs[0]
    area=0.0
    for L in Ls:
        # area element in (x=-L, y=sqrt(Q)) plane: dy = dQ/(2 sqrt Q); count y>0 half and double
        R=(rr[None,:]**2+a*a-a*L)**2-(rr[None,:]**2-2*rr[None,:]+a*a)*(rr[None,:]**2+(L-a)**2+Qs[:,None])
        cap=np.all(R>0,axis=1)
        y=np.sqrt(Qs)
        # integrate indicator over y using Q-grid -> convert: capture set is y<y_max(L); take max y captured
        if cap.any():
            ymax=y[cap].max()
            area+=2*ymax*dL
    return area
if __name__=="__main__":
    from mpmath import pi, sqrt, ellipe, mpf
    for a in [0.5,0.9]:
        A=7*pi+pi*sqrt(1-a*a)+16*sqrt(1+a)*ellipe(2*a/(1+a))
        print(a, capture_area(a), float(A), flush=True)
