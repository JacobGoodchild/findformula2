import numpy as np
def captured(a,L,Q,rr):
    R=(rr**2+a*a-a*L)**2-(rr**2-2*rr+a*a)*(rr**2+(L-a)**2+Q)
    return np.all(R>0)
def capture_area(a,n=4000):
    rp=1+np.sqrt(1-a*a); rr=rp+np.geomspace(1e-10,80,6000)
    Ls=np.linspace(-8,8,n+1); Ls=(Ls[:-1]+Ls[1:])/2; dL=Ls[1]-Ls[0]
    area=0.0
    for L in Ls:
        if not captured(a,L,0.0,rr): continue
        lo,hi=0.0,40.0
        for _ in range(50):
            m=(lo+hi)/2
            if captured(a,L,m,rr): lo=m
            else: hi=m
        area+=2*np.sqrt(lo)*dL
    return area
if __name__=="__main__":
    from mpmath import pi, sqrt, ellipe
    for a in [0.5,0.9,0.999]:
        A=7*pi+pi*sqrt(1-a*a)+16*sqrt(1+a)*ellipe(2*a/(1+a))
        print(a, capture_area(a), float(A), flush=True)
