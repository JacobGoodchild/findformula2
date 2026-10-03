import numpy as np
def capture_area(v,n=3000):
    E=1/np.sqrt(1-v*v); p2=E*E-1
    rr=1+np.geomspace(1e-9,80,6000)
    def cap(L,Q):
        R=(E*(rr**2+1)-L)**2-(rr-1)**2*(rr**2+(L-E)**2+Q)
        return np.all(R>0)
    Ls=np.linspace(-12*E,12*E,n+1); Ls=(Ls[:-1]+Ls[1:])/2; dL=Ls[1]-Ls[0]
    area=0
    for L in Ls:
        if not cap(L,0.0): continue
        lo,hi=0.0,100*E*E
        for _ in range(50):
            m=(lo+hi)/2
            if cap(L,m): lo=m
            else: hi=m
        area+=2*np.sqrt(lo)*dL
    return area/p2
if __name__=="__main__":
    from anyv_ext_area import area
    for v in [0.6,0.9]:
        print(v, capture_area(v), float(area(v)[0]), flush=True)
