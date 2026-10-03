import numpy as np
from mpmath import mp, mpf, sqrt, findroot
from mb_avg_sim import Lc
mp.dps=30
def Lc_param(a,mu):
    a=mpf(a); mu=mpf(mu)
    m=lambda y: -(y**4-2*y**3+a*a)/(a*sqrt(3*y**4-4*y**3+a*a))
    y1=1+sqrt(1-a); y2=1+sqrt(1+a)
    # bisection on y for m(y)=mu (m decreasing from +1 to -1)
    lo,hi=y1,y2
    for _ in range(100):
        mid=(lo+hi)/2
        if m(mid)>mu: lo=mid
        else: hi=mid
    y=(lo+hi)/2
    return sqrt((3*y**4-4*y**3+a*a)/(y-1)**2)
for a in [0.5,0.9,0.99]:
    rp=1+np.sqrt(1-a*a); rr=rp+np.geomspace(1e-10,80,6000)
    for mu in [-0.9,0.0,0.5,0.95]:
        print(a,mu, Lc(a,mu,rr), float(Lc_param(a,mu)))
