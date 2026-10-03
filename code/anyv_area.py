# Capture area (in impact-parameter plane, times p^2 = E^2-1) for Kerr, energy E per unit mass, incidence angle theta
import sympy as sp
from mpmath import mp, mpf, sqrt, quad, pi, cos, sin, findroot
r,a,E=sp.symbols('r a E',positive=True)
exec(open('anyv.py').read().split("t,p=sp.symbols")[0].replace("timeout",""))  # re-solve
mp.dps=40
Lf=sp.lambdify((r,a,E),sol[1][L],'mpmath'); Qf=sp.lambdify((r,a,E),sol[1][Q],'mpmath')
def area(Ev,av,th):
    Ev=mpf(Ev); av=mpf(av); p2=Ev*Ev-1; c2=cos(th)**2; s2=sin(th)**2
    g=lambda rr: Qf(rr,av,Ev)+av*av*p2*c2-Lf(rr,av,Ev)**2*c2/s2
    # find range of r where g>=0 : scan
    import numpy as np
    from mpmath import sqrt as msq
    E2=Ev*Ev
    if abs(E2-1)<mpf(10)**-30: rc=mpf(4)
    else:
        disc=(3*E2-4)**2-16*(1-E2); cands=[(-(3*E2-4)+sg*msq(disc))/(2*(1-E2)) for sg in (1,-1)]
        rc=[x for x in cands if 3<x.real<=4.0001][0]
    width=6*av+mpf(10)**-6
    rs=[rc-width+2*width*k/4000 for k in range(4001)]
    vals=[g(x) for x in rs]
    idx=[i for i,v in enumerate(vals) if v.real>0 and abs(v.imag)<1e-20]
    lo=findroot(g,(rs[idx[0]-1],rs[idx[0]]),solver='bisect'); hi=findroot(g,(rs[idx[-1]],rs[idx[-1]+1]),solver='bisect')
    from mpmath import diff
    f=lambda t:(lambda x: sqrt(max(g(x).real,0))*abs(diff(lambda z: Lf(z,av,Ev),x))*(hi-lo)*sin(t)/2)(lo+(hi-lo)*(1-cos(t))/2)
    return 2/sin(th)*quad(f,[0,pi/2,pi])/p2
if __name__=="__main__":
    pass #print(area('1.0000001','0.5',pi/2)*mpf('0.0000002'))
