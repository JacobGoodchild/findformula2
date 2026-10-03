import numpy as np
from mpmath import mp, mpf, sqrt, quad, pi, asin, polyroots, cos, sin
mp.dps=30
def W_coeffs(a,b):
    # W = -(r^4 + (6b-6) r^3 + (13b^2-22b+9) r^2 + (12b^3-24b^2+12b-4a^2) r + 4b^4-8b^3+4b^2-4a^2 b)
    return [-1,-(6*b-6),-(13*b*b-22*b+9),-(12*b**3-24*b*b+12*b-4*a*a),-(4*b**4-8*b**3+4*b*b-4*a*a*b)]
def A_reduced(a,b):
    a=mpf(a); b=mpf(b); co=W_coeffs(a,b); W=lambda r: sum(co[i]*r**(4-i) for i in range(5))
    rts=sorted([x.real for x in polyroots(co,maxsteps=500,extraprec=300) if abs(x.imag)<mpf(10)**-20])
    lo,hi=rts[-2],rts[-1]
    P=lambda r: -2*(b-1)*(b*b-8*b+3)+(-3*b*b+32*b-21)*r+(15-b)*r*r-2*(b+3)*(a-b+1)*(a+b-1)/(r+b-1)
    f=lambda t:(lambda r: P(r)/sqrt(W(r))*(hi-lo)*sin(t)/2)(lo+(hi-lo)*(1-cos(t))/2)
    return quad(f,[0,pi/2,pi]).real
def A_ext(a):
    a=mpf(a); return sqrt(a)*sqrt(4-a)*(a+14)+8*(2*a+1)*asin(sqrt(a)/2)+4*pi*(2*a+1)
def shoelace(a,b,N=1500000,ext=False):
    co=[float(x) for x in W_coeffs(a,b)]; rts=np.sort(np.roots(co).real)
    lo,hi=rts[-2],rts[-1]
    if ext: lo=max(lo,1-b+1e-5)
    t=np.linspace(0,np.pi,N); r=lo+(hi-lo)*(1-np.cos(t))/2
    xi=-(a*a*b+a*a*r+a*a+2*b*b*r+3*b*r*r-2*b*r+r**3-3*r*r)/(a*(b+r-1))
    W=np.polyval(co,r); eta=r*r*W/(a*a*(b+r-1)**2)
    X=np.concatenate([-xi,(-xi)[::-1]]); Y=np.concatenate([np.sqrt(np.clip(eta,0,None)),-np.sqrt(np.clip(eta,0,None))[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
for a,b in [(0.5,0.2),(0.3,0.4),(0.8,0.1)]:
    print('general',a,b, A_reduced(a,b), shoelace(a,b))
for a in [0.5,0.8,0.2]:
    print('extremal',a, A_ext(a), shoelace(a,1-a,ext=True))
