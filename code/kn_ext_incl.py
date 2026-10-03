from mpmath import mp, mpf, sqrt, quad, sin, cos, pi
import numpy as np
mp.dps=30
def A_formula(a,th):
    a=mpf(a); s=sin(th); c2=cos(th)**2
    N=lambda r: -(r*r-2*(1+a*s)*r-a*a*c2)*(r*r-2*(1-a*s)*r-a*a*c2)
    rD=1+a*s+sqrt(1+a*a+2*a*s); rC=1-a*s+sqrt(1+a*a-2*a*s)
    lo=max(rC,mpf(1))
    I=4*quad(lambda r:(3*r*r+2*(a*a-1)*r-a*a*c2)/sqrt(N(r)),[lo,(lo+rD)/2,rD])
    B=0 if rC>=1 else (1+a*a+a*a*s*s)*sqrt(N(mpf(1)))/(a*a*s*s)
    return I+B
def A_direct(a,th):
    # general KN critical curve, e^2=1-a^2, observer inclination th
    a=mpf(a); e2=1-a*a; s=sin(th); c=cos(th)
    xi=lambda r: -(r**3-3*r**2+a*a*r+a*a+2*e2*r)/(a*(r-1))
    eta=lambda r: r**2*(4*a*a*(r-e2)-(r*r-3*r+2*e2)**2)/(a*a*(r-1)**2)
    b2=lambda r: eta(r)+a*a*c*c-xi(r)**2*c*c/(s*s)
    rD=1+a*s+sqrt(1+a*a+2*a*s); rC=1-a*s+sqrt(1+a*a-2*a*s)
    lo=max(rC,1+mpf(10)**-12)
    return 2*quad(lambda r: sqrt(max(b2(r),0))*2*abs(r-1)/(a*s),[lo,(lo+rD)/2,rD])
for a,deg in [(0.3,40),(0.6,30),(0.6,80),(0.8,60),(1,60),(0.4,89.9)]:
    th=mpf(deg)*pi/180
    print(a,deg, A_formula(a,th), A_direct(a,th))
