from mpmath import mp, mpf, sqrt, quad, pi, sin, cos, findroot, diff
from mb_incl import A_mb
mp.dps=30
def sigma_master(a,th):
    a=mpf(a); s=sin(th)
    C=lambda y:(3*y**4-4*y**3+a*a)/(y-1)**2
    m=lambda y:-(y**4-2*y**3+a*a)/(a*sqrt(3*y**4-4*y**3+a*a))
    y1=1+sqrt(1-a); y2=1+sqrt(1+a)
    # find y range where |mu|<=s  (mu decreasing in y)
    ya=findroot(lambda y: m(y)-s,(y1,y2),solver='bisect') if s<1 else y1
    yb=findroot(lambda y: m(y)+s,(y1,y2),solver='bisect') if s<1 else y2
    def g(y):
        w=s*s-m(y)**2
        return C(y)*(-diff(m,y))/sqrt(w) if w>0 else mpf(0)
    f=lambda t: g(ya+(yb-ya)*(1-cos(t))/2)*(yb-ya)*sin(t)/2
    return quad(f,[0,pi/2,pi])
for a,deg in [('0.5',90),('0.5',60),('0.9',30),('0.99',75)]:
    th=mpf(deg)*pi/180
    print(a,deg, sigma_master(a,th), A_mb(a,th))
