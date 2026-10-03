from mpmath import mp, mpf, sqrt, quad, cos, acos, pi
from shadow_area import area, roots
mp.dps=30
def Ared(a):
    r1,r2=roots(a)
    W=lambda r: r*(4*a**2-r*(r-3)**2)
    f=lambda r:(6-21*r+15*r**2+6*(1-a**2)/(r-1))/sqrt(W(r))
    return quad(f,[r1,(r1+r2)/2,r2])
for a in [mpf('0.1'),mpf('0.5'),mpf('0.9'),mpf('0.99')]:
    print(a, Ared(a), -area(a))
