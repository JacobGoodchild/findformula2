# First moment of the slow-particle capture region (edge-on): M1 = Int Int L_z dL_z d(sqrt Q)  (over capture region)
import sympy as sp
from mpmath import mp, mpf, sqrt, quad, pi, ellipe, ellipk, cos, sin
u,a=sp.symbols('u a',positive=True)
y=1+u
L=-(y**4-2*y**3+a**2)/(a*u)
sqQ_over_sqrtP=y**2/(a*u)          # sqrt(Q) = y^2 sqrt(P)/(a u),  P=a^2-(u^2-1)^2
dL=sp.diff(L,u)
P=a**2-(u**2-1)**2
# area element: 2*sqrt(Q)*|dL|  -> moment integrand L*2*sqrt(Q)*dL  (orientation: L decreasing/increasing; take abs at end)
integrand=sp.together(sp.expand(L*2*sqQ_over_sqrtP*dL*P))   # times 1/sqrt(P)
print('numerator degree', sp.degree(sp.numer(integrand),u), 'denominator', sp.factor(sp.denom(integrand)))
f=sp.lambdify((u,a),integrand,'mpmath')
mp.dps=30
def M1(av):
    av=mpf(av); lo,hi=sqrt(1-av),sqrt(1+av)
    g=lambda t:(lambda uu: f(uu,av)/sqrt(av*av-(uu*uu-1)**2)*(hi-lo)*sin(t)/2)(lo+(hi-lo)*(1-cos(t))/2)
    return quad(g,[0,pi/2,pi])
for av in ['0.1','0.5','0.9']:
    print(av, M1(av))
