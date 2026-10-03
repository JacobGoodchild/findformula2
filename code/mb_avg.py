# Direction-averaged slow-particle capture cross-section: <sigma> v^2 = (pi/2) Int_{-1}^{1} C dmu along marginally bound orbits
import sympy as sp
from mpmath import mp, mpf, sqrt, quad, pi, diff
y,a=sp.symbols('y a',positive=True)
C=(3*y**4-4*y**3+a**2)/(y-1)**2
mu=-(y**4-2*y**3+a**2)/(a*sp.sqrt(3*y**4-4*y**3+a**2))
integrand=sp.simplify(C*sp.diff(mu,y))
print(sp.factor(integrand))
f=sp.lambdify((y,a),integrand,'mpmath')
mp.dps=30
def avg(av):
    av=mpf(av); y1=1+sqrt(1-av); y2=1+sqrt(1+av)
    I=quad(lambda yy: f(yy,av),[y1,(y1+y2)/2,y2])
    return pi/2*abs(I)
for av in ['0.001','0.5','0.9','1']:
    print(av, avg(av))
