# High-precision capture area for Kerr, energy E, incidence theta (analytic dL/dr)
import sympy as sp
from mpmath import mp, mpf, sqrt, quad, pi, cos, sin, findroot
r,a,E,L,Q=sp.symbols('r a E L Q',positive=True)
D=r**2-2*r+a**2
R=sp.expand((E*(r**2+a**2)-a*L)**2-D*(r**2+(L-a*E)**2+Q))
sol=sp.solve([R,sp.diff(R,r)],[L,Q],dict=True)[1]
Lf=sp.lambdify((r,a,E),sol[L],'mpmath'); Qf=sp.lambdify((r,a,E),sol[Q],'mpmath')
dLf=sp.lambdify((r,a,E),sp.diff(sol[L],r),'mpmath')
def area(Ev,av,th):
    Ev=mpf(Ev); av=mpf(av); p2=Ev*Ev-1; c2=cos(th)**2; s2=sin(th)**2
    g=lambda rr: Qf(rr,av,Ev)+av*av*p2*c2-Lf(rr,av,Ev)**2*c2/s2
    E2=Ev*Ev; disc=(3*E2-4)**2-16*(1-E2)
    rc=[x for x in [(-(3*E2-4)+sg*sqrt(disc))/(2*(1-E2)) for sg in (1,-1)] if 3<x.real<=4][0]
    w=6*av; rs=[rc-w+2*w*k/2000 for k in range(2001)]
    idx=[i for i,x in enumerate(rs) if g(x).real>0]
    lo=findroot(g,(rs[idx[0]-1],rs[idx[0]]),solver='anderson'); hi=findroot(g,(rs[idx[-1]],rs[idx[-1]+1]),solver='anderson')
    f=lambda t:(lambda x: sqrt(max(g(x).real,0))*abs(dLf(x,av,Ev))*(hi-lo)*sin(t)/2)(lo+(hi-lo)*(1-cos(t))/2)
    return 2/sin(th)*quad(f,[0,pi/2,pi])/p2
