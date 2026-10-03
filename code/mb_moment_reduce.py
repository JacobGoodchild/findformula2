import sympy as sp
from reduce_even import reduce, u, a, P
y=1+u
L=-(y**4-2*y**3+a**2)/(a*u)
dL=sp.diff(L,u)
R=sp.together(-L*2*(y**2/(a*u))*dL*P)    # moment integrand times sqrt(P):  Int L * 2 sqrt(Q) * (-L') du
res=reduce(R,pmax=13,nmax=6)
print(sp.simplify(res))
