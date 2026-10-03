import sympy as sp
from reduce_even2 import reduce, u, a, B, P
b=1-B; y=1+u
L=-(y**4-2*y**3+2*b*y+a**2-b**2)/(a*u)
R=L*2*(y**2-b)*P*(3*u**4+4*u**3+B**2-a**2)/(a**2*u**3)
res=reduce(R,pmax=14,nmax=6)
print(sp.simplify(res))
