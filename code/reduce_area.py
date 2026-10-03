# Reduce A(a) = (4/a^2) * Int_{r1}^{r2} r sqrt(W) [1 + (1-a^2)/(r-1)^3] dr,  W = r(4a^2 - r(r-3)^2)
# to basis I0=Int dr/sqrtW, I1=Int r dr/sqrtW, I2=Int r^2 dr/sqrtW, J=Int dr/((r-1) sqrtW)
import sympy as sp
r,a=sp.symbols('r a',positive=True)
W=sp.expand(r*(4*a**2-r*(r-3)**2))
R=sp.together(4/a**2*r*W*(1+(1-a**2)/(r-1)**3))   # integrand = R/sqrtW
c0,c1,c2,d,s0,s1,s2,s3,t1,t2=sp.symbols('c0 c1 c2 d s0 s1 s2 s3 t1 t2')
S=s0+s1*r+s2*r**2+s3*r**3+t1/(r-1)+t2/(r-1)**2
rhs=c0+c1*r+c2*r**2+d/(r-1)+sp.diff(S,r)*W+S*sp.diff(W,r)/2
eq=sp.expand(sp.cancel((R-rhs)*(r-1)**3))
sol=sp.solve(sp.Poly(sp.numer(sp.together(eq)),r).coeffs(),[c0,c1,c2,d,s0,s1,s2,s3,t1,t2],dict=True)
print(sol)
