import sympy as sp
r,a,s=sp.symbols('r a s',positive=True)
N4=sp.expand(-((r-a*(1+s))**2-2*a*(1+s))*((r-a*(1-s))**2-2*a*(1-s)))
c0,c1,c2,s0,s1,s2=sp.symbols('c0 c1 c2 s0 s1 s2')
S=s0+s1*r+s2*r**2
eq=sp.expand((r-a)*N4-(c0+c1*r+c2*r**2+sp.diff(S,r)*N4+S*sp.diff(N4,r)/2))
sol=sp.solve(sp.Poly(eq,r).coeffs(),[c0,c1,c2,s0,s1,s2],dict=True)[0]
print({k:sp.factor(v) for k,v in sol.items()})
