import sympy as sp
x,s=sp.symbols('x s',positive=True)
P=sp.expand(-((x+s)**2-2+2*s)*((x-s)**2-2-2*s))
print('P=',sp.collect(P,x))
c0,c1,c2,s0,s1,s2=sp.symbols('c0 c1 c2 s0 s1 s2')
S=s0+s1*x+s2*x**2
eq=sp.expand(x*P-(c0+c1*x+c2*x**2+sp.diff(S,x)*P+S*sp.diff(P,x)/2))
sol=sp.solve(sp.Poly(eq,x).coeffs(),[c0,c1,c2,s0,s1,s2],dict=True)[0]
print({k:sp.factor(v) for k,v in sol.items()})
