import sympy as sp
r,a,e,xi,eta=sp.symbols('r a e xi eta',real=True)
D=r**2-2*r+a**2+e**2
R=(r**2+a**2-a*xi)**2-D*(eta+(xi-a)**2)
sol=sp.solve([R,sp.diff(R,r)],[xi,eta],dict=True)
for s_ in sol: print(sp.factor(s_[xi]),'|',sp.factor(s_[eta]))
