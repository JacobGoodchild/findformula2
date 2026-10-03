import sympy as sp
r,a,b,xi,eta=sp.symbols('r a b xi eta',real=True)
rho=r*(r+2*b)
D=rho-2*r+a**2
R=(rho+a**2-a*xi)**2-D*(eta+(xi-a)**2)
sol=sp.solve([R,sp.diff(R,r)],[xi,eta],dict=True)
for s_ in sol:
    X=sp.factor(sp.simplify(s_[xi])); H=sp.factor(sp.simplify(s_[eta]))
    print('xi=',X); print('eta=',H); print("xi'=",sp.factor(sp.diff(X,r))); print('---')
