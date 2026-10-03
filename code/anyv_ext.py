import sympy as sp
r,L,Q=sp.symbols('r L Q',real=True)
E,p,t=sp.symbols('E p t',positive=True)
a=1
D=r**2-2*r+a**2
R=sp.expand((E*(r**2+a**2)-a*L)**2-D*(r**2+(L-a*E)**2+Q))
sol=sp.solve([R,sp.diff(R,r)],[L,Q],dict=True)
for s_ in sol:
    Ls=sp.simplify(s_[L]); Qs=sp.simplify(s_[Q])
    print('L=',sp.factor(Ls)); print('Q=',sp.factor(Qs)); print('--')
