import sympy as sp
r,a,b,L,Q,y=sp.symbols('r a b L Q y',positive=True)
rho=r*(r+2*b); D=rho-2*r+a**2
R=sp.expand((rho+a**2-a*L)**2-D*(rho+(L-a)**2+Q))
sol=sp.solve([R,sp.diff(R,r)],[L,Q],dict=True)
for s_ in sol:
    Ly=sp.factor(sp.simplify(s_[L].subs(r,y**2-b))); Qy=sp.factor(sp.simplify(s_[Q].subs(r,y**2-b)))
    print('L=',Ly); print('Q=',Qy); print("L'=",sp.factor(sp.diff(Ly,y))); print('---')
