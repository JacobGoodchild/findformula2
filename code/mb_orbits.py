import sympy as sp
r,a,L,Q=sp.symbols('r a L Q',real=True)
D=r**2-2*r+a**2
R=sp.expand((r**2+a**2-a*L)**2-D*(r**2+(L-a)**2+Q))
sol=sp.solve([R,sp.diff(R,r)],[L,Q],dict=True)
for s_ in sol: print(sp.factor(sp.simplify(s_[L])),' | ',sp.factor(sp.simplify(s_[Q])))
x=sp.symbols('x',real=True)
Lx=sp.factor(sp.simplify(sol[0][L].subs(sp.sqrt(r**5),x**5).subs(r,x**2)))
Qx=sp.factor(sp.simplify(sol[0][Q].subs(sp.sqrt(r**5),x**5).subs(r,x**2)))
print('L(x)=',Lx); print('Q(x)=',Qx)
print("dL/dx=",sp.factor(sp.diff(Lx,x)))
