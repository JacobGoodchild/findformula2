import sympy as sp
r,a,L,Q,q=sp.symbols('r a L Q q',positive=True)
D=r**2-2*r+a**2+q
R=sp.expand((r**2+a**2-a*L)**2-D*(r**2+(L-a)**2+Q))
sol=sp.solve([R,sp.diff(R,r)],[L,Q],dict=True)
for s_ in sol:
    print(sp.factor(sp.simplify(s_[L]))); print(sp.factor(sp.simplify(s_[Q]))); print('--')
y=sp.symbols('y',positive=True)
s1=sol[1]
Ly=sp.factor(sp.simplify(s1[L].subs(r,y**2)))
Qy=sp.factor(sp.simplify(s1[Q].subs(r,y**2)))
print('L(y)=',Ly); print('Q(y)=',Qy)
print("L'(y)=",sp.factor(sp.diff(Ly,y)))
