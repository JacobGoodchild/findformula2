import sympy as sp
r,a=sp.symbols('r a',positive=True)
e2=1-a**2
xi=-(r**3-3*r**2+a**2*r+a**2+2*e2*r)/(a*(r-1))
etaN=(4*a**2*(r-e2)-(r**2-3*r+2*e2)**2)
print('xi =',sp.factor(sp.simplify(xi)))
print('eta numerator =',sp.factor(etaN))
print("xi' =",sp.factor(sp.diff(sp.simplify(xi),r)))
