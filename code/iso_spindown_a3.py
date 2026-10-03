import sympy as sp
a,r,c=sp.symbols('a r c')
C=c**2; s2=1-C
sig=1-a**2*(r+3*C*(4-r))/(2*r**2*(6-r))
ell=-2*a*s2/(r-2)+a**3*s2*(sp.Rational(3,2)*(r-5)-C*(7*r**2-63*r+144)/(2*r))/((r-2)*(6-r)**3)
num=sp.integrate(sp.expand(sig*ell),(c,0,1)); den=sp.integrate(sig,(c,0,1))
L=sp.series(num/den,a,0,4).removeO()
L1=sp.factor(L.coeff(a,1)); L3=sp.factor(L.coeff(a,3))
print('L1',L1); print('L3',L3)
print('da/dlnM = a*(-2+L1) + a^3*L3')
for rv in [4,3]: print(rv, sp.nsimplify(L1.subs(r,rv)), sp.nsimplify(L3.subs(r,rv)))
