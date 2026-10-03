import sympy as sp
r,a,L,Q,E=sp.symbols('r a L Q E',positive=True)
D=r**2-2*r+a**2
R=sp.expand((E*(r**2+a**2)-a*L)**2-D*(r**2+(L-a*E)**2+Q))
sol=sp.solve([R,sp.diff(R,r)],[L,Q],dict=True)
pass
t,p=sp.symbols('t p',positive=True)
s0=sol[1]
rr=1/(t**2-p**2)
def sub(ex):
    ex=ex.subs(sp.sqrt(E**2*r-r+1),sp.sqrt(r)*t*sp.sqrt(r)/sp.sqrt(r))  # placeholder
    return ex
# do it manually: r^(5/2)*sqrt(E^2 r - r + 1) = r^2 * sqrt(r(E^2 r - r+1)) = r^2 * r t = r^3 t
Lx=s0[L]; Qx=s0[Q]
Lx=Lx.subs(r**sp.Rational(5,2)*sp.sqrt(E**2*r-r+1), r**3*t)
Qx=Qx.subs(r**sp.Rational(5,2)*sp.sqrt(E**2*r-r+1), r**3*t)
Lt=sp.factor(sp.simplify(Lx.subs(r,rr).subs(E,sp.sqrt(1+p**2))))
Qt=sp.factor(sp.simplify(Qx.subs(r,rr).subs(E,sp.sqrt(1+p**2))))
print('L(t)=',Lt); print('Q(t)=',Qt)
num=sp.numer(sp.together(Qt))
for pv in [sp.Rational(3,4), sp.Rational(12,5)]:
    nn=sp.expand(num.subs(p,pv))
    print(pv, sp.factor(nn))
