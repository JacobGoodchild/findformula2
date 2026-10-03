# Kerr-Newman edge-on shadow area: Hermite reduction
import sympy as sp
r,a,q=sp.symbols('r a q',positive=True)   # q = e^2
xi=-(r**3-3*r**2+a**2*r+a**2+2*q*r)/(a*(r-1))
W=sp.expand(4*a**2*(r-q)-(r**2-3*r+2*q)**2)       # eta = r^2 W/(a^2 (r-1)^2)
dxi=sp.factor(sp.diff(xi,r)); print('dxi=',dxi)
# area = 2 Int sqrt(eta)|dxi| dr = 2 Int r sqrt(W)/(a (r-1)) * |dxi| dr
num=sp.factor(sp.numer(sp.together(-dxi))); print('-dxi numerator', num)
R=sp.together(2*r/(a*(r-1))*(-dxi)*W)    # integrand = R/sqrt(W)
cs=sp.symbols('c0:3'); ds=sp.symbols('d1:4'); ss=sp.symbols('s0:3'); ts=sp.symbols('t1:3')
S=ss[0]+ss[1]*r+ss[2]*r**2+ts[0]/(r-1)+ts[1]/(r-1)**2
rhs=cs[0]+cs[1]*r+cs[2]*r**2+ds[0]/(r-1)+ds[1]/(r-1)**2+ds[2]/(r-1)**3+sp.diff(S,r)*W+S*sp.diff(W,r)/2
eq=sp.numer(sp.together(R-rhs))
sol=sp.solve(sp.Poly(sp.expand(eq),r).coeffs(),list(cs)+list(ds)+list(ss)+list(ts),dict=True)
for so in sol:
    for k_,v in so.items(): print(k_, sp.factor(v))
so=sol[0]; t1,t2=ts
t2v=sp.solve(so[ds[2]],t2)[0]
t1v=sp.solve(so[ds[1]].subs(t2,t2v),t1)[0]
print('t2=',sp.factor(t2v),' t1=',sp.factor(t1v))
for k_ in list(cs)+[ds[0]]:
    print(k_,'=',sp.factor(sp.simplify(so[k_].subs(t2,t2v).subs(t1,t1v))))
