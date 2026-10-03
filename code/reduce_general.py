# Hermite-reduce Kerr shadow area integrand for general a, c=cos(theta):
# A = (4/(a^2 s^2)) Int sqrt(N) ((r-1)^3+1-a^2)/(r-1)^3 dr
import sympy as sp
r,a,c=sp.symbols('r a c',positive=True)
s2=1-c**2
N=sp.expand(s2*r**3*(4*a**2-r*(r-3)**2)+a**4*c**2*s2*(r-1)**2-c**2*(r**2*(3-r)-a**2*(r+1))**2)
R=4/(a**2*s2)*N*((r-1)**3+1-a**2)/(r-1)**3      # integrand = R/sqrt(N)
cs=sp.symbols('c0:5'); ds=sp.symbols('d1:4'); ss=sp.symbols('s0:3'); ts=sp.symbols('t1:3')
S=sum(ss[i]*r**i for i in range(3))+ts[0]/(r-1)+ts[1]/(r-1)**2
rhs=sum(cs[i]*r**i for i in range(5))+ds[0]/(r-1)+ds[1]/(r-1)**2+ds[2]/(r-1)**3+sp.diff(S,r)*N+S*sp.diff(N,r)/2
num=sp.numer(sp.together(R-rhs))
eqs=sp.Poly(sp.expand(num),r).coeffs()
unk=list(cs)+list(ds)+list(ss)+list(ts)
sol=sp.solve(eqs,unk,dict=True)
for so in sol:
    for k in unk: print(k, sp.factor(so.get(k,k)))
so=sol[0]
t1,t2=ts
t2v=sp.solve(sp.Eq(so[ds[2]],0),t2)[0]
t1v=sp.solve(sp.Eq(so[ds[1]].subs(t2,t2v),0),t1)[0]
print('t2=',sp.factor(t2v),' t1=',sp.factor(t1v))
for k in list(cs)+[ds[0]]:
    print(k,'=',sp.factor(sp.simplify(so[k].subs(t2,t2v).subs(t1,t1v))))
