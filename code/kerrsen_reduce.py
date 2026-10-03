import sympy as sp
r,a,b=sp.symbols('r a b',positive=True)
xi=-(a**2*b+a**2*r+a**2+2*b**2*r+3*b*r**2-2*b*r+r**3-3*r**2)/(a*(b+r-1))
W=sp.expand(-(-4*a**2*b-4*a**2*r+4*b**4+12*b**3*r-8*b**3+13*b**2*r**2-24*b**2*r+4*b**2+6*b*r**3-22*b*r**2+12*b*r+r**4-6*r**3+9*r**2))
# eta = r^2 W/(a^2 (r+b-1)^2); area = 2 Int r sqrt(W)/(a (r+b-1)) |xi'| dr
dxi=sp.diff(xi,r)
R=sp.together(2*r/(a*(r+b-1))*(-dxi)*W)
s_=r+b-1
cs=sp.symbols('c0:3'); ds=sp.symbols('d1:4'); ss=sp.symbols('s0:3'); ts=sp.symbols('t1:3')
S=ss[0]+ss[1]*r+ss[2]*r**2+ts[0]/s_+ts[1]/s_**2
rhs=cs[0]+cs[1]*r+cs[2]*r**2+ds[0]/s_+ds[1]/s_**2+ds[2]/s_**3+sp.diff(S,r)*W+S*sp.diff(W,r)/2
num=sp.numer(sp.together(R-rhs))
sol=sp.solve(sp.Poly(sp.expand(num),r).coeffs(),list(cs)+list(ds)+list(ss)+list(ts),dict=True)[0]
t1,t2=ts
t2v=sp.solve(sol[ds[2]],t2)[0]; t1v=sp.solve(sol[ds[1]].subs(t2,t2v),t1)[0]
for k in list(cs)+[ds[0]]:
    print(k,'=',sp.factor(sp.simplify(sol[k].subs(t2,t2v).subs(t1,t1v))))
