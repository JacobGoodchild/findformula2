import sympy as sp
r,a,c,q=sp.symbols('r a c q',positive=True)
s2=1-c**2
xi=-(r**3-3*r**2+a**2*r+a**2+2*q*r)/(a*(r-1))
etaN=r**2*(4*a**2*(r-q)-(r**2-3*r+2*q)**2)      # eta = etaN/(a^2 (r-1)^2)
# beta^2 * a^2 (r-1)^2 s^2 = N
N=sp.expand(s2*etaN+a**4*c**2*s2*(r-1)**2-c**2*(r**3-3*r**2+a**2*r+a**2+2*q*r)**2)
print('N =',sp.collect(N,r))
# area = 2 Int beta |dxi|/s dr, beta = sqrt(N)/(a s (r-1)), |dxi| = 2((r-1)^3+1-a^2-q)/(a (r-1)^2)
R=4/(a**2*s2)*N*((r-1)**3+1-a**2-q)/(r-1)**3
cs=sp.symbols('c0:5'); ds=sp.symbols('d1:4'); ss=sp.symbols('s0:3'); ts=sp.symbols('t1:3')
S=sum(ss[i]*r**i for i in range(3))+ts[0]/(r-1)+ts[1]/(r-1)**2
rhs=sum(cs[i]*r**i for i in range(5))+ds[0]/(r-1)+ds[1]/(r-1)**2+ds[2]/(r-1)**3+sp.diff(S,r)*N+S*sp.diff(N,r)/2
num=sp.numer(sp.together(R-rhs))
unk=list(cs)+list(ds)+list(ss)+list(ts)
sol=sp.solve(sp.Poly(sp.expand(num),r).coeffs(),unk,dict=True)[0]
t1,t2=ts
t2v=sp.solve(sol[ds[2]],t2)[0]
t1v=sp.solve(sol[ds[1]].subs(t2,t2v),t1)[0]
print('t2=',sp.factor(t2v),' t1=',sp.factor(t1v))
for k_ in list(cs)+[ds[0]]+list(ss):
    print(k_,'=',sp.factor(sp.simplify(sol[k_].subs(t2,t2v).subs(t1,t1v))))
