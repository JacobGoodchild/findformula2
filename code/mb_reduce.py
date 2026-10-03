import sympy as sp
u,a=sp.symbols('u a',positive=True)
P=sp.expand(a**2-(u**2-1)**2)
R=sp.expand(2/a**2*(1+u)**2*(3*u**4+4*u**3+1-a**2)*P)/u**3   # integrand = R/sqrt(P)
cs=sp.symbols('c0:4'); ds=sp.symbols('d1:4'); ss=sp.symbols('s0:5'); ts=sp.symbols('t1:3')
S=sum(ss[i]*u**i for i in range(5))+ts[0]/u+ts[1]/u**2
rhs=sum(cs[i]*u**i for i in range(4))+ds[0]/u+ds[1]/u**2+ds[2]/u**3+sp.diff(S,u)*P+S*sp.diff(P,u)/2
num=sp.numer(sp.together(R-rhs))
unk=list(cs)+list(ds)+list(ss)+list(ts)
sol=sp.solve(sp.Poly(sp.expand(num),u).coeffs(),unk,dict=True)
for so in sol:
    print({k:sp.factor(v) for k,v in so.items()})
