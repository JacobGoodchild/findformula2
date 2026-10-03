# Hermite reduction for Int R(u) du/sqrt(P), P = a^2 - (u^2 - B)^2 over [sqrt(B-a), sqrt(B+a)]
import sympy as sp
u,a,B=sp.symbols('u a B',positive=True)
P=sp.expand(a**2-(u**2-B)**2)
Ap=sp.sqrt(B+a); Bm=sp.sqrt(B-a)
K,E=sp.symbols('K E')   # complete elliptic integrals with parameter m = 2a/(B+a)
# basis values: Int u^k du/sqrt(P)
basis_val={0:K/Ap, 2:Ap*E, -2:E/(Ap*Bm**2), 1:sp.pi/2, 3:B*sp.pi/2, -1:sp.pi/(2*sp.sqrt(B**2-a**2)), -3:B*sp.pi/(2*(B**2-a**2)**sp.Rational(3,2))}
def reduce(R,pmax=10,nmax=6):
    R=sp.together(R)
    cs={k:sp.Symbol('c%d'%(k+10)) for k in basis_val}
    ss=[sp.Symbol('s%d'%k) for k in range(pmax)]; ts=[sp.Symbol('t%d'%k) for k in range(1,nmax)]
    S=sum(ss[k]*u**k for k in range(pmax))+sum(ts[j-1]*u**(-j) for j in range(1,nmax))
    rhs=sum(cs[k]*u**k for k in cs)+sp.diff(S,u)*P+S*sp.diff(P,u)/2
    num=sp.numer(sp.together(R-rhs))
    sol=sp.solve(sp.Poly(sp.expand(num),u).coeffs(),list(cs.values())+ss+ts,dict=True)[0]
    tot=0
    for k,c in cs.items():
        val=sol.get(c,c); val=val.subs({x:0 for x in ss+ts if x not in sol}); tot+=val*basis_val[k]
    return sp.simplify(tot)
