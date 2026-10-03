# Generic Hermite reduction for Int R(u) du / sqrt(P), P = a^2 - (u^2-1)^2, over [sqrt(1-a), sqrt(1+a)]
import sympy as sp
u,a=sp.symbols('u a',positive=True)
P=sp.expand(a**2-(u**2-1)**2)
A=sp.sqrt(1+a); B=sp.sqrt(1-a); m=2*a/(1+a)
K,E=sp.symbols('K E')   # complete elliptic integrals of parameter m
basis_val={0:K/A, 1:sp.pi/2, 2:A*E, 3:sp.pi/2, -1:sp.pi/(2*sp.sqrt(1-a**2)), -2:E/(A*B**2), -3:sp.pi/(2*(1-a**2)**sp.Rational(3,2))}
def reduce(R, pmax=10, nmax=8):
    R=sp.together(R)
    cs={k:sp.Symbol('c%d'%k) for k in [0,1,2,3,-1,-2,-3]}
    ss=[sp.Symbol('s%d'%k) for k in range(pmax)]; ts=[sp.Symbol('t%d'%k) for k in range(1,nmax)]
    S=sum(ss[k]*u**k for k in range(pmax))+sum(ts[j-1]*u**(-j) for j in range(1,nmax))
    rhs=sum(cs[k]*u**k for k in cs)+sp.diff(S,u)*P+S*sp.diff(P,u)/2
    num=sp.numer(sp.together(R-rhs))
    sol=sp.solve(sp.Poly(sp.expand(num),u).coeffs(),list(cs.values())+ss+ts,dict=True)[0]
    tot=0
    for k,c in cs.items():
        val=sol.get(c,c)
        val=val.subs({x:0 for x in ss+ts if x not in sol})   # free params -> 0
        tot+=val*basis_val[k]
    # boundary terms vanish at both ends (P=0) since u>0
    return sp.simplify(tot)
if __name__=="__main__":
    # test on the area integrand
    R=2/a**2*(1+u)**2*(3*u**4+4*u**3+1-a**2)*P/u**3
    print(sp.simplify(reduce(R)))
