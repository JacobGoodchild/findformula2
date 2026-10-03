# Kerr-Newman capture cross-section at any speed: order a^2 correction around Reissner-Nordstrom.
import sympy as sp
a,ci,phi,S,C,q,e,l,r=sp.symbols('a c_i phi S C q e l r')
rc=sp.symbols('r_c')
e2=(rc**2-2*rc+q**2)**2/(rc**2*(rc**2-3*rc+2*q**2)); l2=rc**2*(rc-q**2)/(rc**2-3*rc+2*q**2)
def red(x):
    x=sp.expand(x); out=0
    for (ie,il),c in sp.Poly(x,e,l).terms():
        out+=c*e**(ie%2)*l**(il%2)*e2**(ie//2)*l2**(il//2)
    return out
N=2
Ls=sp.symbols('L1:%d'%(N+1)); Rs=sp.symbols('R1:%d'%(N+1))
Lt=l+sum(Ls[k]*a**(k+1) for k in range(N)); rt=rc+sum(Rs[k]*a**(k+1) for k in range(N))
Lz=Lt*ci; Qc=Lt**2*(1-ci**2); Dl=r**2-2*r+a**2+q**2
Rpot=(e*(r**2+a**2)-a*Lz)**2-Dl*(r**2+(Lz-a*e)**2+Qc)
dR=sp.diff(Rpot,r)
E1=sp.expand(Rpot.subs(r,rt)); E2=sp.expand(dR.subs(r,rt))
print('order0', sp.simplify(red(E1.coeff(a,0))), sp.simplify(red(E2.coeff(a,0))))
sol={}
for k in range(1,N+1):
    c1=red(sp.expand(E1.subs(sol)).coeff(a,k)); c2=red(sp.expand(E2.subs(sol)).coeff(a,k))
    s_=sp.solve([c1,c2],[Ls[k-1],Rs[k-1]],dict=True)[0]
    sol.update({kk:sp.factor(v) for kk,v in s_.items()})
    print('L%d ='%k, sol[Ls[k-1]], flush=True)
Lc=Lt.subs(sol)
L2s=sp.expand(sp.series(Lc**2,a,0,N+1).removeO())
co=sp.integrate(sp.expand(L2s.coeff(a,2).subs(ci,S*sp.cos(phi))),(phi,0,2*sp.pi))/(2*sp.pi)
co=red(sp.together(co))
co=sp.factor(sp.simplify(sp.expand(co/l2).subs(S**2,1-C)))
print('a^2 relative coefficient:', co)
print('q=0:', sp.factor(co.subs(q,0)))
