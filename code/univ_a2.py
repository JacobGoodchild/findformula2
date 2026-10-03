# Universal a^2 correction to capture cross-section / shadow area for rotating metrics of
# Kerr-like (Newman-Janis / Azreg-Ainou) form: Delta = r^2 f(r) + a^2 (Kerr: f = 1 - 2/r).
# Radial potential: R = (E(r^2+a^2) - a Lz)^2 - Delta (mu r^2 + (Lz - aE)^2 + Q).
import sympy as sp, sys
mu=int(sys.argv[1]) if len(sys.argv)>1 else 1
a,ci,phi,S,C,e,l,r,rc=sp.symbols('a c_i phi S C e l r r_c')
f0,f1,f2,f3=sp.symbols('f0 f1 f2 f3')
F=f0+f1*(r-rc)+f2*(r-rc)**2/2+f3*(r-rc)**3/6
if mu==1:
    e2=2*f0**2/(2*f0-rc*f1); l2=f1*rc**3/(2*f0-rc*f1); rel={}
else:
    e2=sp.Integer(1); l2=rc**2/f0; rel={f1:2*f0/rc}   # photon sphere: r f' = 2 f, E=1, L=b
def red(x):
    x=sp.expand(x); out=0
    for (ie,il),c in sp.Poly(x,e,l).terms():
        out+=c*e**(ie%2)*l**(il%2)*e2**(ie//2)*l2**(il//2)
    return sp.expand(out.subs(rel)) if rel else out
N=2
Ls=sp.symbols('L1:%d'%(N+1)); Rs=sp.symbols('R1:%d'%(N+1))
Lt=l+sum(Ls[k]*a**(k+1) for k in range(N)); rt=rc+sum(Rs[k]*a**(k+1) for k in range(N))
Lz=Lt*ci; Qc=Lt**2*(1-ci**2); Dl=r**2*F+a**2
Rpot=(e*(r**2+a**2)-a*Lz)**2-Dl*(mu*r**2+(Lz-a*e)**2+Qc)
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
p2=e2-1 if mu==1 else sp.Integer(1)
tot=sp.factor(sp.simplify(sp.expand(co/l2).subs(S**2,1-C).subs(rel) + C*p2/l2.subs(rel)))
print('TOTAL relative a^2 coefficient:', tot)
print('collected:', sp.collect(sp.expand(sp.numer(tot)),C), ' / ', sp.factor(sp.denom(tot)))
