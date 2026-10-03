# Analytic perturbation with rational parametrisation r_c=(u^2+3)^2/(4u^2):  sqrt(r_c)=(u^2+3)/(2u), sqrt(r_c-3)=(3-u^2)/(2u)
import sympy as sp
a,ci,phi,S,u=sp.symbols('a c_i phi S u',positive=True)
r=sp.symbols('r')
rc=(u**2+3)**2/(4*u**2); src=(u**2+3)/(2*u); srm3=(3-u**2)/(2*u)
E=(rc-2)/(src*srm3); L0=rc/srm3
N=4
Ls=sp.symbols('L1:%d'%(N+1)); Rs=sp.symbols('R1:%d'%(N+1))
Lt=L0+sum(Ls[k]*a**(k+1) for k in range(N)); rt=rc+sum(Rs[k]*a**(k+1) for k in range(N))
Lz=Lt*ci; Q=Lt**2*(1-ci**2); Dl=r**2-2*r+a**2
Rpot=(E*(r**2+a**2)-a*Lz)**2-Dl*(r**2+(Lz-a*E)**2+Q)
dR=sp.diff(Rpot,r)
e1=sp.expand(Rpot.subs(r,rt)); e2=sp.expand(dR.subs(r,rt))
sol={}
for k in range(1,N+1):
    c1=sp.together(sp.expand(e1.subs(sol)).coeff(a,k)); c2=sp.together(sp.expand(e2.subs(sol)).coeff(a,k))
    s_=sp.solve([sp.numer(c1),sp.numer(c2)],[Ls[k-1],Rs[k-1]],dict=True)[0]
    sol.update({kk:sp.factor(v) for kk,v in s_.items()})
    print('order',k,'L%d ='%k, sol[Ls[k-1]], flush=True)
Lc=Lt.subs(sol)
L2=sp.expand(sp.series(Lc**2,a,0,N+1).removeO())
sig=sp.integrate(sp.expand(L2.subs(ci,S*sp.cos(phi))),(phi,0,2*sp.pi))/(2*sp.pi*L0**2)
C=sp.symbols('C')
sig=sp.expand(sig)
R_=sp.symbols('R')   # express in r_c afterwards numerically
for k in range(0,N+1,2):
    co=sp.factor(sp.simplify(sig.coeff(a,k).subs(S,sp.sqrt(1-C))))
    print('a^%d:'%k, co, flush=True)
