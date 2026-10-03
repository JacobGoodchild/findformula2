# Analytic perturbation: critical angular momentum L_c(i) for Kerr, energy E (parametrised by Schwarzschild r_c), to O(a^4)
import sympy as sp
a,ci,rc,phi,S=sp.symbols('a c_i r_c phi S')
r=sp.symbols('r')
E2=(rc-2)**2/(rc*(rc-3)); E=sp.sqrt(E2)
N=4
Ls=sp.symbols('L1:%d'%(N+1)); Rs=sp.symbols('R1:%d'%(N+1))
L0=sp.sqrt(rc**2/(rc-3))
Lt=L0+sum(Ls[k]*a**(k+1) for k in range(N))
rt=rc+sum(Rs[k]*a**(k+1) for k in range(N))
Lz=Lt*ci; Q=Lt**2*(1-ci**2)
Dl=r**2-2*r+a**2
Rpot=(E*(r**2+a**2)-a*Lz)**2-Dl*(r**2+(Lz-a*E)**2+Q)
dR=sp.diff(Rpot,r)
e1=sp.expand(Rpot.subs(r,rt)); e2=sp.expand(dR.subs(r,rt))
sol={}
for k in range(1,N+1):
    c1=sp.simplify(sp.series(e1.subs(sol),a,0,k+1).removeO().coeff(a,k))
    c2=sp.simplify(sp.series(e2.subs(sol),a,0,k+1).removeO().coeff(a,k))
    s_=sp.solve([c1,c2],[Ls[k-1],Rs[k-1]],dict=True)[0]
    sol.update({kk:sp.simplify(v) for kk,v in s_.items()})
    print('order',k,'L%d ='%k, sp.factor(sol[Ls[k-1]]), flush=True)
Lc=Lt.subs(sol)
L2=sp.expand(sp.series(Lc**2,a,0,N+1).removeO())
sig=sp.integrate(sp.expand(L2.subs(ci,S*sp.cos(phi))),(phi,0,2*sp.pi))/(2*sp.pi*L0**2)
C=sp.symbols('C')
sig=sp.expand(sp.simplify(sig))
for k in range(0,N+1,2):
    co=sp.factor(sp.simplify(sig.coeff(a,k).subs(S,sp.sqrt(1-C))))
    print('a^%d:'%k, co)
