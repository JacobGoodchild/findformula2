# Polar (C = 1) capture cross-section, any speed, universal Delta = r^2 f + a^2, to order a^4 (u = a^2).
# Exact on the axis: sigma = pi (Q/p^2 + a^2) with R = R' = 0 at L_z = 0.
import sympy as sp
u,R,x1,x2,Q0,Q1,Q2=sp.symbols('u R x1 x2 Q0 Q1 Q2')
fs=sp.symbols('f0:7'); r=sp.Symbol('r')
F=sum(fs[k]*(r-R)**k/sp.factorial(k) for k in range(7))
E2=2*fs[0]**2/(2*fs[0]-R*fs[1]); L2=fs[1]*R**3/(2*fs[0]-R*fs[1])   # static critical orbit
Rp=E2*(r**2+u)**2-(r**2*F+u)*(r**2+u*E2+sp.Symbol('Q'))
dR=sp.diff(Rp,r)
rs=R+x1*u+x2*u**2; Qs=L2+Q1*u+Q2*u**2
e1=sp.expand(Rp.subs(sp.Symbol('Q'),Qs).subs(r,rs))
e2=sp.expand(dR.subs(sp.Symbol('Q'),Qs).subs(r,rs))
print('order0:',sp.simplify(e1.coeff(u,0)),sp.simplify(e2.coeff(u,0)),flush=True)
s1=sp.solve([e1.coeff(u,1),e2.coeff(u,1)],[x1,Q1],dict=True)[0]
s2=sp.solve([sp.expand(e1.coeff(u,2).subs(s1)),sp.expand(e2.coeff(u,2).subs(s1))],[x2,Q2],dict=True)[0]
p2=E2-1
b2=(Qs/p2+u).subs(s1).subs(s2)
c2=sp.factor(sp.simplify(sp.expand(b2).coeff(u,1)/(L2/p2)))
c4=sp.factor(sp.simplify(sp.expand(b2).coeff(u,2)/(L2/p2)))
print('a^2:',c2); print('a^4:',c4)
rr=sp.Symbol('r')
fK=1-2/rr
sub={fs[k]:sp.diff(fK,rr,k).subs(rr,R) for k in range(7)}
print('Kerr a^2:',sp.factor(c2.subs(sub)),' a^4:',sp.factor(c4.subs(sub)))
open('univ_polar.out','w').write('a^2: %s\na^4: %s\n'%(c2,c4))
