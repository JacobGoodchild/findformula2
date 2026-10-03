# Capture area to O(a^2) and L_z moment to O(a^3) for any speed, with exact impact-plane mapping
import sympy as sp
a,u,ph,S,C,b,mu,p=sp.symbols('a u phi S C b mu p',positive=True)
rc=(u**2+3)**2/(4*u**2); L0=rc/((3-u**2)/(2*u))
c=sp.symbols('c')
L1=8*c*u**2/((u**2-3)*(u**2+3))
L2=-4*u**3*(c**2*(u**8+28*u**6-234*u**4+252*u**2+81)-16*u**6+288*u**4-144*u**2)/((u**2-3)*(u**2+3)**4*(u**4-18*u**2+9))
L3=128*c*u**6*(c**2*(2*u**8-44*u**6+276*u**4-396*u**2+162)-u**8+36*u**6-342*u**4+324*u**2-81)/((u**2-3)*(u**2+3)**3*(u**4-18*u**2+9)**3)
Lc=lambda m: L0+a*L1.subs(c,m)+a**2*L2.subs(c,m)+a**3*L3.subs(c,m)
p2=4*u**2*(9-u**2)*(u**2-1)/((u**2+3)**2*(u**2-3)**2)
# boundary: p^2 (b^2 - a^2 C) = Lc(mu)^2, mu = b S cos(phi)/sqrt(b^2 - a^2 C)
bs=sp.symbols('b0:4')
bt=bs[0]+a*bs[1]+a**2*bs[2]+a**3*bs[3]
m_=bt*S*sp.cos(ph)/sp.sqrt(bt**2-a**2*C)
eq=sp.series(p2*(bt**2-a**2*C)-Lc(m_)**2,a,0,4).removeO()
sol={}
for k in range(4):
    co=sp.expand(eq).coeff(a,k).subs(sol)
    v=sp.solve(co,bs[k]); v=[x for x in v if (k>0 or x.is_positive is not False)]
    sol[bs[k]]=sp.simplify(v[-1] if k==0 else v[0])
    print('b%d done'%k, flush=True)
bc=bt.subs(sol)
area=sp.integrate(sp.expand(sp.series(bc**2/2,a,0,3).removeO()),(ph,0,2*sp.pi))
A0=sp.simplify(area.coeff(a,0))
print('area a^2 rel:', sp.factor(sp.simplify(area.coeff(a,2)/A0)))
# moment: L_z = p b S cos(phi)  -> Int Int L_z b db dphi = p S Int cos(phi) bc^3/3 dphi
mom=sp.integrate(sp.expand(sp.series(S*sp.cos(ph)*bc**3/3,a,0,4).removeO()),(ph,0,2*sp.pi))
# mean L_z per captured particle = p*mom/area ; per unit energy divide by E
E=(rc-2)/(((u**2+3)/(2*u))*((3-u**2)/(2*u)))
for k in (1,3):
    print('mean Lz/(p) a^%d:'%k, sp.factor(sp.simplify(sp.series(mom/area,a,0,4).removeO().coeff(a,k))))
