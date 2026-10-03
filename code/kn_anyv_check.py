# Direct check of the KN capture cross-section a^2 coefficient: exact critical impact parameter b(phi)
# from R = R' = 0 in the impact plane, area = (1/2) Int b^2 dphi, then extract the a^2 coefficient.
from mpmath import mp, mpf, sqrt, cos, sin, pi, findroot, acos
mp.dps=40
def b0_rc(r,q): return sqrt(r**4*(r-q**2)/(4*r**2-r**3-4*q**2*r+q**4))
def E_rc(r,q): return (r**2-2*r+q**2)/(r*sqrt(r**2-3*r+2*q**2))
def area(r_c,q,a,th,n=48):
    E=E_rc(r_c,q); p=sqrt(E*E-1); st,ct=sin(th),cos(th)
    tot=0; b0=b0_rc(r_c,q)
    guess=(r_c,b0)
    for k in range(n):
        ph=2*pi*(k+mpf(1)/2)/n
        def eqs(r,b):
            al=b*cos(ph); be=b*sin(ph)
            Lz=p*al*st; Q=p**2*(be**2+(al**2-a**2)*ct**2)
            D=r*r-2*r+a*a+q*q
            R=(E*(r*r+a*a)-a*Lz)**2-D*(r*r+(Lz-a*E)**2+Q)
            dR=4*r*E*(E*(r*r+a*a)-a*Lz)-(2*r-2)*(r*r+(Lz-a*E)**2+Q)-D*2*r
            return [R,dR]
        rr,bb=findroot(eqs,guess)
        tot+=bb**2
    return pi*tot/n   # = (1/2) Int b^2 dphi
def K_formula(r,q,C):
    return (r**2*(r-q**2)-C*(3*r**3-12*r**2-q**2*r**2+18*q**2*r-8*q**4))/(2*r*(r-q**2)*(r**3-6*r**2+9*q**2*r-4*q**4))
for r_c,q,thdeg in [(mpf('3.5'),mpf('0.5'),60),(mpf('3.2'),mpf('0.7'),30),(mpf('3.8'),mpf('0.3'),90),(mpf('3.0'),mpf('0.8'),45)]:
    th=thdeg*pi/180; C=cos(th)**2
    A0=pi*b0_rc(r_c,q)**2
    hs=[mpf('0.004'),mpf('0.008'),mpf('0.012')]
    ys=[(area(r_c,q,h,th)/A0-1)/h**2 for h in hs]
    # Richardson in h^2: y = K + c h^2 + d h^4
    from mpmath import matrix, lu_solve
    M=matrix([[1,h**2,h**4] for h in hs]); K=lu_solve(M,matrix(ys))[0]
    print(r_c,q,thdeg,'numeric K =',mp.nstr(K,15),' formula =',mp.nstr(K_formula(r_c,q,C),15),flush=True)
