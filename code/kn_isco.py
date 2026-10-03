# Extremal Kerr-Newman (M=1, Q^2=1-a^2): prograde equatorial circular orbits, ISCO and ZAMO-relative speed
from mpmath import mp, mpf, sqrt, findroot, diff
mp.dps=30
def orbit(r,a):
    Q2=1-a*a; D=r*r-2*r+a*a+Q2
    gtt=-(1-(2*r-Q2)/r**2); gtp=-a*(2*r-Q2)/r**2; gpp=r*r+a*a+a*a*(2*r-Q2)/r**2
    s=sqrt(r-Q2); Om=s/(r*r+a*s)
    ut=1/sqrt(-(gtt+2*gtp*Om+gpp*Om*Om))
    E=-(gtt+gtp*Om)*ut
    A=(r*r+a*a)**2-D*a*a; om=a*(r*r+a*a-D)/A; al=sqrt(r*r*D/A); sg=sqrt(A)/r
    v=(Om-om)*sg/al
    return E,v
def isco(a):
    a=mpf(a)
    dE=lambda r: diff(lambda x: orbit(x,a)[0], r)
    # scan for sign change of dE above horizon
    rs=[1+mpf(10)**(-k) for k in range(8,0,-1)]+[1+mpf(i)/10 for i in range(2,60)]
    prev=None
    for r in rs:
        try: val=dE(r)
        except Exception: continue
        if abs(val.imag)>1e-20 if hasattr(val,'imag') else False: prev=None; continue
        val=val.real
        if prev is not None and prev[1]<0<val: return findroot(dE,(prev[0],r),solver='bisect')
        prev=(r,val)
    return None
for a in ['0.3','0.5','0.55','0.577350269','0.6','0.7','0.9','1']:
    r=isco(a); a_=mpf(a)
    print(a, r, [orbit(1+mpf(10)**-k,a_)[1] for k in (4,6,8)], [orbit(1+mpf(10)**-k,a_)[0] for k in (4,8)])
