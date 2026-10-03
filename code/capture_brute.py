# Brute-force capture cross-section (impact-plane polar grid + bisection) vs the O(a^4) formula of [17]
import numpy as np
def sigma_brute(a,E,th,nphi=720):
    p2=E*E-1; p=np.sqrt(p2); s=np.sin(th); c=np.cos(th)
    rp=1+np.sqrt(1-a*a); rr=rp+np.geomspace(1e-10,200,8000); D=rr**2-2*rr+a*a
    def cap(b,phi):
        x=b*np.cos(phi); y=b*np.sin(phi)
        Lz=-p*s*x; Q=p2*y*y-a*a*p2*c*c+(Lz*c/s)**2 if s>0 else p2*y*y
        R=(E*(rr**2+a*a)-a*Lz)**2-D*(rr**2+(Lz-a*E)**2+Q)
        return np.all(R>0)
    area=0.0
    phis=(np.arange(nphi)+0.5)*2*np.pi/nphi
    for ph in phis:
        lo,hi=0.0,30.0
        for _ in range(55):
            m=(lo+hi)/2
            if cap(m,ph): lo=m
            else: hi=m
        area+=lo*lo/2*(2*np.pi/nphi)
    return area
def sigma_series(a,E,th):
    r=None
    E2=E*E; disc=(3*E2-4)**2-16*(1-E2)
    for sg in (1,-1):
        x=(-(3*E2-4)+sg*np.sqrt(disc))/(2*(1-E2))
        if 3<x<=4: r=x
    C=np.cos(th)**2; p2=E2-1
    sS=np.pi*r*r/((r-3)*p2)
    A8=-r**2*(2*r**3+42*r**2-549*r+1458); B8=6*r*(2*r**4+2*r**3-239*r**2+1302*r-2016)
    G8=-10*r**5-2*r**4+1653*r**3-13266*r**2+39744*r-41472
    f=1-a*a*(r+3*C*(4-r))/(2*r*r*(6-r))+a**4*(A8+B8*C+G8*C*C)/(8*r**4*(6-r)**5)
    return sS*f, sS*(1-a*a*(r+3*C*(4-r))/(2*r*r*(6-r))), sS
for a,E,deg in [(0.3,1.25,60),(0.3,2.0,30),(0.2,1.25,90)]:
    th=np.radians(deg); b=sigma_brute(a,E,th); s4,s2,s0=sigma_series(a,E,th)
    print(a,E,deg,'brute',b,'series a^4',s4,'(a^2 only',s2,', a=0',s0,')',flush=True)
