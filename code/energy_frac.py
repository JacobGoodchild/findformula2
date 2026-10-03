# Energy-weighted escape fraction (photon energy at infinity / emitted energy, averaged over isotropic emission)
# near-horizon extremal KN, equator; ZAMO lamp (v=0): E_inf/E_emit = a n_phi ; ISCO lamp (v=k): a (t+k)/sqrt(1-k^2)
from mpmath import mp, mpf, sqrt, pi, quad, asin, acos, pslq, identify, atan, log
mp.dps=30
def eps_zamo(a):
    k=1/(2*a); gam=acos(k)
    def inner(psi):
        m=k/cos(psi) if False else k/mp.cos(psi)
        # n_phi = sin(chi) cos(psi); weight a*n_phi*sin(chi) dchi
        f=lambda chi: a*mp.sin(chi)**2*mp.cos(psi)
        full=quad(f,[0,pi])
        if m>=1: return full/2
        c0=asin(m)
        mid=quad(f,[c0,pi-c0])
        return mid+(full-mid)/2
    return quad(inner,[-gam,0,gam])/(4*pi)
def eps_isco(a):
    k=1/(2*a); w=lambda t: a*(t+k)/sqrt(1-k*k)
    J1=quad(lambda t: w(t),[0,1])*2*pi
    J2=quad(lambda t: 2*asin((t+k)/(k*sqrt(1-t*t)))*w(t),[-k,0])
    return (J1+J2)/(4*pi)
for a in [mpf(1)]:
    ez=eps_zamo(a); ei=eps_isco(a)
    print('zamo',ez, identify(ez,['pi','sqrt(3)']), pslq([ez,1,sqrt(3)/pi,1/pi,sqrt(3)],maxcoeff=10**5,maxsteps=10**6))
    print('isco',ei, identify(ei,['pi','sqrt(3)','sqrt(5)']), pslq([ei,1,sqrt(3),1/pi,sqrt(3)/pi,sqrt(5)/pi,atan(sqrt(mpf(5)/3))/pi,sqrt(5)*atan(sqrt(mpf(5)/3))/pi,sqrt(3)*atan(sqrt(mpf(5)/3))/pi,sqrt(15)*atan(sqrt(mpf(5)/3))/pi],maxcoeff=10**5,maxsteps=10**6))

def eps_isco_closed(a):
    k=1/(2*a); b=1+k*k; ph=asin(sqrt(b/2)); ac=acos(k)
    S=(2*k*k/b**1.5)*ph - k*k*sqrt(1-k**4)/b**1.5 - (2*(1+2*k*k)/sqrt(b))*ph + (1-k)**2*ac/2 + (1+k)**2*(pi/2-ac/2)
    K2=k*k*pi/4-S/2
    return a*(mpf(1)/2+k)/(2*sqrt(1-k*k)) + a*K2/(2*pi*sqrt(1-k*k))
for a in [mpf(1),mpf('0.8'),1/sqrt(2)]:
    print('check', a, eps_isco(a), eps_isco_closed(a))
