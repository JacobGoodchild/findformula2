# Massive particles launched isotropically (speed beta in the orbiting frame) from the near-horizon ISCO of extremal Kerr
from mpmath import mp, mpf, sqrt, pi, quad, asin
mp.dps=25
def P(beta):
    beta=mpf(beta); gb=1/sqrt(1-beta**2); gv=2/sqrt(3); v=mpf(1)/2
    # emitter frame: n_phi = t, (n_r, n_theta) = sqrt(1-t^2)(cos psi, sin psi)
    def frac(t):
        E=gv*gb*(beta*t+v)
        if E<1: return mpf(0)
        K=3*E*E-1                         # need gb^2 beta^2 n_theta^2 < K
        rho=sqrt(1-t*t)
        q=sqrt(K)/(gb*beta*rho) if rho>0 else mpf(10)
        w=2*pi if q>=1 else 4*asin(q)
        return w if t>0 else w/2
    t0=(sqrt(3)/2*sqrt(1-beta**2)-mpf(1)/2)/beta
    lo=max(t0,mpf(-1))
    if lo>=1: return mpf(0)
    pts=[lo,1] if lo>=0 else [lo,0,1]
    return quad(frac,pts)/(4*pi)
if __name__=="__main__":
    print('photon-limit check', P('0.9999999'), mpf(5)/12+__import__('mpmath').atan(sqrt(mpf(5)/3))/(sqrt(5)*pi))
    for b in ['0.32','0.4','0.5','0.6','0.7','0.8','0.9']:
        print(b, P(b))
    print('threshold', (3*sqrt(2)-2)/7)
