# Escape probability for emitter moving with velocity v (along phi) relative to ZAMO, near-horizon extremal,
# using the escape region from [6]: n_phi>0, n_th^2 < (1/k^2-1) n_phi^2, (n_phi>k or n_r>0)
from mpmath import mp, mpf, sqrt, pi, quad, asin, atan, cos, sin
mp.dps=25
def P(k, v):
    g=1/sqrt(1-v*v); K=1/k**2-1
    # emitter-frame direction: n'_phi = t (in [-1,1]), (n'_r, n'_th) = sqrt(1-t^2)(cos psi, sin psi), uniform measure dt dpsi/(4pi)
    def frac(t):
        D=1+v*t; nph=(t+v)/D
        if nph<=0: return mpf(0)
        rho=sqrt(1-t*t)/(g*D)       # |(n_r,n_th)| in ZAMO frame
        # need n_th^2 < K nph^2 : |sin psi| < sqrt(K) nph / rho
        q=sqrt(K)*nph/rho if rho>0 else mpf(10)
        wedge = 2*pi if q>=1 else 4*asin(q)   # measure of psi with |sin psi|<q
        if nph>k: return wedge
        return wedge/2                      # only n_r>0 half (symmetric)
    tk=(k-v)/(1-k*v)        # n_phi>k  <=> t>tk
    t0=-v                   # n_phi>0  <=> t>-v
    pts=[t0,tk,1] if t0<tk<1 else [t0,1]
    return quad(frac,pts)/(4*pi)
if __name__=="__main__":
    print('ZAMO (v=0), k=1/2:', P(mpf(1)/2,0), mpf(7)/24)
    print('ISCO (v=1/2), k=1/2:', P(mpf(1)/2,mpf(1)/2), mpf(5)/12+atan(sqrt(mpf(5)/3))/(sqrt(5)*pi))
