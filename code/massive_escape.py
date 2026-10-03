# Escape probability for massive particles (local speed beta, isotropic in ZAMO frame) from the near-horizon of extremal Kerr (equator)
from mpmath import mp, mpf, sqrt, pi, quad, asin
mp.dps=30
def P(beta):
    beta=mpf(beta)
    if beta**2<mpf(1)/2: return mpf(0)
    t0=sqrt(1-beta**2)/beta; t1=1/(2*beta)
    if t0>=t1: return (1-t0)/2
    J=quad(lambda t: asin(sqrt((3*t*t+1-1/beta**2)/(1-t*t))),[t0,t1])
    return (1-t1)/2+J/(2*pi)
if __name__=="__main__":
    for b in ['1','0.95','0.9','0.866','0.8','0.75','0.7071068']:
        print(b, P(b))
