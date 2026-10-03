from mpmath import mp, mpf, sqrt, pi, quad, asin, atan, acos, log, pslq, identify
from boosted import P
mp.dps=30
def P_isco_int(k):
    J=quad(lambda t: asin((t+k)/(k*sqrt(1-t*t))),[-k,0])
    return mpf(1)/2+J/(2*pi)
for k in [mpf(1)/2, mpf(5)/9, mpf(2)/3]:
    print(k, P_isco_int(k), P(k,k))
k=mpf(1)/2
print(P_isco_int(k), mpf(5)/12+atan(sqrt(mpf(5)/3))/(sqrt(5)*pi))
