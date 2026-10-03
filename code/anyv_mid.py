from mpmath import mp, mpf, pi, sqrt, nstr, acos
from anyv_f2 import fcoef
mp.dps=40
for r,C in [(mpf(7)/2,mpf(1)/4),(mpf(13)/4,mpf(3)/4),(mpf(15)/4,mpf(1)/2)]:
    E=sqrt((r-2)**2/(r*(r-3))); th=acos(sqrt(C))
    f=fcoef(E,th); pred=(r+3*C*(4-r))/(2*r*r*(6-r))
    print(nstr(r,4),nstr(C,3), nstr(f,15), nstr(pred,15))
