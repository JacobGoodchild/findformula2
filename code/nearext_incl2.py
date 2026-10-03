from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
import ext_closed
from nearext_incl import A_master_hp
mp.dps=50
for deg in [50,85,60]:
    th=mpf(deg)*pi/180; A1=ext_closed.A_closed(th); c2=cos(th)**2; s=sin(th)
    pred=pi*(3-6*c2-c2*c2)/(sqrt(2)*s**3)
    v=[(A_master_hp(1-e,th)-A1)/sqrt(e) for e in (mpf(10)**-10,mpf(10)**-12)]
    ext=v[1]-(v[0]-v[1])/(10**(0.9)-1)
    print(deg, [nstr(x,10) for x in v], 'extrap',nstr(ext,8),'pred',nstr(pred,10))
