"""Quick regression: recompute headline numbers of formulas.txt and compare with independent methods.
Run:  python3 verify_all.py   (takes a few minutes)"""
import numpy as np
from mpmath import mp, mpf, pi, sqrt, asin, atan
def ok(name, x, y, tol):
    d=abs(float(x)-float(y)); print(f"{'OK ' if d<tol else 'FAIL'} {name}: {float(x):.12g} vs {float(y):.12g} (diff {d:.1e})")
import area_closed, verify_shoelace
mp.dps=30
ok('[1] edge-on Kerr a=0.8 closed vs polygon', area_closed.A_closed(mpf('0.8')), verify_shoelace.shoelace(0.8), 1e-6)
ok('[2] a->1 limit', area_closed.A_closed(1-mpf(10)**-20), 16*pi+15*sqrt(3), 1e-6)
import ext_closed, ext_formula
mp.dps=30
th=pi/3
ok('[3] extremal 60deg closed vs direct', ext_closed.A_closed(th), ext_formula.A_direct(th).real, 1e-15)
import master
mp.dps=30
ok('[4] master a=0.9 17deg vs polygon', master.A_master(mpf('0.9'),17*pi/180)[0].real, master.shoelace(0.9,float(17*np.pi/180)), 1e-6)
C=np.cos(np.radians(17))**2; a=0.3
ser=1-(1+C)*a**2/18-(1/72+5*C/972-5*C**2/1944)*a**4
ok('[5] series a=0.3 17deg (to a^4)', ser, master.A_master(mpf(a),17*pi/180)[0].real/(27*np.pi), 2e-5)
import escape_kn2
ok('[6] escape a=0.7 equator, sim vs formula', escape_kn2.P_sim(0.7,np.pi/2,1e-7,N=120), escape_kn2.P_formula(0.7,np.pi/2), 3e-3)
import kn_ext_check
ok('[7] extremal KN a=0.8 edge-on, formula vs polygon', kn_ext_check.A_formula(0.8), kn_ext_check.A_shoelace(0.8), 1e-6)
import boosted
mp.dps=25
ok('[8] boosted cone reproduces known ISCO value', boosted.P(mpf(1)/2,mpf(1)/2), mpf(5)/12+atan(sqrt(mpf(5)/3))/(sqrt(5)*pi), 1e-20)
import kn_master_check as K
mp.dps=30
ok('[12] KN master vs polygon (0.7,45deg,Q^2=0.2)', K.A_kn(0.7,mpf(np.radians(45)),0.2)[0], K.shoelace(0.7,np.radians(45),0.2), 1e-6)
import bv_equator, escape
mp.dps=30
ok('[13] equatorial solid angle z=0.3: 1D formula vs Phi-integral', bv_equator.Omega_eq(mpf('0.3'))[0], escape.Omega(mpf('0.3'),pi/2), 1e-20)
