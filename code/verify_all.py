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
# ---- session 2 entries ----
from mpmath import ellipe, ellipk
import mb_area, mb_knext, bv_kn, massive_escape, anyv_ext_area
mp.dps=30
a_=mpf('0.5')
ok('[15] slow capture edge-on closed vs curve integral', 7*pi+pi*sqrt(1-a_**2)+16*sqrt(1+a_)*ellipe(2*a_/(1+a_)), mb_area.A_mb_edge('0.5'), 1e-20)
ok('[15] a=1 value', 7*pi+16*sqrt(2), mb_knext.A_formula(1)[0], 1e-15)
ok('[14] extremal KN solid angle far-field -> [7](1-2z)', bv_kn.Omega_kn('1e-4','0.8')/mpf('1e-8'), kn_ext_check.A_formula(0.8)*(1-2e-4), 1e-3)
ok('[19] massive escape beta=0.8 is 1/8', massive_escape.P('0.8'), mpf(1)/8, 1e-20)
ok('[19] massive escape beta=1 is 7/24', massive_escape.P('1').real, mpf(7)/24, 1e-20)
ok('[21] extremal any-speed capture: slow limit', anyv_ext_area.area('0.0001')[0]*mpf('0.0001')**2, 7*pi+16*sqrt(2), 1e-3)
ok('[22] extremal KN slow capture a->0 -> pi phi^5', mb_knext.A_formula('0.001')[0], pi*((1+sqrt(5))/2)**5, 1e-3)
import mb_master, isco_massive, kerrsen_check, ks_slow_sim, perimeter, kerrsen_incl_check
mp.dps=30
import mb_incl
ok('[23] one-integral slow capture any direction vs [16] integral (a=0.9, 30deg)', mb_master.sigma_master('0.9',mpf(30)*pi/180), mb_incl.A_mb('0.9',mpf(30)*pi/180), 1e-12)
ok('[24] ISCO massive escape at beta=1/2 is 1/4', isco_massive.P('0.5').real, mpf(1)/4, 1e-15)
ok('[25] Kerr-Sen extremal a=1 -> 16pi+15sqrt3', kerrsen_check.A_ext(1), 16*pi+15*sqrt(3), 1e-20)
ok('[29] Kerr-Sen slow capture b=0 -> [15]', ks_slow_sim.A_formula(0.5,0.0), 7*pi+pi*sqrt(1-mpf('0.25'))+16*sqrt(mpf('1.5'))*ellipe(mpf(2)/3), 1e-10)
ok('[30] perimeter a->1 -> 18 sqrt3', perimeter.P_closed('0.99999999'), 18*sqrt(3), 2e-3)
ok('[30] perimeter a->0 -> 2 pi sqrt27', perimeter.P_closed('0.0000001'), 2*pi*sqrt(27), 1e-6)

# [17](b') torque a^3 coefficient vs stored high-precision values (code/torque_a3_c4.out)
def _t3(r,C): return (1-C)*(1.5*(r-5)-C*(7*r*r-63*r+144)/(2*r))/((r-2)*(6-r)**3)
ok("[17b'] r=7/2 C=1/4", _t3(3.5,0.25), -0.0825714285714285714, 1e-14)
ok("[17b'] r=10/3 C=4/5", _t3(10/3,0.8), -0.030955078125, 1e-14)
ok("[17b'] r=13/4 C=1/2", _t3(3.25,0.5), -0.069999422065537768017, 1e-14)

# [31] KN capture a^2 coefficient vs direct numerics (code/kn_anyv_check.out)
def _K31(r,q,C): return (r**2*(r-q**2)-C*(3*r**3-12*r**2-q**2*r**2+18*q**2*r-8*q**4))/(2*r*(r-q**2)*(r**3-6*r**2+9*q**2*r-4*q**4))
import math
ok("[31] KN r=3.5 q=0.5 60deg", _K31(3.5,0.5,0.25), -0.0790432393693279, 1e-12)
ok("[31] KN r=3.0 q=0.8 45deg", _K31(3.0,0.8,0.5), -0.136651895747247, 1e-12)
