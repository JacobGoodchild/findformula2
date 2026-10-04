# Verification of the escape-probability formulas of paperB.tex (Section on extremal horizons).
# The formulas are re-typed from the paper; they are compared with (i) their integral representations
# at high precision and (ii) the brute-force ray simulations in code/ (isotropic grids of emission
# directions at r = M(1 + 1e-7), exact Kerr(-Newman) constants of motion, escape decided from the
# exact radial potential including inward-launched rays that bounce).  Run from the repository root.
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
import numpy as np
from mpmath import mp, mpf, sqrt, pi, asin, atan, quad, nstr
mp.dps = 30
FAIL = []
def check(name, x, y, tol):
    d = abs(mpf(x) - mpf(y)); ok = d < tol
    print(('OK  ' if ok else 'FAIL') + ' %s: %s vs %s (diff %.1e)' % (name, nstr(mpf(x), 12), nstr(mpf(y), 12), float(d)))
    if not ok: FAIL.append(name)

# ---- formulas as in the paper ----
def P_zamo(a, th):                    # photons, emitter at rest (ZAMO), extremal KN, polar angle th
    k = (1 + a**2*mp.cos(th)**2)/(2*a*mp.sin(th))
    return mpf(0) if k >= 1 else mpf(1)/2 - asin(k)/(2*pi) - k/4
def P_isco(a):                        # photons, ISCO emitter, extremal KN, a >= 1/sqrt2
    k = 1/(2*a)
    return mpf(1)/2 - asin(k)/(2*pi) + k*atan(sqrt((1 + k**2)/(1 - k**2)))/(pi*sqrt(1 + k**2))
def eps_zamo(a):                      # energy fraction, emitter at rest, equator
    return (4*a**2 - 1)/(32*a) + sqrt(4*a**2 - 1)/16
def eps_isco_kerr(): return 2*sqrt(3)/9 + 1/(20*pi) + 13*atan(sqrt(mpf(5)/3))/(5*sqrt(15)*pi)
def P_massive(beta, a=1):             # massive, ZAMO emitter, equator, elementary branch
    return (1 - sqrt(1 - beta**2)/(a*beta))/2
def P_massive_isco(beta):             # massive, ISCO emitter of extremal Kerr, elementary branch
    return (2*beta + 1 - sqrt(3*(1 - beta**2)))/(4*beta)

print('== exact values / integral representations ==')
check('equator Kerr = 7/24', P_zamo(mpf(1), pi/2), mpf(7)/24, 1e-25)
check('ISCO Kerr = 5/12 + atan(sqrt(5/3))/(sqrt5 pi)', P_isco(mpf(1)), mpf(5)/12 + atan(sqrt(mpf(5)/3))/(sqrt(5)*pi), 1e-25)
check('ISCO a=1/sqrt2 = 3/8 + sqrt3/9', P_isco(1/sqrt(2)), mpf(3)/8 + sqrt(3)/9, 1e-25)
for k in [mpf(1)/2, mpf(5)/9, mpf(2)/3, mpf('0.7')]:
    integ = mpf(1)/2 + quad(lambda t: asin((t + k)/(k*sqrt(1 - t**2))), [-k, 0])/(2*pi)
    check('ISCO closed form vs integral, k=%s' % nstr(k, 5), P_isco(1/(2*k)), integ, 1e-20)
check('energy fraction ZAMO Kerr = 3/32 + sqrt3/16', eps_zamo(mpf(1)), mpf(3)/32 + sqrt(3)/16, 1e-25)
# energy fraction ZAMO by direct integration over the escape wedge (Kerr, equator): g = n_phi
k = mpf(1)/2
def frac_dir(t):        # fraction of the circle of directions with n_phi = t that escapes (photons, ZAMO)
    if t >= k: return mpf(1)
    if t <= 0: return mpf(0)
    return asin(sqrt((1/k**2 - 1)*t**2/(1 - t**2)))/pi if (1/k**2 - 1)*t**2 < 1 - t**2 else mpf(1)/2
check('photon P (ZAMO) by integration', quad(lambda t: frac_dir(t)/2, [0, k, 1]), mpf(7)/24, 1e-20)
check('energy fraction (ZAMO) by integration', quad(lambda t: t*frac_dir(t)/2, [0, k, 1]), eps_zamo(mpf(1)), 1e-20)
# ISCO energy fraction: emitter moving at v = k; in emitter frame t = cos(angle to motion);
# escapes: all t > 0; fraction (1/pi) asin((t+k)/(k sqrt(1-t^2))) for -k < t < 0; g = a (t + k)/sqrt(1 - k^2)
a = mpf(1); k = 1/(2*a)
def fr(t):
    if t >= 0: return mpf(1)
    if t <= -k: return mpf(0)
    return asin((t + k)/(k*sqrt(1 - t**2)))/pi
check('ISCO photon P by integration', quad(lambda t: fr(t)/2, [-k, 0, 1]), P_isco(a), 1e-20)
check('ISCO energy fraction by integration', quad(lambda t: a*(t + k)/sqrt(1 - k**2)*fr(t)/2, [-k, 0, 1]), eps_isco_kerr(), 1e-20)
check('massive ZAMO exact 1/8 at beta=0.8', P_massive(mpf('0.8')), mpf(1)/8, 1e-25)
check('massive ISCO exact 1/4 at beta=1/2', P_massive_isco(mpf(1)/2), mpf(1)/4, 1e-25)
check('massive ISCO exact 1/2 at beta=sqrt(2/3)', P_massive_isco(sqrt(mpf(2)/3)), mpf(1)/2, 1e-25)
check('massive ISCO threshold (3sqrt2-2)/7 gives P=0', P_massive_isco((3*sqrt(2) - 2)/7), 0, 1e-25)

print('== brute-force ray simulations (grid noise ~1e-3) ==')
import escape_kn2, escape_isco_sim, energy_sim_zamo, energy_sim_isco, massive_sim, massive_kn, isco_massive_sim
for a, thdeg in [(1.0, 90), (0.7, 90), (0.8, 60), (1.0, 70)]:
    th = thdeg*np.pi/180
    check('photon ZAMO sim a=%s %d deg' % (a, thdeg), escape_kn2.P_sim(a, th, 1e-7, N=200), P_zamo(mpf(a), mpf(thdeg)*pi/180), 3e-3)
for a in [1.0, 0.8, 1/np.sqrt(2)]:
    check('photon ISCO sim a=%.4f' % a, escape_isco_sim.P_sim(a, 1e-7, 1/(2*a), N=200), P_isco(mpf(a)), 3e-3)
for a in [1.0, 0.8]:
    check('energy ZAMO sim a=%s' % a, energy_sim_zamo.P_sim(a, np.pi/2, 1e-7, N=200), eps_zamo(mpf(a)), 3e-3)
check('energy ISCO sim a=1', energy_sim_isco.P_sim(1.0, 1e-7, 0.5, N=200), eps_isco_kerr(), 3e-3)
for b in [0.8, 0.75]:
    check('massive ZAMO sim beta=%s' % b, massive_sim.P_sim(b, 1e-7, N=200), P_massive(mpf(b)), 3e-3)
check('massive ZAMO KN sim a=0.8 beta=0.85', massive_kn.P_sim(0.8, 0.85, N=200), P_massive(mpf('0.85'), mpf('0.8')), 3e-3)
for b in [0.5, 0.7]:
    check('massive ISCO sim beta=%s' % b, isco_massive_sim.P_sim(b, N=200), P_massive_isco(mpf(b)), 3e-3)
print('\nSUMMARY: %d failures' % len(FAIL), FAIL)
