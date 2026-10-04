# Verification of the displayed formulas of paperB.tex (Kerr shadow areas/perimeter, slow capture).
# Formulas are re-typed from the paper and compared with direct numerical computations:
#   - direct quadrature of the exact parametric shadow edge / capture boundary,
#   - polygon (shoelace) areas, and brute-force bisection with the exact radial potential.
# Run from the repository root:  python3 papers/paperB/verify_paperB.py
import sympy as sp
from mpmath import (mp, mpf, sqrt, pi, quad, ellipk, ellipe, ellippi, cos, sin, acos, asin, log,
                    findroot, nstr, diff, polyroots, ellipf)
mp.dps = 40
FAIL = []
def check(name, x, y, tol):
    from mpmath import re as _re, im as _im
    if abs(_im(x)) < mpf(10)**-15: x = _re(x)
    if abs(_im(y)) < mpf(10)**-15: y = _re(y)
    d = abs(x - y); ok = d < tol
    print(('OK  ' if ok else 'FAIL') + ' %s: %s vs %s (diff %.1e)' % (name, nstr(x, 18), nstr(y, 18), float(d)))
    if not ok: FAIL.append(name)
def scheck(name, e):
    z = sp.simplify(e); ok = (z == 0)
    print(('OK  ' if ok else 'FAIL') + ' %s (symbolic difference = %s)' % (name, z))
    if not ok: FAIL.append(name)

# ---------------- Kerr photon-orbit data (direct) ----------------
def xi(r, a): return (r**2*(3 - r) - a**2*(r + 1))/(a*(r - 1))
def eta(r, a): return r**3*(4*a**2 - r*(r - 3)**2)/(a**2*(r - 1)**2)
def dxi(r, a): return diff(lambda x: xi(x, a), r)
def photon_roots(a):
    rp = 2*(1 + cos(mpf(2)/3*acos(-a))); rr = 2*(1 + cos(mpf(2)/3*acos(a)))
    return rp, rr
def direct_area(a, th):
    s, c = sin(th), cos(th)
    b2 = lambda r: eta(r, a) + a*a*c*c - xi(r, a)**2*c*c/(s*s)
    rp, rr = photon_roots(a)
    r1 = findroot(b2, rp + mpf('1e-6')) if th != pi/2 else rp
    r2 = findroot(b2, rr - mpf('1e-6')) if th != pi/2 else rr
    g = lambda t: (lambda r: sqrt(max(b2(r), 0))*abs(dxi(r, a))/s*(r2 - r1)*sin(t)/2)(r1 + (r2 - r1)*(1 - cos(t))/2)
    return 2*quad(g, [0, pi/2, pi])
def polygon_area(a, th, n=20000):
    s, c = sin(th), cos(th)
    b2 = lambda r: eta(r, a) + a*a*c*c - xi(r, a)**2*c*c/(s*s)
    rp, rr = photon_roots(a)
    r1 = findroot(b2, rp + mpf('1e-6')) if th != pi/2 else rp
    r2 = findroot(b2, rr - mpf('1e-6')) if th != pi/2 else rr
    pts = []
    for k in range(n + 1):
        r = r1 + (r2 - r1)*(1 - cos(pi*k/n))/2
        pts.append((-xi(r, a)/s, sqrt(max(b2(r), 0))))
    up = sum(pts[i][0]*pts[i + 1][1] - pts[i + 1][0]*pts[i][1] for i in range(n))
    return abs(up)  # upper half closed by the alpha-axis, doubled: 2 * |up|/2

print('== edge-on Kerr shadow area ==')
def area_edge_closed(a):        # paper eq. (edge closed form)
    rp, rr = photon_roots(a); r0 = 6 - rp - rr
    m = (rr - rp)*r0/((rr - r0)*rp); n = (rr - rp)/(rr - r0); n2 = n*(r0 - 1)/(rp - 1); g = 2/sqrt((rr - r0)*rp)
    K, E, P = ellipk(m), ellipe(m), ellippi(n, m)
    V = (n*E + (m - n)*K + (2*n*m + 2*n - n**2 - 3*m)*P)/(2*(n - 1)*(m - n))
    I0 = g*K; I1 = g*(r0*K + (rp - r0)*P); I2 = g*(r0**2*K + 2*r0*(rp - r0)*P + (rp - r0)**2*V)
    J = g*(K/(r0 - 1) + (r0 - rp)/((rp - 1)*(r0 - 1))*ellippi(n2, m))
    return 6*I0 - 21*I1 + 15*I2 + 6*(1 - a**2)*J
def area_edge_integral(a):      # paper eq. (one-line integral)
    rp, rr = photon_roots(a)
    f = lambda t: (lambda r: 3*((r - 1)*(5*r - 2) + 2*(1 - a**2)/(r - 1))/sqrt(r*(4*a**2 - r*(r - 3)**2))*(rr - rp)*sin(t)/2)(rp + (rr - rp)*(1 - cos(t))/2)
    return quad(f, [0, pi/2, pi])
for a in ['0.001', '0.3', '0.5', '0.9', '0.998']:
    a = mpf(a)
    cf = area_edge_closed(a)
    check('edge closed vs one-line integral a=%s' % a, cf, area_edge_integral(a), mpf('1e-15'))
    check('edge closed vs direct contour a=%s' % a, cf, direct_area(a, pi/2), mpf('1e-15'))
for a in ['0.2', '0.8']:
    a = mpf(a); mp.dps = 20
    check('edge closed vs polygon (20000 pts) a=%s' % a, area_edge_closed(a), polygon_area(a, pi/2), mpf('1e-5'))
    mp.dps = 40
check('a->0 limit 27pi', area_edge_closed(mpf('1e-12')), 27*pi, mpf('1e-15'))
# series
a = mpf('0.05')
ser = 27*pi*(1 - a**2/18 - a**4/72 - 2059*a**6/314928 - 88175*a**8/22674816 - 355999*a**10/136048896)
check('small-spin series (a=0.05, through a^10)', area_edge_closed(a), ser, mpf('1e-14'))
# near-extremal
mp.dps = 80
eps = mpf('1e-40'); A1 = 16*pi + 15*sqrt(3)
Aeps = area_edge_closed(1 - eps)
check('near-extremal sqrt(eps) coefficient', (Aeps - A1)/sqrt(eps), 3*pi/sqrt(2), mpf('1e-15'))
b = (Aeps - A1 - 3*pi/sqrt(2)*sqrt(eps))/eps
check('near-extremal eps coefficient', b, 8/sqrt(3)*log(27/(2*eps)) - 6*sqrt(3), mpf('1e-15'))
mp.dps = 40
check('a=1 limit 16pi+15sqrt3', area_edge_closed(1 - mpf('1e-30')), A1, mpf('1e-12'))

print('== extremal Kerr, any inclination ==')
r, s, c = sp.symbols('r s c', positive=True)
lhs = s**2*(r**3*(4 - r) + c**2 - (c**2/s**2)*(r**2 - 2*r - 1)**2)
rhs = -(r**2 - 2*(1 + s)*r - c**2)*(r**2 - 2*(1 - s)*r - c**2)
scheck('factorisation (with c^2 = 1 - s^2)', sp.expand((lhs - rhs).subs(c**2, 1 - s**2)))
def ext_area(th):               # paper: 4 Int (3r^2 - c^2)/sqrt(R) + B
    s, c = sin(th), cos(th)
    rA = 1 - s - sqrt(2 - 2*s); rB = 1 + s - sqrt(2 + 2*s); rC = 1 - s + sqrt(2 - 2*s); rD = 1 + s + sqrt(2 + 2*s)
    thc = asin(sqrt(3) - 1)
    if th <= thc: lo, B = rC, 0
    else: lo, B = mpf(1), (2 + s**2)*sqrt(3*s**2 + s**2*c**2 - 4*c**2)/s**2
    Rf = lambda r: (r - rA)*(r - rB)*(r - rC)*(rD - r)
    f = lambda t: (lambda r: 4*(3*r**2 - c**2)/sqrt(Rf(r))*(rD - lo)*sin(t)/2 if Rf(r) > 0 else mpf(0))(lo + (rD - lo)*(1 - cos(t))/2)
    return quad(f, [0, pi/2, pi]) + B
def ext_closed(th):              # closed form for theta <= theta_c
    s, c = sin(th), cos(th)
    d_ = 1 - s - sqrt(2 - 2*s); cc = 1 + s - sqrt(2 + 2*s); b_ = 1 - s + sqrt(2 - 2*s); a_ = 1 + s + sqrt(2 + 2*s)
    m = (a_ - b_)*(cc - d_)/((a_ - cc)*(b_ - d_)); n = (a_ - b_)/(a_ - cc); g = 2/sqrt((a_ - cc)*(b_ - d_))
    K, E, P = ellipk(m), ellipe(m), ellippi(n, m)
    V = (n*E + (m - n)*K + (2*n*m + 2*n - n**2 - 3*m)*P)/(2*(n - 1)*(m - n))
    I0 = g*K; I2 = g*(cc**2*K + 2*cc*(b_ - cc)*P + (b_ - cc)**2*V)
    return 4*(3*I2 - c**2*I0)
def direct_ext(th):
    s, c = sin(th), cos(th)
    b2 = lambda r: r**3*(4 - r) + c*c - (1 + 2*r - r*r)**2*c*c/(s*s)
    rD = 1 + s + sqrt(2 + 2*s); thc = asin(sqrt(3) - 1)
    lo = 1 - s + sqrt(2 - 2*s) if th <= thc else mpf(1)
    g = lambda t: (lambda r: sqrt(max(b2(r), 0))*abs(2 - 2*r)/s*(rD - lo)*sin(t)/2)(lo + (rD - lo)*(1 - cos(t))/2)
    A = 2*quad(g, [0, pi/2, pi])
    if th > thc:   # straight NHEK segment at alpha = -xi(1)/s = -2/s, |beta| <= beta_m
        A += 0     # its contribution is the polygon between r=1 endpoint and the segment; handled by B in formula
    return A
for deg in [10, 30, 45]:
    th = deg*pi/180
    check('extremal closed vs integral %d deg' % deg, ext_closed(th), ext_area(th), mpf('1e-15'))
    check('extremal integral vs direct contour %d deg' % deg, ext_area(th), direct_ext(th), mpf('1e-15'))
check('extremal theta=0 limit pi(12+8sqrt2)', ext_closed(mpf('1e-8')), pi*(12 + 8*sqrt(2)), mpf('1e-12'))
check('extremal theta=90 deg = 16pi+15sqrt3', ext_area(pi/2), 16*pi + 15*sqrt(3), mpf('1e-15'))
check('theta_c', asin(sqrt(3) - 1)*180/pi, mpf('47.0585971351200897'), mpf('1e-14'))
s_ = mpf('0.2')
ser = pi*(12 + 8*sqrt(2)) + sqrt(2)*pi*(s_**2/2 + s_**4/128 - 7*s_**6/2048 - 1735*s_**8/524288)
check('near-axis series (s=0.2)', ext_closed(asin(s_)), ser, mpf('1e-6'))

print('== master formula (any a, theta) ==')
def master(a, th):
    u = a**2*cos(th)**2 if th != pi/2 else mpf(0)
    N = lambda r: -r**6 + 6*r**5 - (9 + 2*u)*r**4 + 4*a**2*r**3 + 6*u*r**2 - 4*a**2*u*r - u**2*(r - 1)**2
    rp, rr = photon_roots(a)
    r1 = findroot(N, rp + mpf('1e-6')); r2 = findroot(N, rr - mpf('1e-6'))
    P = lambda r: (2*(u**2 - 5*a**2*u + 6*u + 3*a**2 - 3) - 2*(2*u**2 + a**2*u - 3*u - 3*a**2 + 3)*r
                   + 2*(u**2 - 6*u - 3)*r**2 + 4*(u + 3)*r**3 + 2*(u - 3)*r**4 - 2*(a**2 - 1)*(u**2 + 6*u - 3)/(r - 1))
    f = lambda t: (lambda r: P(r)/sqrt(N(r))*(r2 - r1)*sin(t)/2 if N(r) > 0 else mpf(0))(r1 + (r2 - r1)*(1 - cos(t))/2)
    return quad(f, [0, pi/2, pi])/(u - 1)
for a, deg in [('0.5', 17), ('0.9', 17), ('0.5', 60), ('0.94', 30), ('0.99', 85), ('0.998', 45)]:
    a = mpf(a); th = deg*pi/180
    check('master vs direct a=%s %d deg' % (a, deg), master(a, th), direct_area(a, th), mpf('1e-15'))
check('master at 90 deg = edge closed form', master(mpf('0.7'), pi/2), area_edge_closed(mpf('0.7')), mpf('1e-15'))
for deg in [60, 80]:
    th = deg*pi/180
    # near-extremal check with the independent master formula: A(1-eps) = A_ext + c(theta) sqrt(eps) + O(eps)
    eps = mpf('1e-12'); s_, c_ = sin(th), cos(th)
    ctheta = pi*(3 - 6*c_**2 - c_**4)/(sqrt(2)*s_**3)
    check('extremal integral + c(theta) sqrt(eps) vs master at a=1-1e-12, %d deg' % deg, ext_area(th) + ctheta*sqrt(eps), master(1 - eps, th), mpf('1e-9'))

print('== perimeter ==')
def perim_edge_closed(a):
    rp, rr = photon_roots(a); r0 = 6 - rp - rr
    m = (rr - rp)/(rr - r0); n = (rr - rp)/(rr - 1)
    K, E, P = ellipk(m), ellipe(m), ellippi(n, m)
    V = (n*E + (m - n)*K + (2*n*m + 2*n - n**2 - 3*m)*P)/(2*(n - 1)*(m - n))
    return 16/sqrt(rr - r0)*(r0*K + (rr - r0)*E - K + (1 - a**2)*V/(rr - 1)**2)
def perim_edge_integral(a):
    rp, rr = photon_roots(a)
    f = lambda t: (lambda r: 8*((r - 1) + (1 - a**2)/(r - 1)**2)/sqrt(4*a**2 - r*(r - 3)**2)*(rr - rp)*sin(t)/2)(rp + (rr - rp)*(1 - cos(t))/2)
    return quad(f, [0, pi/2, pi])
def perim_direct(a, th):
    s, c = sin(th), cos(th)
    b2 = lambda r: eta(r, a) + a*a*c*c - xi(r, a)**2*c*c/(s*s)
    rp, rr = photon_roots(a)
    r1 = findroot(b2, rp + mpf('1e-6')) if th != pi/2 else rp
    r2 = findroot(b2, rr - mpf('1e-6')) if th != pi/2 else rr
    al = lambda r: -xi(r, a)/s; be = lambda r: sqrt(max(b2(r), 0))
    g = lambda t: (lambda r: sqrt(diff(al, r)**2 + diff(lambda x: be(x)**2, r)**2/(4*max(b2(r), mpf(10)**-60)))*(r2 - r1)*sin(t)/2)(r1 + (r2 - r1)*(1 - cos(t))/2)
    return 2*quad(g, [0, pi/2, pi])
def perim_general(a, th):       # paper: one-integral formula for any (a, theta)
    u = a**2*cos(th)**2
    N = lambda r: -r**6 + 6*r**5 - (9 + 2*u)*r**4 + 4*a**2*r**3 + 6*u*r**2 - 4*a**2*u*r - u**2*(r - 1)**2
    rp, rr = photon_roots(a)
    r1 = findroot(N, rp + mpf('1e-6')) if th != pi/2 else rp; r2 = findroot(N, rr - mpf('1e-6')) if th != pi/2 else rr
    f = lambda t: (lambda r: 8*((r - 1)**3 + 1 - a**2)*sqrt(r*(r**2 - u))/((r - 1)**2*sqrt(N(r)))*(r2 - r1)*sin(t)/2)(r1 + (r2 - r1)*(1 - cos(t))/2)
    return quad(f, [0, pi/2, pi])
for a in ['0.3', '0.7', '0.9', '0.999']:
    a = mpf(a)
    check('perimeter closed vs integral a=%s' % a, perim_edge_closed(a), perim_edge_integral(a), mpf('1e-15'))
    check('perimeter general formula at 90 deg a=%s' % a, perim_general(a, pi/2), perim_edge_closed(a), mpf('1e-15'))
mp.dps = 25
for a, deg in [('0.5', 60), ('0.9', 17)]:
    a = mpf(a); th = deg*pi/180
    check('perimeter general vs direct arc length a=%s %d deg' % (a, deg), perim_general(a, th), perim_direct(a, th), mpf('1e-8'))
mp.dps = 40
check('perimeter a=1 edge-on 18sqrt3', perim_edge_closed(1 - mpf('1e-30')), 18*sqrt(3), mpf('1e-10'))
a = mpf('0.05'); C = cos(mpf(60)*pi/180)**2
check('perimeter small-spin series (a=0.05, 60 deg)', perim_general(a, mpf(60)*pi/180)/(6*sqrt(3)*pi),
      1 - (1 + C)*a**2/36 - (mpf(35)/5184 + 35*C/7776 - 23*C**2/15552)*a**4, mpf('1e-9'))

print('== slow-particle capture (E = 1), edge-on ==')
def mb_L(y, a): return -(y**4 - 2*y**3 + a**2)/(a*(y - 1))
def mb_Q(y, a): return y**4*(a**2 - (y**2 - 2*y)**2)/(a**2*(y - 1)**2)
def sigma_closed(a): return 7*pi + pi*sqrt(1 - a**2) + 16*sqrt(1 + a)*ellipe(2*a/(1 + a))
def sigma_direct(a):
    lo, hi = 1 + sqrt(1 - a), 1 + sqrt(1 + a)
    g = lambda t: (lambda y: 2*sqrt(max(mb_Q(y, a), 0))*abs(diff(lambda z: mb_L(z, a), y))*(hi - lo)*sin(t)/2)(lo + (hi - lo)*(1 - cos(t))/2)
    return quad(g, [0, pi/2, pi])
def sigma_reduced(a):
    lo, hi = sqrt(1 - a), sqrt(1 + a)
    g = lambda t: (lambda u: (16*u**2 + 14*u + 2*(1 - a**2)/u)/sqrt(a**2 - (u**2 - 1)**2)*(hi - lo)*sin(t)/2)(lo + (hi - lo)*(1 - cos(t))/2)
    return quad(g, [0, pi/2, pi])
def captured(a, L, Q):          # E = 1 radial potential has no turning point outside the horizon
    rh = 1 + sqrt(1 - a**2)
    R = lambda x: (x*x + a*a - a*L)**2 - (x*x - 2*x + a*a)*(x*x + (L - a)**2 + Q)
    coeffs = sp.Poly(sp.expand(((sp.Symbol('x')**2 + float(a)**2 - float(a)*float(L))**2
                      - (sp.Symbol('x')**2 - 2*sp.Symbol('x') + float(a)**2)*(sp.Symbol('x')**2 + (float(L) - float(a))**2 + float(Q)))), sp.Symbol('x')).all_coeffs()
    import numpy as np
    rts = np.roots([float(cc) for cc in coeffs])
    return not any(abs(z.imag) < 1e-9 and z.real > float(rh) for z in rts)
def sigma_brute(a, nL=1500):
    import numpy as np
    a = float(a); tot = 0.0
    Ls = np.linspace(-6, 6, nL + 1); dL = Ls[1] - Ls[0]
    for L in Ls[:-1] + dL/2:
        if not captured(a, L, 0.0): continue
        lo, hi = 0.0, 40.0
        for _ in range(40):
            mid = (lo + hi)/2
            if captured(a, L, mid): lo = mid
            else: hi = mid
        tot += 2*np.sqrt(lo)*dL
    return mpf(tot)
for a in ['0.3', '0.5', '0.9', '0.999']:
    a = mpf(a)
    check('slow capture closed vs boundary integral a=%s' % a, sigma_closed(a), sigma_direct(a), mpf('1e-15'))
    check('slow capture closed vs reduced integral a=%s' % a, sigma_closed(a), sigma_reduced(a), mpf('1e-15'))
for a in ['0.5', '0.9']:
    check('slow capture closed vs brute force a=%s' % a, sigma_closed(mpf(a)), sigma_brute(mpf(a)), mpf('2e-2'))
check('slow capture a=0: 16pi', sigma_closed(mpf(0)), 16*pi, mpf('1e-30'))
check('slow capture a=1: 7pi+16sqrt2', sigma_closed(mpf(1)), 7*pi + 16*sqrt(2), mpf('1e-30'))
a = mpf('0.05')
check('slow capture series', sigma_closed(a)/(16*pi), 1 - a**2/16 - 31*a**4/2048 - 233*a**6/32768, mpf('1e-10'))

print('== torque moment, edge-on slow capture ==')
def M1_closed(a):
    m = 2*a/(1 + a); K, E = ellipk(m), ellipe(m)
    return -(32*sqrt(1 + a)*((a**2 + 2)*E - 2*(1 - a)*K) + 35*pi*a**2 + 10*pi*(1 - sqrt(1 - a**2)))/(5*a)
def M1_direct(a):
    lo, hi = 1 + sqrt(1 - a), 1 + sqrt(1 + a)
    g = lambda t: (lambda y: mb_L(y, a)*2*sqrt(max(mb_Q(y, a), 0))*abs(diff(lambda z: mb_L(z, a), y))*(hi - lo)*sin(t)/2)(lo + (hi - lo)*(1 - cos(t))/2)
    return quad(g, [0, pi/2, pi])
for a in ['0.1', '0.5', '0.9']:
    a = mpf(a); check('torque M1 closed vs direct a=%s' % a, M1_closed(a), M1_direct(a), mpf('1e-15'))
check('M1 at a->1 (limit; correction ~ sqrt(1-a) ~ 1e-15)', M1_closed(1 - mpf('1e-30')), -(96*sqrt(2) + 45*pi)/5, mpf('1e-13'))
a = mpf('0.05')
check('mean L_z series', M1_closed(a)/sigma_closed(a), -a*(1 + 3*a**2/32 + 37*a**4/1024), mpf('1e-8'))

print('== isotropic slow capture ==')
y = sp.symbols('y', positive=True); A = sp.symbols('a', positive=True)
Ly = -(y**4 - 2*y**3 + A**2)/(A*(y - 1)); Qy = y**4*(A**2 - (y**2 - 2*y)**2)/(A**2*(y - 1)**2)
scheck('Carter C = Q + L^2', sp.factor(Qy + Ly**2) - (3*y**4 - 4*y**3 + A**2)/(y - 1)**2)
def iso_direct(a, n=0):
    lo, hi = 1 + sqrt(1 - a), 1 + sqrt(1 + a)
    Lc2 = lambda y: (3*y**4 - 4*y**3 + a**2)/(y - 1)**2
    mu = lambda y: -(y**4 - 2*y**3 + a**2)/(a*sqrt(3*y**4 - 4*y**3 + a**2))
    g = lambda t: (lambda y: Lc2(y)*abs(diff(mu, y))*(hi - lo)*sin(t)/2)(lo + (hi - lo)*(1 - cos(t))/2)
    return pi/2*quad(g, [0, pi/2, pi])
def iso_formula_int(a):
    lo, hi = 1 + sqrt(1 - a), 1 + sqrt(1 + a)
    g = lambda t: (lambda y: 2*y**3*(3*y**4 - 8*y**3 + 6*y**2 - a**2)/(a*(y - 1)**2*sqrt(3*y**4 - 4*y**3 + a**2))*(hi - lo)*sin(t)/2)(lo + (hi - lo)*(1 - cos(t))/2)
    return pi/2*abs(quad(g, [0, pi/2, pi]))
for a in ['0.5', '0.9']:
    a = mpf(a); check('isotropic integral formula vs direct a=%s' % a, iso_formula_int(a), iso_direct(a), mpf('1e-15'))
iso1 = pi*(mpf(47)/9 + 31*sqrt(2)/6 + 7*sqrt(6)/9 - sqrt(3)/54*log((1 + sqrt(2) + sqrt(6))/(1 + sqrt(2))))
# a = 1: curve on 1 <= y <= 1+sqrt2 (mu from -1 to sqrt(2/3)) plus straight piece C = 4/mu^2 for mu in [sqrt(2/3), 1]
Lc2 = lambda y: 3*y**2 + 2*y + 1
mu1 = lambda y: -(y**3 - y**2 - y - 1)/sqrt(3*y**2 + 2*y + 1)
part1 = quad(lambda y: Lc2(y)*abs(diff(mu1, y)), [1, 1 + sqrt(2)])
part2 = quad(lambda m: 4/m**2, [sqrt(mpf(2)/3), 1])
check('isotropic a=1 closed form vs direct', iso1, pi/2*(part1 + part2), mpf('1e-15'))
check('isotropic a=1 value', iso1, mpf('45.27564306526540429854285'), mpf('1e-15'))
check('isotropic a=1 vs a->1 integral formula', iso_formula_int(1 - mpf('1e-12')), iso1, mpf('1e-4'))
num = quad(lambda y: mu1(y)*Lc2(y)**mpf(1.5)*abs(diff(mu1, y)), [1, 1 + sqrt(2)]) + quad(lambda m: m*8/m**3, [sqrt(mpf(2)/3), 1])
numclosed = -mpf(4876)/189 - mpf(790)/63*sqrt(2) + mpf(628)/189*sqrt(6) - 2*sqrt(3)/81*log((1 + sqrt(2) + sqrt(6))/(1 + sqrt(2)))
check('isotropic a=1 L_z moment closed form', num, numclosed, mpf('1e-15'))
lz = mpf(2)/3*numclosed/(2*iso1/pi)
check('isotropic a=1 mean l_z', lz, mpf('-0.81932684180255246848'), mpf('1e-18'))

print('\nSUMMARY: %d failures' % len(FAIL), FAIL)
