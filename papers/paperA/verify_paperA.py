# Verification of every displayed formula in paperA.tex.
# Formulas are re-typed here from the paper (not imported from the derivation scripts) and
# compared with (i) symbolic reductions and (ii) direct high-precision numerics of the exact
# shadow contour / exact capture boundary.  Run from the repository root:
#     python3 papers/paperA/verify_paperA.py
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
import sympy as sp
from mpmath import mp, mpf, sqrt, pi, cos, sin, findroot, diff, matrix, lu_solve, nstr, quad
mp.dps = 40
FAIL = []
def check(name, x, y, tol):
    d = abs(complex(x) - complex(y))
    ok = d < tol
    print(('OK  ' if ok else 'FAIL') + ' %s: %s vs %s (diff %.1e)' % (name, nstr(x, 16), nstr(y, 16), d))
    if not ok: FAIL.append(name)
def scheck(name, expr):
    z = sp.simplify(expr)
    ok = (z == 0)
    print(('OK  ' if ok else 'FAIL') + ' %s (symbolic difference = %s)' % (name, z))
    if not ok: FAIL.append(name)

# ---------- formulas as written in the paper ----------
def K2_shadow(f, kappa, R, C):            # eq (A2) bracket
    return (1 - C) * (1/(2*f*R**2) - 4/(R**2*kappa)) + C*(2*f - 1)/(f*R**2)
def K4_shadow(f, kappa, f3, f4, R, C):    # eqs (A4)-(N2)
    k = kappa
    N0 = -k**5 + 16*f*k**4 - 160*f**2*k**3 + 32*f**2*f3*k**2 - 384*f**3*k**2 + 256*f**3*f3*k - 32*f**3*f4*k - 96*f**3*f3**2
    N1 = 2*(-3*k**5 + 4*f*k**5 + 32*f*k**4 - 32*f**2*k**4 - 160*f**2*k**3 + 32*f**2*f3*k**2 + 192*f**3*k**3
            + 384*f**3*k**2 - 64*f**3*f3*k**2 - 256*f**3*f3*k + 32*f**3*f4*k + 96*f**3*f3**2)
    N2 = (15*k**5 - 24*f*k**5 - 144*f*k**4 + 8*f**2*k**5 + 192*f**2*k**4 + 480*f**2*k**3 - 96*f**2*f3*k**2
          - 64*f**3*k**4 - 384*f**3*k**3 - 384*f**3*k**2 + 128*f**3*f3*k**2 + 256*f**3*f3*k - 32*f**3*f4*k - 96*f**3*f3**2)
    return (N0 + N1*C + N2*C**2)/(8*f**2*k**5)/R**4
def K2_capture(f, fp, fpp, r, C):         # eq (sig2)
    J = 3*f*fp + r*f*fpp - 2*r*fp**2
    return (1 - C)*(1/(2*f*r**2) + fp**2/(r*f*J)) + C*(2*f + r*fp - 2)/(r**3*fp)
def K4_pole(f, fp, fpp, r):               # eq (sigpole), a^4 coefficient
    J = 3*f*fp + r*f*fpp - 2*r*fp**2
    return 2*(1 - f)**2*(3*fp + r*fpp)/(r**5*fp*J)
def torque1(f): return -(1 - f)/f         # eq (torque): <L>/E = a sin^2 * this

# ---------- symbolic checks ----------
print('== symbolic checks ==')
r, q, C, R = sp.symbols('r q C R', positive=True)
fS = sp.Function('f')
# capture bracket reduces to shadow bracket at the photon sphere (r f' = 2f)
f0, f1, f2 = sp.symbols('f0 f1 f2')
kap = 2*f0 - r**2*f2
scheck('capture -> shadow at photon sphere', (K2_capture(f0, 2*f0/r, f2, r, C) - K2_shadow(f0, kap, r, C)))
# pole a^4 light limit equals A4 at C=1 is checked numerically below (needs f3,f4).
# Kerr specialisations
fK = 1 - 2/r
def ders(fexpr, x): return [sp.diff(fexpr, r, k).subs(r, x) for k in range(5)]
d = ders(fK, 3)
scheck('Kerr shadow a^2', K2_shadow(d[0], 2*d[0] - 9*d[2], 3, C) + (1 + C)/18)
scheck('Kerr shadow a^4', K4_shadow(d[0], 2*d[0] - 9*d[2], 27*d[3], 81*d[4], 3, C)
       + (sp.Rational(1, 72) + 5*C/972 - 5*C**2/1944))
f_, fp_, fpp_ = fK, sp.diff(fK, r), sp.diff(fK, r, 2)
scheck('Kerr capture a^2', K2_capture(f_, fp_, fpp_, r, C) + (r + 3*C*(4 - r))/(2*r**2*(6 - r)))
scheck('Kerr pole a^2', K2_capture(f_, fp_, fpp_, r, 1) + 1/r**2)
scheck('Kerr pole a^4', K4_pole(f_, fp_, fpp_, r) + 4/(r**4*(6 - r)))
scheck('Kerr torque', torque1(f_) + 2/(r - 2))
scheck('Kerr sigma0', sp.pi*r**3*fp_/(2*f_**2 - 2*f_ + r*fp_) - sp.pi*r**3/(4 - r))
v = sp.symbols('v', positive=True)
rv = (4*v**2 - 1 + sp.sqrt(1 + 8*v**2))/(2*v**2)
scheck('Kerr r(v): E^2 = 1/(1-v^2)', sp.simplify((2*fK**2/(2*fK - r*sp.diff(fK, r))).subs(r, rv) - 1/(1 - v**2)))
# Kerr-Newman
fKN = 1 - 2/r + q**2/r**2
f_, fp_, fpp_ = fKN, sp.diff(fKN, r), sp.diff(fKN, r, 2)
KNsig = -(r**2*(r - q**2) - C*(3*r**3 - 12*r**2 - q**2*r**2 + 18*q**2*r - 8*q**4))/(2*r*(r - q**2)*(6*r**2 - r**3 - 9*q**2*r + 4*q**4))
scheck('KN capture a^2', K2_capture(f_, fp_, fpp_, r, C) - KNsig)
scheck('KN sigma0', sp.pi*r**3*fp_/(2*f_**2 - 2*f_ + r*fp_) - sp.pi*r**4*(r - q**2)/(4*r**2 - r**3 - 4*q**2*r + q**4))
scheck('KN E^2', 2*f_**2/(2*f_ - r*fp_) - (r**2 - 2*r + q**2)**2/(r**2*(r**2 - 3*r + 2*q**2)))
scheck('KN torque', torque1(f_) + (2*r - q**2)/(r**2 - 2*r + q**2))
# KN shadow at the photon ring: substitute q^2 = (3R - R^2)/2
q2 = (3*R - R**2)/2
fR = (1 - 2/r + q2/r**2)
dd = [sp.diff(fR, r, k).subs(r, R) for k in range(5)]
scheck('KN shadow a^2', K2_shadow(dd[0], 2*dd[0] - R**2*dd[2], R, C) + (R + 3*C*(R - 2))/(R**2*(R - 1)*(2*R - 3)))
K4KN = (C**2*(R**6 - 38*R**5 + 286*R**4 - 768*R**3 + 774*R**2 - 108*R - 162)
        - C*(18*R**6 - 108*R**5 + 288*R**4 - 492*R**3 + 504*R**2 - 216*R)
        - 15*R**6 + 58*R**5 - 78*R**4 + 36*R**3)/(2*R**4*(R - 1)**2*(2*R - 3)**5)
scheck('KN shadow a^4', K4_shadow(dd[0], 2*dd[0] - R**2*dd[2], R**3*dd[3], R**4*dd[4], R, C) - K4KN)
scheck('KN shadow A0', sp.pi*R**2/dd[0] - sp.pi*R**4/(R**2 - 2*R + q2))
# KN double series in q (to q^2 a^2)
Rq = (3 + sp.sqrt(9 - 8*q**2))/2
A0q = sp.pi*Rq**4/(Rq**2 - 2*Rq + q**2)
Aq = A0q*(1 - sp.Symbol('a')**2*(Rq + 3*C*(Rq - 2))/(Rq**2*(Rq - 1)*(2*Rq - 3)))
ser = sp.series(Aq/(27*sp.pi), q, 0, 5).removeO()
aa = sp.Symbol('a')
target = 1 - q**2/3 - (1 + C)*aa**2/18 - q**4/27 - (3 + C)*aa**2*q**2/81
diffser = sp.expand(ser - target)
# keep only terms up to a^2 and q^4 (series truncated); a^2 q^4 terms are beyond the stated order
diffser = sum(t for t in sp.Add.make_args(diffser) if sp.degree(t, q) <= 2 or sp.degree(t, aa) == 0)
scheck('KN double series', diffser)
# Kerr isotropic a^3 torque (from eqs kerrsig, kerrL3)
cc = sp.symbols('c')
Cc = cc**2
sig = 1 - aa**2*(r + 3*Cc*(4 - r))/(2*r**2*(6 - r))
ell = -2*aa*(1 - Cc)/(r - 2) + aa**3*(1 - Cc)*(sp.Rational(3, 2)*(r - 5) - Cc*(7*r**2 - 63*r + 144)/(2*r))/((r - 2)*(6 - r)**3)
num = sp.integrate(sp.expand(sig*ell), (cc, 0, 1)); den = sp.integrate(sig, (cc, 0, 1))
Liso = sp.series(num/den, aa, 0, 4).removeO()
scheck('Kerr isotropic torque', Liso - (-4*aa/(3*(r - 2)) + 4*aa**3*(3*r**3 - 19*r**2 + 48*r - 144)/(15*r**2*(r - 2)*(6 - r)**3)))
scheck('isotropic slow limit', sp.expand(Liso.subs(r, 4)) + 2*aa/3 + aa**3/15)
scheck('isotropic light limit', sp.expand(Liso.subs(r, 3)) + 4*aa/3 + 8*aa**3/81)
# Kerr full a^4 capture: consistency
a8 = -r**2*(2*r**3 + 42*r**2 - 549*r + 1458); b8 = 6*r*(2*r**4 + 2*r**3 - 239*r**2 + 1302*r - 2016)
g8 = -10*r**5 - 2*r**4 + 1653*r**3 - 13266*r**2 + 39744*r - 41472
scheck('alpha8+beta8+gamma8 = -32(6-r)^4', a8 + b8 + g8 + 32*(6 - r)**4)
scheck('Kerr a^4 capture at r=3 = shadow a^4', ((a8 + b8*C + g8*C**2)/(8*r**4*(6 - r)**5)).subs(r, 3)
       + (sp.Rational(1, 72) + 5*C/972 - 5*C**2/1944))
# Kerr a^6, a^8 from the stored universal coefficients
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code')
for n, fn, kmax in [(6, 'univ_a6_fast.out', 10), (8, 'univ_a8_fast.out', 14)]:
    fs = sp.symbols('f0:16'); Csym = sp.Symbol('C')
    e = sp.sympify(open(os.path.join(base, fn)).read().split('a^%d relative (R=1):' % n)[1].split('\n')[0],
                   locals={'C': Csym, **{'f%d' % i: fs[i] for i in range(16)}})
    x = sp.Symbol('x')
    val = sp.expand(sp.nsimplify(e.subs({fs[k]: sp.diff(1 - 2/x, x, k).subs(x, 3)*3**k for k in range(kmax)})/3**n))
    want = {6: -(sp.Rational(2059, 314928) + 167*Csym/104976 - 185*Csym**2/104976 - 5*Csym**3/3888),
            8: -(sp.Rational(88175, 22674816) + 3869*Csym/5668704 - 4865*Csym**2/3779136 - 83*Csym**3/69984 - 667*Csym**4/7558272)}[n]
    scheck('Kerr shadow a^%d from stored universal coefficient' % n, val - want)

# ---------- numerical checks against exact geometry ----------
print('== numerical checks (exact contour / exact capture boundary) ==')
def shadow_area(fun, a, th):
    D = lambda x: x*x*fun(x) + a*a
    Dp = lambda x: diff(D, x)
    xi = lambda x: ((x*x + a*a)*Dp(x) - 4*x*D(x))/(a*Dp(x))
    eta = lambda x: x*x*(16*a*a*D(x) - (4*D(x) - x*Dp(x))**2)/(a*a*Dp(x)**2)
    s, c = sin(th), cos(th)
    b2 = lambda x: eta(x) + a*a*c*c - xi(x)**2*c*c/(s*s)
    Rp = findroot(lambda x: 2*fun(x) - x*diff(fun, x), 3)
    r1 = findroot(b2, Rp - 1.8*a); r2 = findroot(b2, Rp + 1.8*a)
    if r1 > r2: r1, r2 = r2, r1
    g = lambda t: (lambda x: sqrt(max(b2(x), 0))*abs(diff(xi, x))/s*(r2 - r1)*sin(t)/2)(r1 + (r2 - r1)*(1 - cos(t))/2)
    return 2*quad(g, [0, pi/2, pi]), Rp
def capture_area(fun, rc, a, th, n=48, moment=False):
    F0 = fun(rc); F1 = diff(fun, rc)
    E = sqrt(2*F0**2/(2*F0 - rc*F1)); L0 = sqrt(F1*rc**3/(2*F0 - rc*F1)); p = sqrt(E*E - 1); b0 = L0/p
    st, ct = sin(th), cos(th); tot = 0; mom = 0
    for k in range(n):
        ph = 2*pi*(k + mpf(1)/2)/n
        def eqs(rr, b):
            al = b*cos(ph); be = b*sin(ph); Lz = p*al*st; Q = p**2*(be**2 + (al**2 - a**2)*ct**2)
            Rf = lambda x: (E*(x*x + a*a) - a*Lz)**2 - (x*x*fun(x) + a*a)*(x*x + (Lz - a*E)**2 + Q)
            return [Rf(rr), diff(Rf, rr)]
        rr, bb = findroot(eqs, (rc, b0))
        tot += bb**2; mom += p*st*cos(ph)*bb**3/3
    if moment: return (mom/(tot/2))/E
    return pi*tot/n, pi*b0**2
def fitcoef(vals, hs, nterms):
    c = lu_solve(matrix([[h**(2*j) for j in range(nterms)] for h in hs]), matrix(vals))
    return c
bard = lambda x: 1 - 2*x**2/(x**2 + mpf('0.16'))**mpf(1.5)
hayw = lambda x: 1 - 2*x**2/(x**3 + mpf('0.5'))
kerr = lambda x: 1 - 2/x
def kn(qq): return lambda x: 1 - 2/x + qq**2/x**2
def shadow_coeffs(fun, thdeg):
    th = thdeg*pi/180; hs = [mpf(k)/100 for k in range(1, 6)]; ys = []
    for h in hs:
        A, Rp = shadow_area(fun, h, th); ys.append((A/(pi*Rp**2/fun(Rp)) - 1)/h**2)
    c = fitcoef(ys, hs, 5)
    f = fun(Rp); kap = 2*f - Rp**2*diff(fun, Rp, 2)
    C = cos(th)**2
    return c[0], K2_shadow(f, kap, Rp, C), c[1], K4_shadow(f, kap, Rp**3*diff(fun, Rp, 3), Rp**4*diff(fun, Rp, 4), Rp, C)
for name, fun, thdeg in [('Bardeen 90', bard, 90), ('Bardeen 50', bard, 50), ('Hayward 90', hayw, 90),
                         ('Kerr 50', kerr, 50), ('KN q=0.6 60', kn(mpf('0.6')), 60), ('KN q=0.8 90', kn(mpf('0.8')), 90)]:
    n2, f2_, n4, f4_ = shadow_coeffs(fun, thdeg)
    check('shadow a^2 ' + name, n2, f2_, 1e-9)
    check('shadow a^4 ' + name, n4, f4_, 1e-8)
def capture_coeffs(fun, rc, thdeg):
    th = thdeg*pi/180; hs = [mpf('0.004'), mpf('0.008'), mpf('0.012')]; ys = []
    for h in hs:
        A, A0 = capture_area(fun, rc, h, th); ys.append((A/A0 - 1)/h**2)
    return fitcoef(ys, hs, 3)[0], K2_capture(fun(rc), diff(fun, rc), diff(fun, rc, 2), rc, cos(th)**2)
for name, fun, rc, thdeg in [('Bardeen r=3.4 60', bard, mpf('3.4'), 60), ('Hayward r=3.1 90', hayw, mpf('3.1'), 90),
                             ('KN q=0.5 r=3.5 60', kn(mpf('0.5')), mpf('3.5'), 60), ('KN q=0.8 r=3.0 45', kn(mpf('0.8')), mpf('3.0'), 45)]:
    nn, ff = capture_coeffs(fun, rc, thdeg); check('capture a^2 ' + name, nn, ff, 1e-11)
for rc in [mpf('3.4'), mpf('3.1')]:
    hs = [mpf(k)/100 for k in range(1, 6)]; ys = []
    for h in hs:
        A, A0 = capture_area(bard, rc, h, mpf(0), n=8); ys.append((A/A0 - 1)/h**2)
    c = fitcoef(ys, hs, 5)
    check('pole a^2 Bardeen r=%s' % rc, c[0], K2_capture(bard(rc), diff(bard, rc), diff(bard, rc, 2), rc, 1), 1e-11)
    check('pole a^4 Bardeen r=%s' % rc, c[1], K4_pole(bard(rc), diff(bard, rc), diff(bard, rc, 2), rc), 1e-9)
for name, fun, rc, thdeg in [('Bardeen r=3.4 60', bard, mpf('3.4'), 60), ('KN q=0.5 r=3.5 60', kn(mpf('0.5')), mpf('3.5'), 60)]:
    a = mpf('1e-8'); th = thdeg*pi/180
    check('torque ' + name, capture_area(fun, rc, a, th, moment=True)/a, sin(th)**2*torque1(fun(rc)), 1e-6)
# light limit of the pole a^4 formula equals the shadow a^4 at C = 1 (Bardeen)
Rp = findroot(lambda x: 2*bard(x) - x*diff(bard, x), 3)
f = bard(Rp); kap = 2*f - Rp**2*diff(bard, Rp, 2)
check('pole a^4 light limit = shadow a^4 at C=1', K4_pole(f, diff(bard, Rp), diff(bard, Rp, 2), Rp),
      K4_shadow(f, kap, Rp**3*diff(bard, Rp, 3), Rp**4*diff(bard, Rp, 4), Rp, 1), 1e-25)
# displacement of the shadow centre (Bardeen, 60 deg)
a = mpf('1e-5'); th = 60*pi/180
D = lambda x: x*x*bard(x) + a*a; Dp = lambda x: diff(D, x)
xi = lambda x: ((x*x + a*a)*Dp(x) - 4*x*D(x))/(a*Dp(x)); eta = lambda x: x*x*(16*a*a*D(x) - (4*D(x) - x*Dp(x))**2)/(a*a*Dp(x)**2)
b2 = lambda x: eta(x) + a*a*cos(th)**2 - xi(x)**2*cos(th)**2/sin(th)**2
r1 = findroot(b2, Rp - 1.8*a); r2 = findroot(b2, Rp + 1.8*a)
check('centre displacement / a', abs(xi(r1) + xi(r2))/(2*sin(th))/a, sin(th)*(1 - f)/f, 1e-6)

# Kerr: numerically identified full a^4 capture term (eq. kerr4) and a^3 torque (eq. kerrL3), new points
def kerr_K4(rr, C):
    a8 = -rr**2*(2*rr**3 + 42*rr**2 - 549*rr + 1458); b8 = 6*rr*(2*rr**4 + 2*rr**3 - 239*rr**2 + 1302*rr - 2016)
    g8 = -10*rr**5 - 2*rr**4 + 1653*rr**3 - 13266*rr**2 + 39744*rr - 41472
    return (a8 + b8*C + g8*C**2)/(8*rr**4*(6 - rr)**5)
for rc, thdeg in [(mpf('3.55'), 40), (mpf('3.15'), 70)]:
    th = thdeg*pi/180; hs = [mpf(k)/1000*3 for k in range(1, 6)]; ys = []
    for h in hs:
        A, A0 = capture_area(kerr, rc, h, th, n=64); ys.append((A/A0 - 1)/h**2)
    c = fitcoef(ys, hs, 5)
    check('Kerr capture a^4 (identified) r=%s %d deg' % (rc, thdeg), c[1], kerr_K4(rc, cos(th)**2), 1e-8)
for rc, thdeg in [(mpf('3.55'), 40), (mpf('3.15'), 70)]:
    th = thdeg*pi/180; C = cos(th)**2; hs = [mpf(k)/1000*3 for k in range(1, 6)]; ys = []
    for h in hs: ys.append(capture_area(kerr, rc, h, th, n=64, moment=True)/h)
    c = fitcoef(ys, hs, 5)
    pred = (1 - C)*(mpf(3)/2*(rc - 5) - C*(7*rc**2 - 63*rc + 144)/(2*rc))/((rc - 2)*(6 - rc)**3)
    check('Kerr torque a^3 (identified) r=%s %d deg' % (rc, thdeg), c[1], pred, 1e-8)
print('\nSUMMARY: %d failures' % len(FAIL), FAIL)
