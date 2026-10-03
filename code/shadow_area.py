# Area of Kerr shadow (critical curve) seen by a distant observer at inclination theta_o.
from mpmath import mp, mpf, sqrt, quad, cos, sin, cot, findroot, pi, acos
mp.dps = 30
def xi(r,a,M=1): return (r**2*(3*M-r)-a**2*(r+M))/(a*(r-M))
def eta(r,a,M=1): return r**3*(4*a**2*M-r*(r-3*M)**2)/(a**2*(r-M)**2)
def roots(a):
    # equatorial photon orbit radii
    rp = 2*(1+cos(2*acos(a)/3)); rm = 2*(1+cos(2*acos(-a)/3))
    return rp, rm
def area(a, th=pi/2):
    rp, rm = roots(a)
    # alpha = -xi/sin th, beta = ±sqrt(eta + a^2 cos^2 - xi^2 cot^2)
    def beta(r): 
        v = eta(r,a)+a**2*cos(th)**2-xi(r,a)**2*cot(th)**2
        return sqrt(v) if v>0 else mpf(0)
    def dalpha(r): return -mp.diff(lambda s: xi(s,a), r)/sin(th)
    # find r-range where beta^2>=0
    return 2*quad(lambda r: beta(r)*abs(dalpha(r)), [rp, (rp+rm)/2, rm])
if __name__=="__main__":
    for a in [mpf('1e-6'), mpf('0.5'), mpf('0.9'), mpf('0.999999')]:
        print(a, area(a), area(a)/(27*pi))
