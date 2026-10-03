# Silhouette solid angle, extremal Kerr-Newman (M=1, a^2+Q^2=1), equatorial Carter-frame observer at r=1/z
from mpmath import mp, mpf, sqrt, pi, quad, asin, atan, sin, cos
mp.dps=30
def cosI(w,z):
    return (1-z-w*z)*sqrt((w*z+1-z)**2+4*z*(1-z))/(1+(w*w-1)*z*z)
def Omega_kn(z,a):
    z=mpf(z); a=mpf(a); k=z/(a*(1+z))
    if a>mpf(1)/2:
        # Case I: sigma in [-1/(2a),1], w = 1+2a sigma ; param sigma=sin(p)
        p0=asin(-1/(2*a))
        I1=2*quad(lambda p: cosI(1+2*a*sin(p),z),[p0,pi/2])
        I2=2*asin(sqrt(4*a*a-1)/(2*a*sqrt(1-k*k)))-2*k*atan(k*sqrt(4*a*a-1)/sqrt(1-4*a*a*k*k))
    else:
        I1=2*quad(lambda p: cosI(1+2*a*sin(p),z),[-pi/2,pi/2]); I2=0
    return 2*pi-I1-I2
if __name__=="__main__":
    import kn_ext_check
    for a in ['0.3','0.8','1']:
        for z in ['1e-4']:
            print('far-field', a, Omega_kn(z,a)/mpf(z)**2, kn_ext_check.A_formula(float(a)))
    import escape
    print('a=1 vs Kerr:', Omega_kn('0.3',1), escape.Omega(mpf('0.3'),pi/2))
