# Exact equatorial silhouette solid angle for extremal Kerr at r* = M/z (Baines-Visser edge), as 1D integrals
from mpmath import mp, mpf, sqrt, pi, quad, asin, atan, nstr
import escape
mp.dps=30
def Omega_eq(z):
    z=mpf(z); k=z/(1+z)
    def cosI(w):
        return (1-z-w*z)*sqrt((w*z+1-z)**2+4*z*(1-z))/(1+(w*w-1)*z*z)
    from mpmath import cos as mc, sin as ms
    # Case I: 2 Int_0^3 cosTheta(w) dw/sqrt((3-w)(1+w)); substitute w = 1 + 2 sin(phi) (phi in [-pi/6, pi/2])
    I1=2*quad(lambda p: cosI(1+2*ms(p)), [-pi/6, pi/2])
    # Case II: 2 Int_{1/2}^{1} sqrt(s^2-k^2)/(s sqrt(1-s^2)) ds  (s=|sin Phi|)
    I2=2*quad(lambda p: sqrt(ms(p)**2-k*k)/ms(p), [pi/6, pi/2])
    return 2*pi-I1-I2, I2
if __name__=='__main__':
 for z in ['0.01','0.1','0.3','0.5','0.9']:
  O,I2=Omega_eq(z)
  print(z, O, escape.Omega(mpf(z),pi/2))
