# Closed form (elliptic K,E,Pi) for extremal Kerr shadow area at inclination th
from mpmath import *
from ext_formula import roots, A_new, A_direct
mp.dps=40
def A_closed(th):
    s=sin(th); c2=cos(th)**2
    d,c,b,a=roots(th)            # rA<rB<rC<rD
    m=(a-b)*(c-d)/((a-c)*(b-d)); n=(a-b)/(a-c); g=2/sqrt((a-c)*(b-d))
    def pieces(phi):  # integrals over u in [0,u(phi)]
        u=ellipf(phi,m); sn=sin(phi); cn=cos(phi); dn=sqrt(1-m*sn**2)
        P=ellippi(n,phi,m); E=ellipe(phi,m)
        V=(n*E+(m-n)*u+(2*n*m+2*n-n*n-3*m)*P-n*n*sn*cn*dn/(1-n*sn**2))/(2*(n-1)*(m-n))
        I0=g*u; I2=g*(c*c*u+2*c*(b-c)*P+(b-c)**2*V)
        return 4*(3*I2-c2*I0)
    full=pieces(pi/2)
    if b>=1: return full
    phi0=asin(sqrt((a-c)*(1-b)/((a-b)*(1-c))))
    return full-pieces(phi0)+(2+s*s)*sqrt(3*s*s+s*s*c2-4*c2)/s**2
if __name__=='__main__':
  for deg in ['5','30','45','50','60','80','89.9']:
    th=mpf(deg)*pi/180; print(deg, nstr(A_closed(th),32), nstr(A_new(th),32))
