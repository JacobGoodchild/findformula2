# Inner shadow radius for a polar observer of extremal Kerr: photon (L_z=0, b^2 = eta + a^2) whose first equatorial crossing is at the horizon
from mpmath import mp, mpf, sqrt, quad, ellipk, pi, findroot, inf, identify
mp.dps=30
a=mpf(1)
def Gth(b):  # Int_0^{pi/2} dtheta/sqrt(eta + a^2 cos^2), eta=b^2-a^2
    return ellipk(a*a/(b*b))/b
def Gr(b):
    R=lambda r:(r*r+a*a)**2-(r-1)**2*b*b
    return quad(lambda r: 1/sqrt(R(r)),[1,2,10,inf])
bin=findroot(lambda b: Gth(b)-Gr(b),3.0)
print('b_in =',bin, ' b_in^2=',bin**2, ' area pi b^2 =',pi*bin**2)
print(identify(bin),identify(bin**2,['sqrt(2)','sqrt(3)','sqrt(5)']))
# also critical (shadow) radius for comparison: pi(12+8sqrt2)
print('shadow radius', sqrt(12+8*sqrt(2)))
