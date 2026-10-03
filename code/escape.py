# Solid angle of the extremal-Kerr silhouette for an observer at r*=M/z, declination theta (Baines-Visser edge)
from mpmath import mp, mpf, sqrt, sin, cos, pi, quad, asin
mp.dps=30
def cosTheta(Phi,z,th):
    s=sin(th); c2=cos(th)**2; sp=sin(Phi)
    w=s*sp+sqrt((1+s*sp)**2+c2)
    if w>=0:
        return (1-z-w*z)*sqrt((w*z+1-z)**2+4*z*(1-z))/(1+(w*w-1)*z*z)
    k=z*(1+c2)/((1+z)*s)
    return sqrt(1-(k/sp)**2)
def Omega(z,th):
    s=sin(th); c2=cos(th)**2
    pts=[0,pi/2,pi,3*pi/2,2*pi]
    if s>sqrt(3)-1:
        p=asin(-(1+c2)/(2*s)); pts=sorted(pts+[pi-p, 2*pi+p])
    return 2*pi-quad(lambda P: cosTheta(P,z,th),pts)
if __name__=="__main__":
    import ext_closed
    for deg in [30,60,90]:
        th=mpf(deg)*pi/180
        A=ext_closed.A_closed(th)
        for z in [mpf('1e-3'),mpf('1e-4')]:
            print(deg, z, Omega(z,th)/z**2, A)
