from mpmath import mp, mpf, sqrt, quad, pi, findroot, diff, cos, sin
mp.dps=30
def area(v):
    v=mpf(v); E=1/sqrt(1-v*v); p2=E*E-1
    w=lambda r: sqrt(r*(p2*r+1))
    L=lambda r: E*(r+1)-(r-1)*w(r)
    Q=lambda r: r*r*(-p2*r*r+(2*p2-1)*r+1+2*E*w(r))
    rmax=findroot(Q,4.0)
    f=lambda t:(lambda r: sqrt(max(Q(r),0))*abs(diff(L,r))*(rmax-1)*sin(t)/2)(1+(rmax-1)*(1-cos(t))/2)
    return 2*quad(f,[0,pi/2,pi])/p2, rmax
if __name__=="__main__":
    for v in ['0.001','0.3','0.6','0.9','0.999','0.99999']:
        A,rm=area(v); print(v, A, 'times v^2 (slow norm):', A*mpf(v)**2, ' times p^2/E^2:', A*(1-mpf(v)**2)**0*mpf(v)**2)
    print('slow limit 7pi+16sqrt2 =',7*pi+16*sqrt(2),'  photon 16pi+15sqrt3 =',16*pi+15*sqrt(3))
