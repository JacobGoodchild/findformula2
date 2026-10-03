from mpmath import mp, mpf, sqrt, pi, findroot, diff, sin, cos, nstr, matrix, lu_solve
mp.dps=22
def perim(fun,a,th,n=6000):
    D=lambda r: r*r*fun(r)+a*a; Dp=lambda r: diff(D,r)
    xi=lambda r:((r*r+a*a)*Dp(r)-4*r*D(r))/(a*Dp(r)); eta=lambda r: r*r*(16*a*a*D(r)-(4*D(r)-r*Dp(r))**2)/(a*a*Dp(r)**2)
    s,c=sin(th),cos(th); b2=lambda r: eta(r)+a*a*c*c-xi(r)**2*c*c/(s*s)
    rph=findroot(lambda r: 2*fun(r)-r*diff(fun,r),3)
    r1=findroot(b2,rph-1.8*a); r2=findroot(b2,rph+1.8*a)
    pts=[]
    for k in range(n+1):
        t=pi*k/n; r=r1+(r2-r1)*(1-cos(t))/2
        pts.append((-xi(r)/s, sqrt(max(b2(r),0))))
    L=sum(sqrt((pts[i+1][0]-pts[i][0])**2+(pts[i+1][1]-pts[i][1])**2) for i in range(n))
    return 2*L, rph
fun=lambda x: 1-2*x**2/(x**2+mpf('0.16'))**1.5
th=60*pi/180; C=cos(th)**2
hs=[mpf('0.02'),mpf('0.03'),mpf('0.04')]; ys=[]
for h in hs:
    P,rph=perim(fun,h,th); P0=2*pi*rph/sqrt(fun(rph)); ys.append((P/P0-1)/h**2)
Kp=lu_solve(matrix([[1,h**2,h**4] for h in hs]),matrix(ys))[0]
f=fun(rph); k=2*f-rph**2*diff(fun,rph,2)
K=(1-C)*(1/(2*f*rph**2)-4/(rph**2*k))+C*(2*f-1)/(f*rph**2)
print('perimeter a^2 coeff numeric',nstr(Kp,10),' K/2 =',nstr(K/2,10))
