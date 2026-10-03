from mpmath import mp, mpf, sqrt, cos, sin, pi, findroot, diff, nstr
mp.dps=40
def moment(fun,r_c,a,th,n=48):
    F0=fun(r_c); F1=diff(fun,r_c)
    E=sqrt(2*F0**2/(2*F0-r_c*F1)); L0=sqrt(F1*r_c**3/(2*F0-r_c*F1)); p=sqrt(E*E-1); b0=L0/p
    st,ct=sin(th),cos(th); m=0; ar=0
    for k in range(n):
        ph=2*pi*(k+mpf(1)/2)/n
        def eqs(rr,b):
            al=b*cos(ph); be=b*sin(ph); Lz=p*al*st; Q=p**2*(be**2+(al**2-a**2)*ct**2)
            R=lambda x:(E*(x*x+a*a)-a*Lz)**2-(x*x*fun(x)+a*a)*(x*x+(Lz-a*E)**2+Q)
            return [R(rr),diff(R,rr)]
        rr,bb=findroot(eqs,(r_c,b0)); m+=p*st*cos(ph)*bb**3/3; ar+=bb**2/2
    return m/ar/E
fun=lambda x: 1-2*x**2/(x**2+mpf('0.16'))**1.5
for r_c,thdeg in [(mpf('3.4'),60),(mpf('3.1'),35)]:
    th=thdeg*pi/180; a=mpf('1e-6')
    print('Bardeen g=0.4',r_c,thdeg,nstr(moment(fun,r_c,a,th)/a,12), nstr(-sin(th)**2*(1-fun(r_c))/fun(r_c),12))
