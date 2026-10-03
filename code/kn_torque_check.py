from kn_anyv_check import *
def moment(r_c,q,a,th,n=48):
    E=E_rc(r_c,q); p=sqrt(E*E-1); st,ct=sin(th),cos(th); m=0; ar=0; b0=b0_rc(r_c,q)
    for k in range(n):
        ph=2*pi*(k+mpf(1)/2)/n
        def eqs(r,b):
            al=b*cos(ph); be=b*sin(ph); Lz=p*al*st; Q=p**2*(be**2+(al**2-a**2)*ct**2)
            D=r*r-2*r+a*a+q*q; R=(E*(r*r+a*a)-a*Lz)**2-D*(r*r+(Lz-a*E)**2+Q)
            dR=4*r*E*(E*(r*r+a*a)-a*Lz)-(2*r-2)*(r*r+(Lz-a*E)**2+Q)-D*2*r
            return [R,dR]
        rr,bb=findroot(eqs,(r_c,b0)); m+=p*st*cos(ph)*bb**3/3; ar+=bb**2/2
    return m/ar/E
for r_c,q,thdeg in [(mpf('3.5'),mpf('0.5'),60),(mpf('3.0'),mpf('0.8'),45)]:
    th=thdeg*pi/180; a=mpf('1e-6')
    print(r_c,q,thdeg, mp.nstr(moment(r_c,q,a,th)/a,12), mp.nstr(-sin(th)**2*(2*r_c-q**2)/(r_c**2-2*r_c+q**2),12))
