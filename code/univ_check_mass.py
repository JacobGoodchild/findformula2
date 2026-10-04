import sympy as sp
from mpmath import mp, mpf, sqrt, cos, sin, pi, findroot, diff, matrix, lu_solve, nstr
r,q,C=sp.symbols('r_c q C')
f0,f1,f2=sp.symbols('f0 f1 f2')
tot=sp.sympify(open('univ_a2_mass.out').read().split('TOTAL relative a^2 coefficient:')[1].split('\n')[0])
fK=1-2/r+q**2/r**2
kn=tot.subs({f0:fK,f1:sp.diff(fK,r),f2:sp.diff(fK,r,2)})
K31=-(r**2*(r-q**2)-C*(3*r**3-12*r**2-q**2*r**2+18*q**2*r-8*q**4))/(2*r*(r-q**2)*(6*r**2-r**3-9*q**2*r+4*q**4))
print('KN reduction minus [31]:', sp.simplify(kn-K31))
Kfun=sp.lambdify((f0,f1,f2,r,C),tot,'mpmath')
mp.dps=40
def area(fun,r_c,a,th,n=48):
    F0=fun(r_c); F1=diff(fun,r_c)
    E=sqrt(2*F0**2/(2*F0-r_c*F1)); L0=sqrt(F1*r_c**3/(2*F0-r_c*F1)); p=sqrt(E*E-1); b0=L0/p
    st,ct=sin(th),cos(th); tot=0
    for k in range(n):
        ph=2*pi*(k+mpf(1)/2)/n
        def eqs(rr,b):
            al=b*cos(ph); be=b*sin(ph); Lz=p*al*st; Q=p**2*(be**2+(al**2-a**2)*ct**2)
            D=lambda x: x*x*fun(x)+a*a
            R=lambda x:(E*(x*x+a*a)-a*Lz)**2-D(x)*(x*x+(Lz-a*E)**2+Q)
            return [R(rr),diff(R,rr)]
        rr,bb=findroot(eqs,(r_c,b0)); tot+=bb**2
    return pi*tot/n, pi*b0**2
if __name__!="__main__": tests={}
else: tests={'Bardeen g=0.4':lambda x: 1-2*x**2/(x**2+mpf('0.16'))**1.5,'Hayward l=0.5':lambda x: 1-2*x**2/(x**3+mpf('0.5'))}
for name,fun in tests.items():
    for r_c,thdeg in [(mpf('3.4'),60),(mpf('3.1'),90)]:
        th=thdeg*pi/180; hs=[mpf('0.004'),mpf('0.008'),mpf('0.012')]; ys=[]
        for h in hs:
            A,A0=area(fun,r_c,h,th); ys.append((A/A0-1)/h**2)
        K=lu_solve(matrix([[1,h**2,h**4] for h in hs]),matrix(ys))[0]
        print(name,r_c,thdeg,'numeric',nstr(K,13),'formula',nstr(Kfun(fun(r_c),diff(fun,r_c),diff(fun,r_c,2),r_c,cos(th)**2),13),flush=True)
