import sympy as sp
from mpmath import mp, mpf, sqrt, quad, pi, findroot, diff, sin, cos, matrix, lu_solve, nstr
expr=open('univ_a4_fast.out').read().split('a^4 relative (R=1):')[1].split('\n')[0]
Cs=sp.Symbol('C'); fs=sp.symbols('f0:8')
E4=sp.sympify(expr,locals={'C':Cs,**{'f%d'%i:fs[i] for i in range(8)}})
def K4(fun,R,C):
    vals={fs[k]:diff(fun,R,k)*R**k for k in range(6)}
    return sp.N(E4.subs(vals).subs(Cs,C),30)/R**4
# Kerr exact
r=sp.Symbol('r'); fK=1-2/r
kk=sp.nsimplify(sp.simplify(E4.subs({fs[k]:sp.diff(fK,r,k).subs(r,3)*3**k for k in range(6)})/81))
print('Kerr a^4:',sp.expand(kk),'  expected -(1/72 + 5C/972 - 5C^2/1944) relative to 27pi')
mp.dps=40
from univ_check_light import make
fun=lambda x: 1-2*x**2/(x**2+mpf('0.16'))**1.5
ar=make(fun)
for thdeg in [90,50]:
    th=thdeg*pi/180; C=cos(th)**2
    hs=[mpf(k)/100 for k in range(1,6)]; ys=[]
    for h in hs:
        A,rph=ar(h,th); A0=pi*rph**2/fun(rph); ys.append((A/A0-1)/h**2)
    c=lu_solve(matrix([[h**(2*j) for j in range(5)] for h in hs]),matrix(ys))
    print('Bardeen',thdeg,'a^4 numeric',nstr(c[1],12),' formula',nstr(K4(fun,rph,C),12),flush=True)
