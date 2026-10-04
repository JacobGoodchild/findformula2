from mpmath import mp, mpf, pi, cos, matrix, lu_solve, nstr, sqrt
from univ_check_light import make
mp.dps=40
def K4(r,C): return (C**2*(r**6-38*r**5+286*r**4-768*r**3+774*r**2-108*r-162)+C*(-18*r**6+108*r**5-288*r**4+492*r**3-504*r**2+216*r)-15*r**6+58*r**5-78*r**4+36*r**3)/(2*r**4*(r-1)**2*(2*r-3)**5)
for q,thdeg in [(mpf('0.6'),60),(mpf('0.8'),90)]:
    fun=lambda x,q=q: 1-2/x+q**2/x**2
    ar=make(fun); th=thdeg*pi/180; C=cos(th)**2
    hs=[mpf(k)/100 for k in range(1,6)]; ys=[]
    for h in hs:
        A,rph=ar(h,th); A0=pi*rph**2/fun(rph); ys.append((A/A0-1)/h**2)
    c=lu_solve(matrix([[h**(2*j) for j in range(5)] for h in hs]),matrix(ys))
    print('KN q=%s %d deg: a^4 numeric %s  formula %s'%(q,thdeg,nstr(c[1],12),nstr(K4(rph,C),12)),flush=True)
