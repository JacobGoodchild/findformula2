from mpmath import mp, mpf, pi, matrix, lu_solve, nstr, diff
from univ_check_mass import area
mp.dps=40
fun=lambda x: 1-2*x**2/(x**2+mpf('0.16'))**1.5
for r_c in [mpf('3.4'),mpf('3.1')]:
    hs=[mpf(k)/100 for k in range(1,6)]; ys=[]
    for h in hs:
        A,A0=area(fun,r_c,h,mpf(0),n=8); ys.append((A/A0-1)/h**2)
    c=lu_solve(matrix([[h**(2*j) for j in range(5)] for h in hs]),matrix(ys))
    f=fun(r_c); g=diff(fun,r_c); h2=diff(fun,r_c,2); r=r_c; J=3*f*g+r*f*h2-2*r*g**2
    print('Bardeen polar r_c=%s: a^2 %s vs %s | a^4 %s vs %s'%(r_c,nstr(c[0],12),nstr((2*f+r*g-2)/(r**3*g),12),nstr(c[1],12),nstr(2*(1-f)**2*(3*g+r*h2)/(r**5*g*J),12)),flush=True)
