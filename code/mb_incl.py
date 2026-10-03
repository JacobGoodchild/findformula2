from mpmath import mp, mpf, sqrt, quad, pi, cos, sin, polyroots, matrix, lu_solve, nstr
import sympy as sp
mp.dps=50
_u,_a,_c=sp.symbols('u a c')
_y=1+_u
_N=sp.Poly(sp.expand((1-_c**2)*_y**4*(_a**2-(_y**2-2*_y)**2)-_c**2*(_y**4-2*_y**3+_a**2)**2),_u)
_coeffs=[sp.lambdify((_a,_c),co,'mpmath') for co in _N.all_coeffs()]
def A_mb(a,th):
    a=mpf(a); c=cos(th); s=sin(th)
    co=[f(a,c) for f in _coeffs]
    N=lambda u: sum(co[i]*u**(len(co)-1-i) for i in range(len(co)))
    rts=sorted([x.real for x in polyroots(co,maxsteps=800,extraprec=600) if abs(x.imag)<mpf(10)**-25 and x.real>0])
    lo,hi=rts[-2],rts[-1]
    g=lambda u: sqrt(N(u))*abs(3*u**4+4*u**3+1-a*a)/(a*a*s*u**3)
    f=lambda t:(lambda u: g(u)*(hi-lo)*sin(t)/2)(lo+(hi-lo)*(1-cos(t))/2)
    return 2/s*quad(f,[0,pi/2,pi]).real
if __name__=="__main__":
    print(A_mb('0.5',pi/2), 7*pi+pi*sqrt(1-mpf('0.25'))+16*sqrt(mpf('1.5'))*__import__('mpmath').ellipe(mpf(2)/3))
    print(A_mb('0.5',pi/3), A_mb('0.5',pi/180))
