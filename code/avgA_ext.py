from mpmath import *
import ext_closed, ext_formula
mp.dps=45
def f(x):
    if x<mpf('1e-40'): return 16*pi+15*sqrt(3)
    if x==1: return pi*(12+8*sqrt(2))
    if x<mpf('0.05'): return ext_formula.A_new(acos(x)).real
    return ext_closed.A_closed(acos(x))
xc=sqrt(2*sqrt(3)-3)
v=quad(f,[0,mpf('0.05'),xc,1])
print(nstr(v,40))
