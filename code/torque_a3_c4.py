from torque_a3_mid import *
from fractions import Fraction
for rc in [Fraction(7,2),Fraction(10,3),Fraction(18,5)]:
    for C in [Fraction(1,4),Fraction(4,5)]:
        r=mpf(rc.numerator)/rc.denominator; E=sqrt((r-2)**2/(r*(r-3)))
        Cm=mpf(C.numerator)/C.denominator
        c1,c3=odd_coeffs(E,acos(sqrt(Cm)))
        pred=(1-Cm)*(mpf(3)*(r-5)/2-Cm*(7*r**2-63*r+144)/(2*r))/((r-2)*(6-r)**3)
        print(rc,C,nstr(c3,18),nstr(pred,18),nstr(c3-pred,3),flush=True)
