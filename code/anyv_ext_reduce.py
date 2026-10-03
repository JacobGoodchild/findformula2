import sympy as sp
t,E=sp.symbols('t E',positive=True)
p2=E**2-1
r=1/(t**2-p2); w=r*t
L=E*(r+1)-(r-1)*w
N4=sp.expand(t**4+2*E*t**3-t**2-2*E*p2*t-p2**2)
dL=sp.factor(sp.diff(L,t))
print('dL/dt =',dL)
# area integrand (before 2/p^2): sqrt(N4)/(t^2-p2)^2 * (+/- dL)  = R / sqrt(N4)
R=sp.together(N4*dL/(t**2-p2)**2)
print('R denominator:', sp.factor(sp.denom(R)), ' numerator degree', sp.degree(sp.numer(R),t))
