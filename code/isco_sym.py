import sympy as sp
t,k=sp.symbols('t k',positive=True)
# J = -Int_{-k}^0 t(1+kt)/((1-t^2) sqrt(-t(2k+(1+k^2)t))) dt ; substitute t=-k*u... use t = -2k/(1+k^2) * sin^2? do numerically-guided: substitute t=-s, s in [0,k]
s=sp.symbols('s',positive=True)
integrand=-((-s)*(1-k*s)/((1-s**2)*sp.sqrt(s*(2*k-(1+k**2)*s))))   # dt=-ds, limits flip -> Int_0^k
# integrand above already includes the outer minus and dt=-ds with flipped limits
# split rational part: s(1-ks)/(1-s^2) = k + (s - k)/(1-s^2)... let sympy apart
R=sp.apart(s*(1-k*s)/(1-s**2),s)
print('apart:',R)
w=sp.symbols('w',positive=True)
F=sp.integrate(s*(1-k*s)/((1-s**2)*sp.sqrt(s*(2*k-(1+k**2)*s))),(s,0,k))
print(sp.simplify(F))
