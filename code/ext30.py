from mpmath import *
mp.dps=60
def Aext(th):
    s=sin(th)
    P=lambda x: -((x+s)**2-2+2*s)*((x-s)**2-2-2*s)
    xD=s+sqrt(2+2*s); xC=-s+sqrt(2-2*s)
    lo=max(xC,0)
    return 4/s**2*quad(lambda x:x*sqrt(P(x)),[lo,(lo+xD)/2,xD])
A=Aext(pi/6).real; print(A)
k2=(2-sqrt(3))/4; K=ellipk(k2); E=ellipe(k2)
print(K, mpf(3)**0.25*gamma(mpf(1)/3)**3/(2**(mpf(7)/3)*pi))
basis=[A,pi,sqrt(3)*pi,K,sqrt(3)*K,E,sqrt(3)*E,1/K, sqrt(3)/K]
print(pslq(basis,maxcoeff=10**5,maxsteps=10**6))
