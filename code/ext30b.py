from mpmath import *
from ext_closed import A_closed
mp.dps=60
A=A_closed(pi/6)
m=(2-sqrt(3))/4; K=ellipk(m); E=ellipe(m); P=ellippi(mpf(1)/2,m)
r3=sqrt(3)
for basis,names in [([A,K,r3*K,E,r3*E,P,r3*P],'A K r3K E r3E P r3P'),([A,K,r3*K,E,r3*E,pi,r3*pi],'A K r3K E r3E pi r3pi'),([A,K,r3*K,pi/K,r3*pi/K,P,r3*P,pi,r3*pi],'A K r3K pi/K r3pi/K P r3P pi r3pi')]:
    print(names, pslq(basis,maxcoeff=10**6,maxsteps=10**6))
