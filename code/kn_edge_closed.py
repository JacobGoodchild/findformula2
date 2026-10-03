from mpmath import mp, mpf, sqrt, quad, ellipk, ellipe, ellippi, polyroots, pi
import numpy as np
mp.dps=30
def Wcoeffs(a,q):  # W = 4a^2(r-q) - (r^2-3r+2q)^2
    return [-1, 6, -(9+4*q), 12*q+4*a*a, -(4*q*q+4*a*a*q)]
def A_closed(a,q):
    a=mpf(a); q=mpf(q)
    rts=sorted([x.real for x in polyroots(Wcoeffs(a,q),maxsteps=500,extraprec=300)])
    d,c,b,A1=rts
    m=(A1-b)*(c-d)/((A1-c)*(b-d)); n=(A1-b)/(A1-c); g=2/sqrt((A1-c)*(b-d))
    K,E,P=ellipk(m),ellipe(m),ellippi(n,m)
    V=(n*E+(m-n)*K+(2*n*m+2*n-n*n-3*m)*P)/(2*(n-1)*(m-n))
    I0=g*K; I1=g*(c*K+(b-c)*P); I2=g*(c*c*K+2*c*(b-c)*P+(b-c)**2*V)
    n2=n*(c-1)/(b-1); J=g*(K/(c-1)+(c-b)/((b-1)*(c-1))*ellippi(n2,m))
    c0=2*(4*q*q-2*q-3)/(q-1); c1=-(8*q*q+10*q-21)/(q-1); c2=(14*q-15)/(q-1); d1=-2*(4*q-3)*(a*a+q-1)/(q-1)
    return c0*I0+c1*I1+c2*I2+d1*J, rts
def A_shoelace(a,q,N=1000000):
    a=float(a); q=float(q)
    rts=np.roots([float(x) for x in Wcoeffs(a,q)]).real; rts.sort(); lo,hi=rts[2],rts[3]
    t=np.linspace(0,np.pi,N); r=lo+(hi-lo)*(1-np.cos(t))/2
    xi=-(r**3-3*r**2+a*a*r+a*a+2*q*r)/(a*(r-1)); eta=r**2*(4*a*a*(r-q)-(r*r-3*r+2*q)**2)/(a*a*(r-1)**2)
    al=-xi; be=np.sqrt(np.clip(eta,0,None))
    X=np.concatenate([al,al[::-1]]); Y=np.concatenate([be,-be[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
if __name__=="__main__":
    for a,q in [(0.5,0.0),(0.5,0.3),(0.3,0.6),(0.8,0.2),(0.9,0.15),(0.1,0.9)]:
        A,rts=A_closed(a,q); print(a,q, A, A_shoelace(a,q), [float(x) for x in rts])

def A_reduced(a,q):
    a=mpf(a); q=mpf(q)
    co=Wcoeffs(a,q); W=lambda r: sum(co[i]*r**(4-i) for i in range(5))
    rts=sorted([x.real for x in polyroots(co,maxsteps=500,extraprec=300) if abs(x.imag)<mpf(10)**-20])
    lo,hi=rts[-2],rts[-1]
    c0=2*(4*q*q-2*q-3)/(q-1); c1=-(8*q*q+10*q-21)/(q-1); c2=(14*q-15)/(q-1); d1=-2*(4*q-3)*(a*a+q-1)/(q-1)
    from mpmath import cos as mc, sin as ms
    f=lambda t: (lambda r: (c0+c1*r+c2*r*r+d1/(r-1))/sqrt(W(r))*(hi-lo)*ms(t)/2)(lo+(hi-lo)*(1-mc(t))/2)
    return quad(f,[0,pi/2,pi])
