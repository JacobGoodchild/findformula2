# Massive-particle escape from near-horizon extremal Kerr-Newman (a^2+Q^2=1), ZAMO isotropic emitter, equator
import numpy as np
from mpmath import mp, mpf, sqrt, pi, quad, asin
def P_sim(a,beta,d=1e-7,N=220):
    r=1+d; Del=d*d; A=(r*r+a*a)**2-Del*a*a; gpp=A/(r*r); om=a*(r*r+a*a-Del)/A; al=d*r/np.sqrt(A)
    num=-(r*r+a*a)*d*(2+d)-d*d
    g=1/np.sqrt(1-beta**2)
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij'); nr=MU; st=np.sqrt(1-MU**2); nth=st*np.cos(PH); nph=st*np.sin(PH)
    E=g*(al+om*np.sqrt(gpp)*beta*nph); L=g*beta*nph*np.sqrt(gpp); Q=(g*beta*nth)**2*r*r
    X0=g*((1+a*a)*al+beta*nph*np.sqrt(gpp)*a*num/A)
    dout=np.geomspace(d,1e4,4000); din=d*np.geomspace(1e-8,1,3000)[:-1]
    esc=0
    for i in range(N):
        for j in range(2*N):
            e=E[i,j]
            if e<1: continue
            l=L[i,j]; q=Q[i,j]; x0=X0[i,j]
            Ro=(x0+e*(2*dout+dout**2))**2-dout**2*((1+dout)**2+(l-a*e)**2+q)
            if not np.all(Ro>0): continue
            if nr[i,j]>0: esc+=1
            else:
                Ri=(x0+e*(2*din+din**2))**2-din**2*((1+din)**2+(l-a*e)**2+q)
                esc+=np.any(Ri<0)
    return esc/(2*N*N)
def P_formula(a,beta):
    a=mpf(a); beta=mpf(beta); k=1/(2*a)
    t0=sqrt(1-beta**2)/(a*beta); t1=k/beta
    if t0>=1: return mpf(0)
    if t0>=t1: return (1-t0)/2
    J=quad(lambda t: asin(sqrt(((1/k**2-1)*t*t+1-1/beta**2)/(1-t*t))),[t0,t1])
    return ((1-t1)/2+J/(2*pi)).real
if __name__=="__main__":
    for a,b in [(0.8,0.95),(0.8,0.85),(0.8,0.7),(0.6,0.9)]:
        print(a,b, P_sim(a,b), float(P_formula(a,b)), flush=True)
