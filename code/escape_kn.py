# Independent ray-by-ray check: escape probability, isotropic ZAMO emitter, extremal Kerr-Newman (M=1, a^2+e^2=1),
# at r*=1+eps, polar angle th.  Conjecture: P = 1/2 - asin(k)/(2pi) - k/4, k=(1+a^2 cos^2 th)/(2 a sin th)
import numpy as np, sys
def P_sim(a, th, eps, N=200):
    e2=1-a*a; r=1+eps; s=np.sin(th); c=np.cos(th)
    Sig=r*r+a*a*c*c; Delta=r*r-2*r+a*a+e2; A=(r*r+a*a)**2-Delta*a*a*s*s
    gpp=A*s*s/Sig; omega=a*(r*r+a*a-Delta)/A; alpha=np.sqrt(Sig*Delta/A)
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij')
    nr=MU; st=np.sqrt(1-MU**2); nth=st*np.cos(PH); nph=st*np.sin(PH)
    E=alpha+omega*np.sqrt(gpp)*nph; L=nph*np.sqrt(gpp)
    lam=L/E; Theta=(nth**2*Sig)/E**2
    eta=Theta-a*a*c*c+lam**2*c*c/(s*s)
    rr=1+np.geomspace(eps,1e4,4000); rin=1+eps*np.geomspace(1e-7,1,3000)[:-1]
    esc=0
    for i in range(N):
        for j in range(2*N):
            if E[i,j]<=0: continue
            l=lam[i,j]; q=eta[i,j]
            R=(rr*rr+a*a-a*l)**2-(rr*rr-2*rr+a*a+e2)*(q+(l-a)**2)
            if not np.all(R>0): continue
            if nr[i,j]>0: esc+=1
            else:
                Rin=(rin*rin+a*a-a*l)**2-(rin*rin-2*rin+a*a+e2)*(q+(l-a)**2)
                esc+=np.any(Rin<0)
    return esc/(2*N*N)
def P_formula(a,th):
    k=(1+a*a*np.cos(th)**2)/(2*a*np.sin(th))
    return 0.0 if k>=1 else 0.5-np.arcsin(k)/(2*np.pi)-k/4
if __name__=="__main__":
    for a,thdeg in [(1,90),(1,70),(0.8,90),(0.8,60),(0.6,90),(0.9,50)]:
        th=np.radians(thdeg)
        print(a,thdeg,[round(P_sim(a,th,eps),5) for eps in (1e-4,1e-6)], round(P_formula(a,th),5), flush=True)
