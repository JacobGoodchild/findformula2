# Cancellation-free ray-by-ray check (extremal Kerr-Newman, M=1, a^2+e^2=1, ZAMO isotropic emitter at r=1+delta)
import numpy as np
def P_sim(a, th, dl, N=300):
    s=np.sin(th); c=np.cos(th); r=1+dl
    Sig=r*r+a*a*c*c; Del=dl*dl; A=(r*r+a*a)**2-Del*a*a*s*s
    gpp=A*s*s/Sig; om=a*(r*r+a*a-Del)/A; al=dl*np.sqrt(Sig/A)
    num=-(r*r+a*a)*dl*(2+dl)-dl*dl*(1+a*a*c*c)          # (1+a^2)(r^2+a^2-Del)-A, exact
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij')
    nr=MU; st=np.sqrt(1-MU**2); nth=st*np.cos(PH); nph=st*np.sin(PH)
    E=al+om*np.sqrt(gpp)*nph; L=nph*np.sqrt(gpp)
    D0=((1+a*a)*al+nph*np.sqrt(gpp)*a*num/A)/E          # = 1+a^2-a*lambda
    lam=L/E; eta=nth**2*Sig/E**2-a*a*c*c+lam**2*c*c/(s*s)
    C=eta+(lam-a)**2
    dout=np.geomspace(dl,1e4,4000); din=dl*np.geomspace(1e-8,1,3000)[:-1]
    esc=0
    for i in range(N):
        for j in range(2*N):
            if E[i,j]<=0: continue
            Ro=(dout*dout+2*dout+D0[i,j])**2-dout*dout*C[i,j]
            if not np.all(Ro>0): continue
            if nr[i,j]>0: esc+=E[i,j]
            else:
                Ri=(din*din+2*din+D0[i,j])**2-din*din*C[i,j]
                esc+=E[i,j]*np.any(Ri<0)
    return esc/(2*N*N)
def P_formula(a,th):
    k=(1+a*a*np.cos(th)**2)/(2*a*np.sin(th))
    return 0.0 if k>=1 else 0.5-np.arcsin(k)/(2*np.pi)-k/4
if __name__=="__main__":
    for a in [1.0,0.8]:
        k=1/(2*a); print('zamo',a, P_sim(a,np.pi/2,1e-7,N=250), (4*a*a-1)/(32*a)+np.sqrt(4*a*a-1)/16, flush=True)
