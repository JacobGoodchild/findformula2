# Brute force: massive particles, speed beta, isotropic in ZAMO frame at r=1+d, extremal Kerr equator
import numpy as np
def P_sim(beta,d,N=250):
    a=1.0; r=1+d; Del=d*d; A=(r*r+1)**2-Del; gpp=A/(r*r); om=(r*r+1-Del)/A; al=d*r/np.sqrt(A)
    num=-(r*r+1)*d*(2+d)-d*d      # (1+a^2)(r^2+a^2-Del)-A  -> used for 2w-1
    g=1/np.sqrt(1-beta**2)
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij')
    nr=MU; st=np.sqrt(1-MU**2); nth=st*np.cos(PH); nph=st*np.sin(PH)
    E=g*(al+om*np.sqrt(gpp)*beta*nph); L=g*beta*nph*np.sqrt(gpp); Q=(g*beta*nth)**2*r*r
    # X0 = 2E - L computed cancellation-free: 2E-L = g(2 al + beta nph sqrt(gpp)(2 om - 1)), 2om-1 = num/A
    X0=g*(2*al+beta*nph*np.sqrt(gpp)*num/A)
    dout=np.geomspace(d,1e4,4000); din=d*np.geomspace(1e-8,1,3000)[:-1]
    esc=0
    for i in range(N):
        for j in range(2*N):
            e=E[i,j]
            if e<1: continue
            l=L[i,j]; q=Q[i,j]; x0=X0[i,j]
            # X(r) = E(r^2+1)-L = x0 + E(2 dd + dd^2)
            Ro=(x0+e*(2*dout+dout**2))**2-dout**2*((1+dout)**2+(l-e)**2+q)
            if not np.all(Ro>0): continue
            if nr[i,j]>0: esc+=1
            else:
                Ri=(x0+e*(2*din+din**2))**2-din**2*((1+din)**2+(l-e)**2+q)
                esc+=np.any(Ri<0)
    return esc/(2*N*N)
if __name__=="__main__":
    from massive_escape import P
    for b in [0.95,0.8,0.75]:
        print(b, P_sim(b,1e-7), float(P(b).real), flush=True)
