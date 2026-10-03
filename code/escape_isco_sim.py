# Brute-force check: isotropic emitter on near-horizon circular orbit (speed v rel. to ZAMO, prograde), extremal KN, equator
import numpy as np
def P_sim(a, dl, v, N=300):
    r=1+dl; c=0.0; s=1.0
    Sig=r*r; Del=dl*dl; A=(r*r+a*a)**2-Del*a*a
    gpp=A/Sig; om=a*(r*r+a*a-Del)/A; al=dl*np.sqrt(Sig/A)
    num=-(r*r+a*a)*dl*(2+dl)-dl*dl*1.0
    g=1/np.sqrt(1-v*v)
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij')
    # emitter-frame direction
    t=MU; st=np.sqrt(1-t*t); nr_e=st*np.cos(PH); nth_e=st*np.sin(PH); nph_e=t
    D=g*(1+v*nph_e)                     # ZAMO energy (emitter energy=1)
    nph=(nph_e+v)/(1+v*nph_e); nr=nr_e/D; nth=nth_e/D
    EZ=D
    E=EZ*(al+om*np.sqrt(gpp)*nph); L=EZ*nph*np.sqrt(gpp)
    D0=(EZ*((1+a*a)*al+nph*np.sqrt(gpp)*a*num/A))/E
    lam=L/E; eta=(EZ*nth)**2*Sig/E**2
    C=eta+(lam-a)**2
    dout=np.geomspace(dl,1e4,4000); din=dl*np.geomspace(1e-8,1,3000)[:-1]
    esc=0
    for i in range(N):
        for j in range(2*N):
            if E[i,j]<=0: continue
            Ro=(dout*dout+2*dout+D0[i,j])**2-dout*dout*C[i,j]
            if not np.all(Ro>0): continue
            if nr[i,j]>0: esc+=1
            else:
                Ri=(din*din+2*din+D0[i,j])**2-din*din*C[i,j]
                esc+=np.any(Ri<0)
    return esc/(2*N*N)
def P_formula(k):
    return 0.5-np.arcsin(k)/(2*np.pi)+k*np.arctan(np.sqrt((1+k*k)/(1-k*k)))/(np.pi*np.sqrt(1+k*k))
if __name__=="__main__":
    for a in [1.0,0.8,1/np.sqrt(2)]:
        k=1/(2*a)
        print(a,[round(P_sim(a,d,k),5) for d in (1e-7,)], 'formula', round(P_formula(k),5), flush=True)
