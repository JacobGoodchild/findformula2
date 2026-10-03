import numpy as np
def P_sim(beta,v=0.5,d=1e-7,N=220):
    a=1.0; r=1+d; Del=d*d; A=(r*r+1)**2-Del; gpp=A/(r*r); om=(r*r+1-Del)/A; al=d*r/np.sqrt(A)
    num=-(r*r+1)*d*(2+d)-d*d
    gb=1/np.sqrt(1-beta**2); gv=1/np.sqrt(1-v*v)
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij'); ner=MU; st=np.sqrt(1-MU**2); neth=st*np.cos(PH); neph=st*np.sin(PH)
    ut=gv*gb*(1+v*beta*neph); uph=gv*gb*(beta*neph+v); ur=gb*beta*ner; uth=gb*beta*neth
    E=ut*al+om*np.sqrt(gpp)*uph; L=np.sqrt(gpp)*uph; Q=uth**2*r*r
    X0=2*ut*al+uph*np.sqrt(gpp)*num/A     # 2E - L  (a=1)
    dout=np.geomspace(d,1e4,4000); din=d*np.geomspace(1e-8,1,3000)[:-1]
    esc=0
    for i in range(N):
        for j in range(2*N):
            e=E[i,j]
            if e<1: continue
            l=L[i,j]; q=Q[i,j]; x0=X0[i,j]
            Ro=(x0+e*(2*dout+dout**2))**2-dout**2*((1+dout)**2+(l-e)**2+q)
            if not np.all(Ro>0): continue
            if ur[i,j]>0: esc+=1
            else:
                Ri=(x0+e*(2*din+din**2))**2-din**2*((1+din)**2+(l-e)**2+q)
                esc+=np.any(Ri<0)
    return esc/(2*N*N)
if __name__=="__main__":
    from isco_massive import P
    for b in [0.4,0.5,0.7,0.9]:
        print(b, P_sim(b), float(P(b).real), flush=True)
