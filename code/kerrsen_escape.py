# Near-horizon escape, extremal Kerr-Sen (a=1-b), ZAMO isotropic emitter at equator: brute force vs formula
import numpy as np
def P_sim(b,d=1e-7,N=220):
    a=1-b; rH=1-b; r=rH+d
    rho=r*(r+2*b); Del=(r-rH)**2; A=(rho+a*a)**2-Del*a*a
    gpp=A/rho; om=a*(rho+a*a-Del)/A; al=np.sqrt(rho*Del/A)
    # 1 + ... : X(r) = rho + a^2 - a*lambda ; need rho_H + a^2 - a*lambda accurately: rho_H + a^2 = 2a (so lambda_H = 2)
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij'); nr=MU; st=np.sqrt(1-MU**2); nth=st*np.cos(PH); nph=st*np.sin(PH)
    E=al+om*np.sqrt(gpp)*nph; L=nph*np.sqrt(gpp)
    # D0 = (rho+a^2) - a*lambda at the emitter = ((rho+a^2)E - aL)/E ; (rho+a^2) om - a = a(... )/A computed exactly:
    numer=(rho+a*a)*(rho+a*a-Del)-A           # = Del*a^2 - (rho+a^2)Del = -Del*rho  (exact)
    numer=-Del*rho
    D0=((rho+a*a)*al+nph*np.sqrt(gpp)*a*numer/A)/E
    lam=L/E; eta=(nth**2*rho)/E**2; C=eta+(lam-a)**2
    dout=np.geomspace(d,1e4,4000); din=d*np.geomspace(1e-8,1,3000)[:-1]
    def X(dd): return D0 + ((r+dd)*(r+dd+2*b)-rho)      # rho(r+dd)+a^2-a lam
    esc=0
    for i in range(N):
        for j in range(2*N):
            if E[i,j]<=0: continue
            Ro=(D0[i,j]+((r+dout)*(r+dout+2*b)-rho))**2-(d+dout)**2*C[i,j]
            if not np.all(Ro>0): continue
            if nr[i,j]>0: esc+=1
            else:
                Ri=(D0[i,j]+((r-din)*(r-din+2*b)-rho))**2-(d-din)**2*C[i,j]
                esc+=np.any(Ri<0)
    return esc/(2*N*N)
def P_formula(b):
    k=(1+b)/2; return 0.5-np.arcsin(k)/(2*np.pi)-k/4
if __name__=="__main__":
    for b in [0.0,0.3,0.6]:
        print(b, P_sim(b), P_formula(b), flush=True)
