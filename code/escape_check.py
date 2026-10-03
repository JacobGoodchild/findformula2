# Independent check: escape probability for isotropic emitter at rest in ZAMO frame, extremal Kerr, equator, r*=1+eps
import numpy as np
from scipy.optimize import brentq
a=1.0
def escape_prob(rs, N=400, frame='ZAMO'):
    r=rs; Sig=r*r; Delta=r*r-2*r+a*a; A=(r*r+a*a)**2-Delta*a*a
    gpp=A/Sig; omega=2*a*r/A; alpha=np.sqrt(Sig*Delta/A)
    # grid over directions: mu=cos(angle from outward radial), phi azimuth in (theta,phi) plane
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij')
    nr=MU; st=np.sqrt(1-MU**2); nth=st*np.cos(PH); nph=st*np.sin(PH)
    if frame=='Carter':
        # boost from ZAMO along phi with velocity v of Carter observer
        OmC=a/(r*r+a*a); v=(OmC-omega)*np.sqrt(gpp)/alpha
        g=1/np.sqrt(1-v*v)
        # photon in Carter frame (E'=1, n'), transform to ZAMO: E=g(1+v nph'), nph=(nph'+v)/(1+v nph'), others /(g(1+v nph'))
        D=g*(1+v*nph); EZ=D; nph2=(nph+v)/(1+v*nph); nr2=nr/D; nth2=nth/D
        nr,nth,nph=nr2,nth2,nph2
    else:
        EZ=np.ones_like(nr)
    E=EZ*(alpha+omega*np.sqrt(gpp)*nph); L=EZ*nph*np.sqrt(gpp)
    lam=L/E; eta=(EZ*nth*np.sqrt(Sig))**2/E**2
    esc=np.zeros_like(nr,dtype=bool)
    rr=1+np.geomspace(rs-1,1e4,6000)
    rin=1+(rs-1)*np.geomspace(1e-6,1,3000)[:-1]
    for i in range(nr.shape[0]):
        for j in range(nr.shape[1]):
            l=lam[i,j]; e=eta[i,j]
            R=(rr*rr+a*a-a*l)**2-(rr*rr-2*rr+a*a)*(e+(l-a)**2)
            if not np.all(R>0): continue
            if nr[i,j]>0: esc[i,j]=True
            else:
                Rin=(rin*rin+a*a-a*l)**2-(rin*rin-2*rin+a*a)*(e+(l-a)**2)
                esc[i,j]=np.any(Rin<0)
    return esc.mean()
if __name__=="__main__":
    import sys
    for eps in [1e-2,1e-3,1e-4,1e-5]:
        print(eps, escape_prob(1+eps,N=150), escape_prob(1+eps,N=150,frame='Carter'), 7/24)
