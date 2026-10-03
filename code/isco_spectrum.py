import numpy as np
from escape_isco_sim import P_sim  # reuse geometry by copying logic with energy histogram
def spectrum_sim(a=1.0,d=1e-7,v=0.5,N=300):
    import numpy as np
    r=1+d; Sig=r*r; Del=d*d; A=(r*r+a*a)**2-Del*a*a
    gpp=A/Sig; om=a*(r*r+a*a-Del)/A; al=d*np.sqrt(Sig/A)
    num=-(r*r+a*a)*d*(2+d)-d*d
    g=1/np.sqrt(1-v*v)
    mus=(np.arange(N)+0.5)/N*2-1; phs=(np.arange(2*N)+0.5)/(2*N)*2*np.pi
    MU,PH=np.meshgrid(mus,phs,indexing='ij')
    t=MU; st=np.sqrt(1-MU**2); nr_e=st*np.cos(PH); nth_e=st*np.sin(PH); nph_e=t
    D=g*(1+v*nph_e); nph=(nph_e+v)/(1+v*nph_e); nr=nr_e/D; nth=nth_e/D
    EZ=D; E=EZ*(al+om*np.sqrt(gpp)*nph); L=EZ*nph*np.sqrt(gpp)
    D0=(EZ*((1+a*a)*al+nph*np.sqrt(gpp)*a*num/A))/E
    lam=L/E; eta=(EZ*nth)**2*Sig/E**2; C=eta+(lam-a)**2
    dout=np.geomspace(d,1e4,4000); din=d*np.geomspace(1e-8,1,3000)[:-1]
    gs=[]
    for i in range(N):
        for j in range(2*N):
            if E[i,j]<=0: continue
            Ro=(dout*dout+2*dout+D0[i,j])**2-dout*dout*C[i,j]
            if not np.all(Ro>0): continue
            ok=nr[i,j]>0
            if not ok:
                Ri=(din*din+2*din+D0[i,j])**2-din*din*C[i,j]; ok=np.any(Ri<0)
            if ok: gs.append(E[i,j])
    return np.array(gs), 2*N*N
gs,tot=spectrum_sim()
bins=np.linspace(0,1.8,19)
h,_=np.histogram(gs,bins)
def pred(g):
    if g>=1/np.sqrt(3) and g<=np.sqrt(3): return np.sqrt(3)/4
    if 0<g<1/np.sqrt(3):
        q=2*np.sqrt(3)*g/np.sqrt(3+2*np.sqrt(3)*g-3*g*g); return np.sqrt(3)/(4*np.pi)*np.arcsin(min(q,1))
    return 0
for k in range(len(h)):
    c=(bins[k]+bins[k+1])/2
    print(round(c,3), round(h[k]/tot/(bins[1]-bins[0]),4), round(pred(c),4))
print('mean g', gs.sum()/tot)
