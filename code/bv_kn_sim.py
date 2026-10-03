import numpy as np
from escape_isco_sim import P_sim
from bv_kn import Omega_kn
def v_carter(a,R):
    dl=R-1; Del=dl*dl; A=(R*R+a*a)**2-Del*a*a
    om=a*(R*R+a*a-Del)/A; OmC=a/(R*R+a*a)
    return (OmC-om)*A/(R*R*dl)
for a,z in [(0.8,0.5),(0.8,0.2),(0.4,0.5),(1.0,0.5)]:
    R=1/z; v=v_carter(a,R)
    P=P_sim(a,R-1,v,N=200)
    print(a,z,'v=',round(v,4),' sim',4*np.pi*(1-P),' formula',float(Omega_kn(z,a)),flush=True)
