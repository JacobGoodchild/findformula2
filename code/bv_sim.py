import numpy as np
from escape_check import escape_prob
from mpmath import mpf
from bv_equator import Omega_eq
for z in [0.5,0.3,0.1]:
    rs=1/z
    Pc=escape_prob(rs,N=200,frame='Carter'); Pz=escape_prob(rs,N=200,frame='ZAMO')
    print(z, 'sim Carter 4pi(1-P)=',4*np.pi*(1-Pc),' ZAMO:',4*np.pi*(1-Pz),' formula:',float(Omega_eq(z)[0]),flush=True)
