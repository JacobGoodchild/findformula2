import numpy as np
from escape_kn import P_sim, P_formula
for a,thdeg in [(1,90),(0.8,60),(0.9,50)]:
    th=np.radians(thdeg)
    print(a,thdeg,'N=400 eps=1e-7:',round(P_sim(a,th,1e-7,N=400),5),'formula',round(P_formula(a,th),5),flush=True)
