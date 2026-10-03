import numpy as np
from escape_kn2 import P_sim, P_formula
for a,thdeg in [(0.9,50),(0.8,60)]:
    th=np.radians(thdeg)
    print(a,thdeg,[ (N,round(P_sim(a,th,1e-7,N=N),5)) for N in (150,300,600)],'formula',round(P_formula(a,th),5),flush=True)
