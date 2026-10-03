import numpy as np
from mb_avg_sim import Lc
def moments(a,n=2000):
    rp=1+np.sqrt(max(1-a*a,0)); rr=rp+np.geomspace(1e-10,80,6000)
    mus=np.linspace(-1,1,n+1); mus=(mus[:-1]+mus[1:])/2
    L=np.array([Lc(a,m,rr) for m in mus])
    return np.sum(mus*L**3)*(2/n), np.sum(L**2)*(2/n)
if __name__=="__main__":
    T,S=moments(1.0); print('a=1', T, S, (2/3)*T/S)
    T,S=moments(0.5); print('a=0.5', T, S, (2/3)*T/S)
