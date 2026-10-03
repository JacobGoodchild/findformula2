import numpy as np
def captured(a,Lz,Q,rr):
    R=(rr**2+a*a-a*Lz)**2-(rr**2-2*rr+a*a)*(rr**2+(Lz-a)**2+Q)
    return np.all(R>0)
def Lc(a,mu,rr):
    lo,hi=0.0,8.0
    for _ in range(50):
        m=(lo+hi)/2
        if captured(a,mu*m,(1-mu*mu)*m*m,rr): lo=m
        else: hi=m
    return lo
def avg_sigma(a,n=2000):
    rp=1+np.sqrt(max(1-a*a,0)); rr=rp+np.geomspace(1e-10,80,6000)
    mus=np.linspace(-1,1,n+1); mus=(mus[:-1]+mus[1:])/2
    return np.pi/2*sum(Lc(a,m,rr)**2 for m in mus)*(2/n)
if __name__=="__main__":
    import math
    exact=math.pi*(47/9+31*math.sqrt(2)/6+7*math.sqrt(6)/9-math.sqrt(3)/54*math.log((1+math.sqrt(2)+math.sqrt(6))/(1+math.sqrt(2))))
    print('a=1  sim',avg_sigma(1.0),' exact',exact)
    print('a=0.5 sim',avg_sigma(0.5),' integral 49.4330447022985735715671195891')
