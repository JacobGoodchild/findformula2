from mpmath import *
mp.dps=40
def roots(th):
    s=sin(th)
    return (1-s-sqrt(2-2*s), 1+s-sqrt(2+2*s), 1-s+sqrt(2-2*s), 1+s+sqrt(2+2*s))
def A_new(th):
    s=sin(th); c2=cos(th)**2
    rA,rB,rC,rD=roots(th)
    R=lambda r:(r-rA)*(r-rB)*(r-rC)*(rD-r)
    lo=max(rC,mpf(1))
    I=4*quad(lambda r:(3*r*r-c2)/sqrt(R(r)),[lo,(lo+rD)/2,rD])
    B=0 if rC>=1 else (2+s*s)*sqrt(3*s*s+s*s*c2-4*c2)/s**2
    return I+B
def A_direct(th):
    s=sin(th); c2=cos(th)**2
    Q=lambda r: r**3*(4-r)+c2-(c2/(1-c2))*(r*r-2*r-1)**2
    rA,rB,rC,rD=roots(th); lo=max(rC,mpf(1))
    return 4/s*quad(lambda r:(r-1)*sqrt(Q(r)),[lo,(lo+rD)/2,rD])
if __name__=="__main__":
  for deg in ['0.5','10','30','45','47.06','50','60','80','90']:
    th=mpf(deg)*pi/180
    print(deg, nstr(A_new(th),30), nstr(A_direct(th),30))
  print('polar', pi*(12+8*sqrt(2)))
