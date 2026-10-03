from mpmath import mp, mpf, sqrt, quad, ellipk, ellipe, ellippi, polyroots, pi
mp.dps=40
def rts(a):
    rs=sorted([x.real for x in polyroots([1,-6,9,-4*a**2],maxsteps=200,extraprec=200)])
    return rs  # r3 < rpro < rretro
def A_closed(a):
    c,b,A1=rts(a); d=0
    m=(A1-b)*(c-d)/((A1-c)*(b-d)); al2=(A1-b)/(A1-c); g=2/sqrt((A1-c)*(b-d))
    K,E,P=ellipk(m),ellipe(m),ellippi(al2,m)
    V2=(al2*E+(m-al2)*K+(2*al2*m+2*al2-al2**2-3*m)*P)/(2*(al2-1)*(m-al2))
    I0=g*K; I1=g*(c*K+(b-c)*P); I2=g*(c**2*K+2*c*(b-c)*P+(b-c)**2*V2)
    n2=al2*(c-1)/(b-1)
    J=g*(K/(c-1)+(c-b)/((b-1)*(c-1))*ellippi(n2,m))
    return 6*I0-21*I1+15*I2+6*(1-a**2)*J
def A_num(a):
    c,b,A1=rts(a)
    W=lambda r: r*(4*a**2-r*(r-3)**2)
    return quad(lambda r:(6-21*r+15*r**2+6*(1-a**2)/(r-1))/sqrt(W(r)),[b,A1])
if __name__=="__main__":
    for a in ['0.001','0.3','0.7','0.9','0.998','0.9999999']:
        a=mpf(a); x=A_closed(a); print(a, x, A_num(a), x/(27*pi))
    print(16*pi+15*sqrt(3))
