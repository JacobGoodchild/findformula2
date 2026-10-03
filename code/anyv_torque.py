from mpmath import mp, mpf, pi, sqrt, cos, sin, quad, findroot, diff, nstr
import anyv_area as AA
mp.dps=30
def mean_Lz_per_E(Ev,av,th):
    Ev=mpf(Ev); av=mpf(av); p2=Ev*Ev-1; c2=cos(th)**2; s2=sin(th)**2
    g=lambda rr: AA.Qf(rr,av,Ev)+av*av*p2*c2-AA.Lf(rr,av,Ev)**2*c2/s2
    E2=Ev*Ev; disc=(3*E2-4)**2-16*(1-E2)
    rc=[x for x in [(-(3*E2-4)+sg*sqrt(disc))/(2*(1-E2)) for sg in (1,-1)] if 3<x.real<=4][0]
    w=6*av; rs=[rc-w+2*w*k/3000 for k in range(3001)]
    idx=[i for i,x in enumerate(rs) if g(x).real>0]
    lo=findroot(g,(rs[idx[0]-1],rs[idx[0]]),solver='bisect'); hi=findroot(g,(rs[idx[-1]],rs[idx[-1]+1]),solver='bisect')
    Lp=lambda x: diff(lambda z: AA.Lf(z,av,Ev),x)
    def integ(t,wt):
        x=lo+(hi-lo)*(1-cos(t))/2
        return wt(x)*sqrt(max(g(x).real,0))*abs(Lp(x))*(hi-lo)*sin(t)/2
    area=quad(lambda t: integ(t,lambda x:1),[0,pi/2,pi])
    mom=quad(lambda t: integ(t,lambda x:AA.Lf(x,av,Ev)),[0,pi/2,pi])
    return mom/area/Ev, rc
for E,th in [('1.25',pi/2),('2',pi/3),('1.0001',pi/2)]:
    v,rc=mean_Lz_per_E(E,'0.001',th)
    print(E, nstr(th,4), nstr(v/mpf('0.001'),12), 'pred', nstr(-2*sin(th)**2/(rc-2),12))
