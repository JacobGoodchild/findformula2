import numpy as np
def A_formula(a):
    if a<=0.5: return 16*np.pi+8*np.pi*a*a
    return 4*np.pi*(a*a+2)+8*(a*a+2)*np.arcsin(1/(2*a))+(14*a*a+1)*np.sqrt(4*a*a-1)/(a*a)
def A_shoelace(a,N=2000000):
    # general KN formulas (not the simplified extremal ones), e^2=1-a^2
    e2=1-a*a
    lo=max(1.0+1e-5,2-2*a); hi=2+2*a
    t=np.linspace(0,np.pi,N); r=lo+(hi-lo)*(1-np.cos(t))/2
    xi=-(r**3-3*r**2+a*a*r+a*a+2*e2*r)/(a*(r-1))
    eta=r**2*((r*r-3*r+2*e2)**2-4*a*a*(r-e2))*(-1)/(a*a*(r-1)**2)
    al=-xi; be=np.sqrt(np.clip(eta,0,None))
    X=np.concatenate([al,al[::-1]]); Y=np.concatenate([be,-be[::-1]])
    return 0.5*abs(np.dot(X,np.roll(Y,-1))-np.dot(Y,np.roll(X,-1)))
for a in [0.1,0.3,0.5,0.6,0.8,0.95,1.0]:
    print(a, A_formula(a), A_shoelace(a))
