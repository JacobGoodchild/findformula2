import numpy as np
from isco_spectrum import spectrum_sim
gs,tot=spectrum_sim(v=0.0)
bins=np.linspace(0,1.1,12); h,_=np.histogram(gs,bins)
def pred(g):
    if 0.5<=g<=1: return 0.5
    if 0<g<0.5: return np.arcsin(min(np.sqrt(3)*g/np.sqrt(1-g*g),1))/(2*np.pi)
    return 0
for k in range(len(h)):
    c=(bins[k]+bins[k+1])/2; print(round(c,3), round(h[k]/tot/(bins[1]-bins[0]),4), round(pred(c),4))
print('total', len(gs)/tot, 7/24, ' mean g', gs.sum()/tot, 3/32+np.sqrt(3)/16)
