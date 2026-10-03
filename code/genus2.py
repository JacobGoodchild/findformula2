# Test whether sextic N(r) (Kerr shadow edge, general a, theta) has an involution splitting its Jacobian
import itertools, numpy as np, sympy as sp
r,a,c=sp.symbols('r a c',positive=True)
s2=1-c**2
N=sp.expand(s2*r**3*(4*a**2-r*(r-3)**2)+a**4*c**2*s2*(r-1)**2-c**2*(r**2*(3-r)-a**2*(r+1))**2)
print(sp.factor(N))
def test(av,cv):
    coeffs=[float(sp.N(x)) for x in sp.Poly(N.subs({a:av,c:cv}),r).all_coeffs()]
    rts=np.roots(coeffs)
    best=1e9
    for perm in itertools.permutations(range(6)):
        if perm[0]>perm[1] or perm[2]>perm[3] or perm[4]>perm[5] or perm[0]>perm[2] or perm[2]>perm[4]: continue
        M=[]
        for i in range(3):
            x,y=rts[perm[2*i]],rts[perm[2*i+1]]
            M.append([1,-(x+y),x*y])
        d=abs(np.linalg.det(np.array(M)))
        best=min(best,d)
    return rts,best
for av,cv in [(0.5,0.5),(0.9,0.3),(0.7,0.8)]:
    rts,b=test(av,cv); print(av,cv,np.round(rts,5),b)
