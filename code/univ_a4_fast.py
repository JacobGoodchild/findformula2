# Fast version of univ_a4_shadow.py: truncated power series in eps with coefficients in Q(f0..fK, S)[x].
# Units: photon-sphere radius R = 1 (restore by f_k -> R^k f_k, a -> a/R, area -> R^2 area).
import sys, sympy as sp
from sympy.polys.fields import field
from sympy.polys.rings import ring
from sympy import QQ
N=int(sys.argv[1]) if len(sys.argv)>1 else 4
KT=N+4
names=','.join(['f%d'%k for k in range(KT+2)])+',S'
Kf,*gens=field(names,QQ)
fsym=gens[:-1]; S=gens[-1]
Rx,x=ring('x',Kf)
one=Rx(1); zero=Rx(0)
def smul(A,B):
    C=[zero]*(N+2)
    for i,a in enumerate(A):
        if a==0: continue
        for j,b in enumerate(B):
            if i+j<=N+1 and b!=0: C[i+j]+=a*b
    return C
def sadd(A,B): return [a+b for a,b in zip(A,B)]
def sscale(A,c): return [a*c for a in A]
def const(c): return [Rx(c)]+[zero]*(N+1)
from math import factorial
f1=2*fsym[0]          # photon sphere with R = 1: f1 = 2 f0
fk=[fsym[0],f1]+list(fsym[2:])
def Tser(k0):        # f^(k0)(1 + eps x) as series in eps
    out=[zero]*(N+2)
    for k in range(N+2):
        if k0+k<len(fk): out[k]=Rx(fk[k0+k]/factorial(k))*x**k
    return out
r=[one,x]+[zero]*N
r2=smul(r,r)
f=Tser(0); fp=Tser(1)
eps2=[zero,zero,one]+[zero]*(N-1)
D=sadd(smul(r2,f),eps2)
Dp=sadd(sscale(smul(r,f),2),smul(r2,fp))
G=sadd(smul(r2,sadd(sscale(f,2),sscale(smul(r,fp),-1))),sscale(eps2,4))
assert G[0]==0
Ge=G[1:]+[zero]
d0=Dp[0]; d0c=d0.coeff(1) if False else Kf(d0.LC)  # d0 is constant
w=[zero]+[c*(1/d0c) for c in Dp[1:]]
iDp=const(1); wp=const(1)
for k in range(1,N+2):
    wp=smul(wp,w); iDp=sadd(iDp,sscale(wp,(-1)**k))
iDp=sscale(iDp,1/d0c)
epsser=[zero,one]+[zero]*N
xi=sadd(sscale(smul(r,smul(Ge,iDp)),-1),epsser)
eta=smul(smul(r2,sadd(sscale(D,16),sscale(smul(Ge,Ge),-1))),smul(iDp,iDp))
C=1-S**2
beta2=sadd(sadd(eta,sscale(eps2,C)),sscale(smul(xi,xi),-C/S**2))
alpha=sscale(xi,-1/S)
k=alpha[0].coeff(x)
b02=beta2[0]+alpha[0]**2
b02c=b02[(0,)]
assert b02-Rx(b02c)==0, b02
print('k =',k,' b0^2 =',b02c,flush=True)
# change variable x = u/k : polynomial in u
Ru,uu=ring('u',Kf)
def tou(p):
    out=Ru(0)
    for (m,),c in p.terms(): out+=c*(uu/k)**m if m else Ru(c)
    return out
dalpha=[p.diff(x) for p in alpha]
g=[tou(p)*(1/k) for p in dalpha]
Qs=[tou(p) for p in beta2]
delta=[Qs[0]-(Ru(b02c)-uu**2)]+Qs[1:]
assert delta[0]==0
Z=Ru(0)
def umul(A,B):
    Cc=[Z]*(N+1)
    for i,a in enumerate(A[:N+1]):
        if a==0: continue
        for j,b in enumerate(B[:N+1]):
            if i+j<=N and b!=0: Cc[i+j]+=a*b
    return Cc
import sympy
def _q(v):
    v=sympy.Rational(v); return QQ(int(v.p),int(v.q))
def FP(P,j):
    out=Kf(0)
    for (m,),c in P.terms():
        if m%2: continue
        n=(m+2-2*j)//2
        if n<0: continue
        out+=c*_q((-1)**j*sympy.binomial(sympy.Rational(1,2)-j,n))*(-b02c)**n
    return out   # times (-pi)
tot=[Kf(0)]*(N+1); dj=[Ru(1)]+[Z]*N
gN=g[:N+1]
for j in range(0,N+1):
    term=umul(dj,gN)
    cj=_q(sympy.binomial(sympy.Rational(1,2),j))
    for n in range(N+1):
        if term[n]!=0: tot[n]+=cj*FP(term[n],j)
    dj=umul(dj,delta[:N+1])
    print('j',j,flush=True)
A0=tot[0]
fsy=sp.symbols(names)
for n in range(1,N+1):
    rel=tot[n]/A0
    e=sp.factor(rel.as_expr().subs(sp.Symbol('S'),sp.sqrt(1-sp.Symbol('C'))))
    print('a^%d relative (R=1):'%n, e, flush=True)
print('A0/(-pi)*2 =', 2*A0, '  (expect -b0^2/2*2 => area pi b0^2)')
