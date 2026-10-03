import sympy as sp
C,q,rc=sp.symbols('C q r_c')
co=sp.sympify(open('kn_anyv_a2.out').read().split('a^2 relative coefficient:')[1].split('\n')[0])
e2=(rc**2-2*rc+q**2)**2/(rc**2*(rc**2-3*rc+2*q**2)); l2=rc**2*(rc-q**2)/(rc**2-3*rc+2*q**2)
b02=sp.factor(l2/(e2-1)); print('b0^2 =',b02)
tot=sp.factor(sp.simplify(co+C/b02))
print('TOTAL a^2 coefficient:',tot)
print('q=0:',sp.factor(tot.subs(q,0)))
num,den=sp.fraction(tot); print('num collected in C:',sp.collect(sp.factor_terms(sp.expand(num)),C))
for cc in [0,1,sp.Rational(1,3)]: print('C=',cc,sp.factor(tot.subs(C,cc)))
# isotropic average = C->1/3
