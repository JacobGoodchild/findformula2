from mpmath import mp, mpf, pi, sqrt, nstr, matrix, lu_solve
import anyv_area as AA
mp.dps=40
def fcoef(E,th):
    hs=[mpf('0.01'),mpf('0.02'),mpf('0.03')]
    ys=[AA.area(E,h,th) for h in hs]
    M=matrix([[1,h**2,h**4] for h in hs]); c=lu_solve(M,matrix(ys))
    return c[0], -c[1]/c[0]
if __name__=="__main__":
    for E in ['1.0001','1.25','2','10']:
        A0e,fe=fcoef(E,pi/2); A0p,fp=fcoef(E,pi/180*1)
        print(E, nstr(A0e,12), 'f_edge',nstr(fe,12),'f_pole',nstr(fp,12), flush=True)
