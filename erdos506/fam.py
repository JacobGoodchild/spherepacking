import itertools, numpy as np, sys, time
from math import sqrt, cos, sin, pi
TOL=1e-7
def score(X):
    n=len(X); keys={}
    for i,j,k in itertools.combinations(range(n),3):
        nr=np.cross(X[j]-X[i],X[k]-X[i]); nn=np.linalg.norm(nr)
        if nn<1e-9: return None
        nr/=nn; d=nr@X[i]; idx=np.argmax(np.abs(nr)>1e-9)
        if nr[idx]<0: nr,d=-nr,-d
        keys.setdefault(tuple(np.round(np.append(nr,d),6)), np.append(nr,d))
    Q=np.array(list(keys.values())); nrm=Q[:,:3]; d=Q[:,3]
    if len(Q)<=1: return None
    maxb=(np.abs(X@nrm.T-d[None,:])<TOL).sum(0).max()
    cand=[]
    for a,b in itertools.combinations(range(len(Q)),2):
        u=np.cross(nrm[a],nrm[b]); nu=np.linalg.norm(u)
        if nu<1e-9: continue
        x0=np.linalg.solve(np.array([nrm[a],nrm[b],u]),[d[a],d[b],0]); uu=u/nu
        bq=x0@uu; disc=bq*bq-(x0@x0-1)
        if disc<-1e-12: continue
        r=sqrt(max(disc,0)); cand.append(x0+(-bq+r)*uu); cand.append(x0+(-bq-r)*uu)
    L=1; pole=None
    if cand:
        C=np.array(cand); C=C[~(np.abs(C@X.T-1)<TOL).any(1)]
        if len(C):
            Ls=(np.abs(C@nrm.T-d[None,:])<TOL).sum(1); i=Ls.argmax(); L=int(Ls[i]); pole=C[i]
    return len(Q)-L, len(Q), L, int(maxb), pole
def ring(z,k,off):
    r=sqrt(max(0,1-z*z)); return [np.array([r*cos(off+2*pi*j/k), r*sin(off+2*pi*j/k), z]) for j in range(k)]
if __name__=="__main__":
    Z=sorted(set([0,1/3,1/2,1/sqrt(2),sqrt(3)/2,1/sqrt(3),0.6,0.8,1/sqrt(5),2/sqrt(5),sqrt(2/3),1/4,3/4,sqrt(5)/3,2/3,1/sqrt(10),3/sqrt(10)]))
    Z=sorted(set(Z+[-z for z in Z]))
    comps=[(3,3,3),(4,4,1),(4,1,4),(3,3,2,1),(4,4),(4,5),(3,6),(2,2,2,2,1),(3,3,1,1,1),(4,2,2,1),(2,4,2,1),(1,4,4)]
    offs=[0,pi/12,pi/6,pi/4,pi/3,pi/2,2*pi/3, pi]
    rng=np.random.default_rng(int(sys.argv[1])); secs=float(sys.argv[2])
    best=99; best5=99; t0=time.time(); cnt=0
    while time.time()-t0<secs:
        comp=comps[rng.integers(len(comps))]
        zs=rng.choice(Z,len(comp),replace=False)
        X=[]
        for k,z in zip(comp,zs):
            X+=ring(z,k,offs[rng.integers(len(offs))] if k>1 else rng.choice(offs)+rng.choice(offs))
        X=np.array(X)
        if len(np.unique(np.round(X,6),axis=0))<9: continue
        s=score(X); cnt+=1
        if s is None: continue
        c,B,L,mb,pole=s
        if c<best or (mb<=5 and c<best5):
            if c<best: best=c
            if mb<=5 and c<best5: best5=c
            print(f"{time.time()-t0:.0f}s circles={c} B={B} L={L} maxblock={mb} comp={comp} | best={best} best(maxblock<=5)={best5}",flush=True)
            np.save(f"fam_{sys.argv[1]}_{c}_{mb}.npy",np.vstack([X,pole if pole is not None else np.zeros(3)]))
    print("FINAL",best,best5,"evaluated",cnt)
