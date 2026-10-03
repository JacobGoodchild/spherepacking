import itertools, numpy as np, sys, time
from math import sqrt
MAXB=int(sys.argv[4]) if len(sys.argv)>4 else 5
rng=np.random.default_rng(int(sys.argv[1])); secs=float(sys.argv[2]); n=int(sys.argv[3]) if len(sys.argv)>3 else 9
phi=(1+sqrt(5))/2
pool=[]
def add(v):
    v=np.array(v,float); v/=np.linalg.norm(v)
    if not any(np.allclose(v,w,atol=1e-9) for w in pool): pool.append(v)
def cyc(v):
    a,b,c=v; return [(a,b,c),(b,c,a),(c,a,b)]
def signs(v):
    out=set()
    for s in itertools.product((1,-1),repeat=3): out.add(tuple(x*y for x,y in zip(s,v)))
    return out
base=[(1,1,1),(1,0,0),(1,1,0),(0,1,phi),(0,1/phi,phi),(phi,0,0),(1/2,phi/2,phi*phi/2)]
for b in base:
    for c in cyc(b):
        for v in signs(c):
            if any(v): add(v)
P=np.array(pool); N=len(P); print("pool",N,flush=True)
def planes_of(S):
    X=P[S]; keys={}
    for i,j,k in itertools.combinations(range(len(S)),3):
        nr=np.cross(X[j]-X[i],X[k]-X[i]); nn=np.linalg.norm(nr)
        nr/=nn; d=nr@X[i]
        idx=np.argmax(np.abs(nr)>1e-9)
        if nr[idx]<0: nr,d=-nr,-d
        keys.setdefault(tuple(np.round(np.append(nr,d),7)), np.append(nr,d))
    return np.array(list(keys.values()))
def score(S):
    Q=planes_of(S); nrm=Q[:,:3]; d=Q[:,3]; X=P[S]
    if len(Q)<=1: return 999,len(Q),0,None
    if (np.abs(X@nrm.T-d[None,:])<1e-7).sum(0).max()>MAXB: return 500,len(Q),0,None
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
        C=np.array(cand); C=C[~(np.abs(C@X.T-1)<1e-7).any(1)]
        if len(C):
            Ls=(np.abs(C@nrm.T-d[None,:])<1e-7).sum(1); i=Ls.argmax(); L=int(Ls[i]); pole=C[i]
    return len(Q)-L, len(Q), L, pole
best=99; t0=time.time(); found=[]
while time.time()-t0<secs:
    S=list(rng.choice(N,n,replace=False)); cur=score(S)[0]; T=2.0
    for step in range(3000):
        S2=S.copy(); out=rng.integers(n); 
        new=rng.integers(N)
        if new in S2: continue
        S2[out]=new; sc=score(S2)[0]
        if sc<=cur or rng.random()<np.exp((cur-sc)/T): S,cur=S2,sc
        T=max(0.05,T*0.998)
        if cur<best:
            best=cur; r=score(S); print(f"{time.time()-t0:.0f}s best {best} B={r[1]} L={r[2]}",flush=True)
            found=[(sorted(S),r[3])]
        elif cur==best and len(found)<50 and sorted(S) not in [f[0] for f in found]:
            found.append((sorted(S),score(S)[3]))
print("FINAL best",best,"distinct optima",len(found))
np.save(f"sa_best_{sys.argv[1]}.npy",np.array([np.concatenate([P[s].ravel(),p if p is not None else np.zeros(3)]) for s,p in found]))
