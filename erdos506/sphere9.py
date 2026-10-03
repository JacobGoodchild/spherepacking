# Search 9-subsets of a symmetric pool of points on the unit sphere.
# B = number of distinct planes spanned by triples (Mobius invariant = #circles+#lines after projection).
# Projecting from pole p (on sphere, not a chosen point) turns planes through p into lines.
# #circles = B - L(p). Float search; winners re-verified exactly later.
import itertools, numpy as np, sys, time
from math import sqrt
s3=sqrt(3)
pool=[]
def add(v):
    v=np.array(v,float); v/=np.linalg.norm(v)
    if not any(np.allclose(v,w) for w in pool): pool.append(v)
for sx in (1,-1):
    for sy in (1,-1):
        for sz in (1,-1): add((sx,sy,sz))                      # cube
for i in range(3):
    for s in (1,-1): e=[0,0,0]; e[i]=s; add(e)                  # octahedron
for i,j in itertools.combinations(range(3),2):
    for a in (1,-1):
        for b in (1,-1): e=[0,0,0]; e[i]=a; e[j]=b; add(e)     # cuboctahedron
P=np.array(pool); N=len(P); print("pool",N)
def plane(a,b,c):
    nrm=np.cross(b-a,c-a); 
    if np.linalg.norm(nrm)<1e-12: return None
    nrm/=np.linalg.norm(nrm); d=nrm@a
    if nrm[np.argmax(np.abs(nrm)>1e-9)]<0: nrm,d=-nrm,-d
    return tuple(np.round(np.append(nrm,d),8))
ids={}; T=np.zeros((N,N,N),np.int32)
for i,j,k in itertools.combinations(range(N),3):
    v=ids.setdefault(plane(P[i],P[j],P[k]),len(ids))
    for a,b,c in itertools.permutations((i,j,k)): T[a,b,c]=v
planes=[None]*len(ids)
for key,v in ids.items(): planes[v]=np.array(key)
PL=np.array(planes)                       # (nplanes, 4): n.x = d
def bestpole(bl, pts):
    Q=PL[bl]; n=Q[:,:3]; d=Q[:,3]; cand=[]
    for a,b in itertools.combinations(range(len(Q)),2):
        u=np.cross(n[a],n[b]); nu=np.linalg.norm(u)
        if nu<1e-9: continue
        x0=np.linalg.solve(np.array([n[a],n[b],u]),[d[a],d[b],0]); uu=u/nu
        bq=x0@uu; disc=bq*bq-(x0@x0-1)
        if disc<-1e-12: continue
        for sg in (1,-1): cand.append(x0+(-bq+sg*sqrt(max(disc,0)))*uu)
    if not cand: return 1,None
    C=np.array(cand)
    C=C[~(np.abs(C@pts.T-1)<1e-7).any(1)]
    if len(C)==0: return 1,None
    L=(np.abs(C@n.T-d[None,:])<1e-7).sum(1)
    i=int(L.argmax()); return int(L[i]),C[i]
tri=np.array(list(itertools.combinations(range(9),3)))
best=int(sys.argv[1]) if len(sys.argv)>1 else 26; t0=time.time(); res=[]
it=itertools.combinations(range(N),9)
while True:
    ch=np.fromiter(itertools.chain.from_iterable(itertools.islice(it,100000)),np.int16)
    if ch.size==0: break
    S=ch.reshape(-1,9).astype(np.int64)
    V=T[S[:,tri[:,0]],S[:,tri[:,1]],S[:,tri[:,2]]]; V.sort(1)
    B=(np.diff(V,1)!=0).sum(1)+1
    for r in np.nonzero(B<=best+12)[0]:
        bl=np.unique(V[r]); 
        if len(bl)<=1: continue
        L,pole=bestpole(bl,P[S[r]]); c=len(bl)-L
        if c<best: best=c; res=[]; print("new best",c,"B",len(bl),"L",L,flush=True)
        if c==best and len(res)<20: res.append((S[r].copy(),pole))
print("best circles",best,"time",time.time()-t0)
np.save("sphere_best.npy",np.array([np.concatenate([P[s].ravel(),p]) for s,p in res]))
