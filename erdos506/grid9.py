import sys, itertools, time, numpy as np
sys.path.insert(0,'/home/user/spherepacking/erdos506')
from circles import circle_key, count_circles, degenerate
G=int(sys.argv[1]); n=int(sys.argv[2]); MAXB=int(sys.argv[3]) if len(sys.argv)>3 else 0
P=[(x,y) for x in range(G) for y in range(G)]; N=len(P)
ids={}; T=-np.ones((N,N,N),np.int32)
for i,j,k in itertools.combinations(range(N),3):
    key=circle_key(P[i],P[j],P[k])
    v=-1 if key is None else ids.setdefault(key,len(ids))
    for a,b,c in itertools.permutations((i,j,k)): T[a,b,c]=v
tri=np.array(list(itertools.combinations(range(n),3)))
best=10**9; hist={}; t0=time.time(); sols=[]
it=itertools.combinations(range(N),n)
while True:
    chunk=np.fromiter(itertools.chain.from_iterable(itertools.islice(it,200000)),np.int16)
    if chunk.size==0: break
    S=chunk.reshape(-1,n).astype(np.int64)
    V=T[S[:,tri[:,0]],S[:,tri[:,1]],S[:,tri[:,2]]]
    V.sort(axis=1)
    cnt=(np.diff(V,axis=1)!=0).sum(1)+1-(V[:,0]==-1)
    cnt[cnt<=1]=10**6
    m=cnt.min()
    for r in np.nonzero(cnt<=min(best,m+4))[0][:400]:
        pts=[P[t] for t in S[r]]
        if degenerate(pts): continue
        if MAXB:
            from extend import count as fc
            if max(len(b[2]) for b in fc([tuple(map(float,p)) for p in pts])[1].values())>MAXB: continue
        c=count_circles(pts)
        if c<best: best=c; sols=[]
        if c==best and len(sols)<5: sols.append(pts)
print(f"grid {G}x{G}, n={n}: min circles = {best}  ({time.time()-t0:.0f}s)")
for s in sols[:3]: print(s)
