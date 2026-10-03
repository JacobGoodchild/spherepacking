import sys, itertools, time, numpy as np
sys.path.insert(0,'/home/user/spherepacking/erdos506')
from circles import circle_key
from extend import mobius_score
G=int(sys.argv[1]); BMAX=int(sys.argv[2]); kind=sys.argv[3] if len(sys.argv)>3 else 'sq'
if kind=='sq': P=[(x,y) for x in range(G) for y in range(G)]
else:  # triangular lattice in integer coords (a,b) -> (2a+b, b*sqrt3) ; use exact keys via line/circle in Q(sqrt3): fall back to float ids
    P=[(a,b) for a in range(G) for b in range(G)]
N=len(P)
def fl(p): return (p[0]+p[1]/2, p[1]*np.sqrt(3)/2) if kind=='tri' else (float(p[0]),float(p[1]))
from extend import blocks
ids={}; T=np.zeros((N,N,N),np.int32)
for i,j,k in itertools.combinations(range(N),3):
    if kind=='sq':
        key=circle_key(P[i],P[j],P[k])
        if key is None:  # line through i,j: canonical by the set of collinear grid points
            (x1,y1),(x2,y2)=P[i],P[j]; a,b=y2-y1,x1-x2; c=-(a*x1+b*y1)
            from math import gcd; g=gcd(gcd(abs(a),abs(b)),abs(c)); a,b,c=a//g,b//g,c//g
            if a<0 or (a==0 and b<0): a,b,c=-a,-b,-c
            key=('L',a,b,c)
    else:
        bl=blocks([fl(P[i]),fl(P[j]),fl(P[k])]); key=list(bl.keys())[0]
    v=ids.setdefault(key,len(ids))
    for a,b,c in itertools.permutations((i,j,k)): T[a,b,c]=v
tri=np.array(list(itertools.combinations(range(9),3)))
best=99; t0=time.time(); seen=0
it=itertools.combinations(range(N),9)
while True:
    ch=np.fromiter(itertools.chain.from_iterable(itertools.islice(it,200000)),np.int16)
    if ch.size==0: break
    S=ch.reshape(-1,9).astype(np.int64)
    V=T[S[:,tri[:,0]],S[:,tri[:,1]],S[:,tri[:,2]]]; V.sort(1)
    B=(np.diff(V,1)!=0).sum(1)+1
    for r in np.nonzero((B<=BMAX)&(B>1))[0]:
        seen+=1
        pts=[fl(P[t]) for t in S[r]]
        c,Bb,L,mb,O=mobius_score(pts)
        if mb>5: continue
        if c<best or c<=24:
            best=min(best,c); print(f"{time.time()-t0:.0f}s circles={c} B={Bb} L={L} maxblock={mb} pts={[P[t] for t in S[r]]} O={O}",flush=True)
print(f"FINAL {kind} G={G}: best (maxblock<=5) = {best}, candidates scored {seen}, {time.time()-t0:.0f}s")
