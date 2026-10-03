import itertools, numpy as np, sys, time
from extend import count, blocks, candidates
MAXB=int(sys.argv[4]) if len(sys.argv)>4 else 5
rng=np.random.default_rng(int(sys.argv[1])); secs=float(sys.argv[2]); kind=sys.argv[3]
def tri_lattice(R):
    s=np.sqrt(3)/2; return [(a+b/2,b*s) for a in range(-R,R+1) for b in range(-R,R+1) if abs(a+b/2)<=R*0.75 and abs(b*s)<=R*0.75]
if kind=='tri': pool=tri_lattice(3)
elif kind=='sq': pool=[(x,y) for x in range(-2,3) for y in range(-3,3)]
elif kind=='closure':
    seed=[(0,0),(1,0),(0,1),(1,1),(2,1)]
    _,B=count(seed); pool=seed+candidates(seed,B)
    pool=[p for p in pool if abs(p[0])<6 and abs(p[1])<6]
pool=[tuple(map(float,p)) for p in pool]; N=len(pool); print(kind,"pool",N,flush=True)
def sc(S):
    P=[pool[i] for i in S]; c,B=count(P)
    if max(len(b[2]) for b in B.values())>MAXB: return 999
    return c
best=99; t0=time.time()
while time.time()-t0<secs:
    S=list(rng.choice(N,9,replace=False)); cur=sc(S); T=1.5
    for step in range(2500):
        S2=S.copy(); new=int(rng.integers(N))
        if new in S2: continue
        S2[rng.integers(9)]=new; s2=sc(S2)
        if s2<=cur or rng.random()<np.exp((cur-s2)/T): S,cur=S2,s2
        T=max(0.05,T*0.998)
        if cur<best:
            best=cur; print(f"{time.time()-t0:.0f}s best {best}: {[pool[i] for i in sorted(S)]}",flush=True)
print("FINAL",best)
