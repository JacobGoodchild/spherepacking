"""Symmetric families + forced incidences. Pick d incidence equations (concyclic 4-sets
or collinear 3-sets), solve for the d family parameters, score with best inversion."""
import itertools, sys, time, numpy as np
from math import cos, sin, pi
from scipy.optimize import least_squares
from extend import mobius_score
rng=np.random.default_rng(int(sys.argv[1])); secs=float(sys.argv[2]); fam=sys.argv[3]
w=np.exp(2j*pi/3)
def pts(p):
    if fam=='C3':      # orbits: 1*w^k, r2 e^{i t2} w^k, r3 e^{i t3} w^k
        r2,t2,r3,t3=p; z=[1,r2*np.exp(1j*t2),r3*np.exp(1j*t3)]
        return [zz*w**k for zz in z for k in range(3)]
    if fam=='C3c':     # 2 triangle orbits + 3 points? -> orbits 3+3+3 with one orbit through pole-free; plus generic translation of centre not needed
        r2,t2,r3,t3=p; z=[1,r2*np.exp(1j*t2),r3*np.exp(1j*t3)]
        return [zz*w**k for zz in z for k in range(3)]
    if fam=='C4':      # two squares + one free point on axis: 4+4+1
        r2,t2,a,b=p; i4=1j
        return [i4**k for k in range(4)]+[r2*np.exp(1j*t2)*i4**k for k in range(4)]+[a+1j*b]
    if fam=='C2':      # 4 pairs symmetric under z->-z plus the centre... 9 = 4*2 + 1 (centre)
        z=[p[0]+1j*p[1],p[2]+1j*p[3],p[4]+1j*p[5],1.0]
        return [s*zz for zz in z for s in (1,-1)]+[0j]
    if fam=='D1':      # mirror symmetry: 4 pairs (x,y),(-x,y) + 1 on axis
        z=[(p[0],p[1]),(p[2],p[3]),(p[4],p[5]),(1.0,0.0)]
        return [complex(s*x,y) for x,y in z for s in (1,-1)]+[complex(0,p[6])]
D={'C3':4,'C3c':4,'C4':4,'C2':6,'D1':7}[fam]
def det4(P,q):
    M=np.array([[abs(z)**2,z.real,z.imag,1] for z in (P[i] for i in q)]); return np.linalg.det(M)
def col3(P,q):
    a,b,c=(P[i] for i in q); return ((b-a).conjugate()*(c-a)).imag
quads=list(itertools.combinations(range(9),4)); trips=list(itertools.combinations(range(9),3))
def res(p,eqs):
    P=pts(p); return np.array([det4(P,q) if len(q)==4 else col3(P,q) for q in eqs])
best=99; t0=time.time(); tries=0
while time.time()-t0<secs:
    k=D
    eqs=[quads[i] for i in rng.choice(len(quads),k,replace=False)]
    if rng.random()<0.4: eqs[-1]=trips[rng.integers(len(trips))]
    p0=rng.uniform(-2.5,2.5,size=D)
    try:
        sol=least_squares(res,p0,args=(eqs,),method='lm',max_nfev=200*D,xtol=1e-15,ftol=1e-15)
    except Exception: continue
    tries+=1
    if np.max(np.abs(sol.fun))>1e-11: continue
    P=pts(sol.x)
    if min(abs(a-b) for a,b in itertools.combinations(P,2))<1e-3 or max(abs(z) for z in P)>1e3: continue
    c,B,L,mb,O=mobius_score([(z.real,z.imag) for z in P])
    if B<=1 or mb>5: continue
    if c<best or c<=24:
        best=min(best,c); print(f"{time.time()-t0:.0f}s try{tries} circles={c} B={B} L={L} maxblock={mb} params={list(sol.x)}",flush=True)
print("FINAL",fam,"best",best,"tries",tries)
