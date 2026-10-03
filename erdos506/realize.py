import itertools, sys, time, numpy as np
from scipy.optimize import least_squares
from ortools.sat.python import cp_model
from extend import mobius_score, count
import designs as D
rng=np.random.default_rng(int(sys.argv[1])); secs=float(sys.argv[2]); TARGET=int(sys.argv[3]) if len(sys.argv)>3 else 24
n=9; blocks=D.blocks; triples=D.triples
def build():
    m,x,y=D.model()
    m.Add(sum(y)<=10)
    pairs=list(itertools.combinations(range(n),2)); ordv=[]
    for p in pairs:
        ov=m.NewBoolVar(""); ordv.append(ov)
        on=[y[i] for i,b in enumerate(blocks) if set(p)<=b]
        m.Add(sum(on)==0).OnlyEnforceIf(ov); m.Add(sum(on)>=1).OnlyEnforceIf(ov.Not())
    m.Add(sum(ordv)>=5)
    w=rng.integers(-5,6,size=len(blocks))
    m.Minimize(sum(int(w[i])*x[i] for i in range(len(blocks))))
    return m,x,y
def residuals(z, circ, lines):
    P=z.reshape(n,2); P=np.vstack([[0,0],[1,0],P[2:]]) if False else P
    r=[]
    for b in lines:
        a=P[b[0]]
        for c in b[2:]:
            u=P[b[1]]-a; v=P[c]-a; r.append(u[0]*v[1]-u[1]*v[0])
    for b in circ:
        M0=[[P[i][0]**2+P[i][1]**2,P[i][0],P[i][1],1] for i in b[:3]]
        for c in b[3:]:
            M=np.array(M0+[[P[c][0]**2+P[c][1]**2,P[c][0],P[c][1],1]]); r.append(np.linalg.det(M))
    # gauge: fix point0=(0,0), point1=(1,0)
    r+= [P[0][0],P[0][1],P[1][0]-1,P[1][1]]
    return np.array(r)
s=cp_model.CpSolver(); s.parameters.num_workers=1; s.parameters.max_time_in_seconds=20
t0=time.time(); tried=0; best=99
while time.time()-t0<secs:
    m,x,y=build(); s.parameters.random_seed=int(rng.integers(1<<30))
    r=s.Solve(m)
    if r not in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        print('no design this round',s.StatusName(r),flush=True); continue
    xb=[i for i in range(len(blocks)) if s.Value(x[i])]; yb=[i for i in xb if s.Value(y[i])]
    lines=[sorted(blocks[i]) for i in yb]; circ=[sorted(blocks[i]) for i in xb if i not in yb and len(blocks[i])>=4]
    tried+=1; ok=False
    for start in range(12):
        z0=rng.normal(size=2*n)*2
        sol=least_squares(residuals,z0,args=(circ,lines),method='trf',max_nfev=400,xtol=1e-14,ftol=1e-14,gtol=1e-14)
        if np.max(np.abs(sol.fun))<1e-10:
            P=[tuple(p) for p in sol.x.reshape(n,2)]
            if min(np.hypot(a[0]-b[0],a[1]-b[1]) for a,b in itertools.combinations(P,2))<1e-4: continue
            c,B,L,mb,O=mobius_score(P)
            if mb>5 or B<=1: continue
            if c<best: best=c; print(f"{time.time()-t0:.0f}s design#{tried}: REALIZED circles={c} B={B} L={L} maxblock={mb} pts={P}",flush=True)
            ok=True; break
    if tried%20==0: print(f"{time.time()-t0:.0f}s designs tried {tried}, best realized {best}",flush=True)
print("FINAL designs",tried,"best realized",best)
