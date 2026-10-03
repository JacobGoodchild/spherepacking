import itertools as it, sys, time
import numpy as np
from ortools.sat.python import cp_model
n=9; SMAX=5; TARGET=24
blocks=[frozenset(c) for s in range(3,SMAX+1) for c in it.combinations(range(n),s)]
bidx={b:i for i,b in enumerate(blocks)}
m=cp_model.CpModel()
x=[m.NewBoolVar("") for _ in blocks]; y=[m.NewBoolVar("") for _ in blocks]
for t in it.combinations(range(n),3):
    m.AddExactlyOne([x[i] for i,b in enumerate(blocks) if set(t)<=b])
for i in range(len(blocks)): m.AddImplication(y[i],x[i])
for i,j in it.combinations(range(len(blocks)),2):
    if len(blocks[i]&blocks[j])>=2: m.AddBoolOr([y[i].Not(),y[j].Not()])
ordv=[]
for p in it.combinations(range(n),2):
    ov=m.NewBoolVar(""); ordv.append(ov)
    on=[y[i] for i,b in enumerate(blocks) if set(p)<=b]
    m.Add(sum(on)==0).OnlyEnforceIf(ov); m.Add(sum(on)>=1).OnlyEnforceIf(ov.Not())
m.Add(sum(ordv)>=5)
m.Add(sum(ordv) >= 3 + sum((len(b)-3)*y[i] for i,b in enumerate(blocks) if len(b)>=4))
m.Add(sum(ordv)+sum(y) >= n)
for q in range(n):
    thr=[(i,b) for i,b in enumerate(blocks) if q in b]
    m.Add(sum(x[i] for i,b in thr if len(b)==3)>=4)
    m.Add(sum(x[i] for i,b in thr if len(b)==3) >= 3 + sum((len(b)-4)*x[i] for i,b in thr if len(b)>=5))
    m.Add(sum(x[i] for i,b in thr) >= n-1)
deg=[sum(x[i] for i,b in enumerate(blocks) if q in b and len(b)>=4) for q in range(n)]
for q in range(n-1): m.Add(deg[q]>=deg[q+1])
m.Add(sum(x)-sum(y)<=TARGET)
def labelled(conf,k):
    out=set()
    for perm in it.permutations(range(k)):
        out.add(frozenset(frozenset(perm[i] for i in L) for L in conf))
    return [list(map(sorted,c)) for c in out]
FANO=labelled([(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)],7)
MK=labelled([(0,1,3),(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,0),(6,7,1),(7,0,2)],8)
PAP=labelled([(0,1,2),(3,4,5),(6,7,8),(0,4,8),(0,5,7),(1,3,8),(1,5,6),(2,3,7),(2,4,6)],9)
CUBE=labelled([(0,1,2,3),(4,5,6,7),(0,1,4,5),(2,3,6,7),(0,2,4,6),(1,3,5,7)],8)
def zx(q,trip): return [x[i] for i,b in enumerate(blocks) if q in b and set(trip)<=b]
def zy(trip): return [y[i] for i,b in enumerate(blocks) if set(trip)<=b]
def zq(Q): return [x[i] for i,b in enumerate(blocks) if set(Q)<=b]
def check(X,Y):
    """Return list of violated-constraint adders for this design."""
    cover=lambda S,blk: any(set(S)<=b for b in blk)
    viol=[]
    for q in range(n):
        thr=[b for b in X if q in b]; oth=[p for p in range(n) if p!=q]
        for S in it.combinations(oth,7):
            for F in FANO:
                if all(cover([S[i] for i in L]+[q],thr) for L in F):
                    viol.append(('fq',q,[[S[i] for i in L] for L in F]))
        for F in MK:
            if all(cover([oth[i] for i in L]+[q],thr) for L in F):
                viol.append(('mq',q,[[oth[i] for i in L] for L in F]))
    for S in it.combinations(range(n),7):
        for F in FANO:
            if all(cover([S[i] for i in L],Y) for L in F): viol.append(('fy',None,[[S[i] for i in L] for L in F]))
    for S in it.combinations(range(n),8):
        for F in MK:
            if all(cover([S[i] for i in L],Y) for L in F): viol.append(('my',None,[[S[i] for i in L] for L in F]))
        for C in CUBE:
            fs=[cover([S[i] for i in f],X) for f in C]
            if sum(fs)==5: viol.append(('mi',None,[[S[i] for i in f] for f in C]))
    for F in PAP:
        ls=[cover(L,Y) for L in F]
        if sum(ls)==8: viol.append(('pa',None,F))
    return viol
def add(v):
    kind,q,F=v
    if kind=='fq': m.Add(sum(sum(zx(q,L)) for L in F)<=6)
    elif kind=='mq': m.Add(sum(sum(zx(q,L)) for L in F)<=7)
    elif kind=='fy': m.Add(sum(sum(zy(L)) for L in F)<=6)
    elif kind=='my': m.Add(sum(sum(zy(L)) for L in F)<=7)
    elif kind=='mi':
        fs=[sum(zq(f)) for f in F]
        for j in range(6): m.Add(sum(fs[k] for k in range(6) if k!=j)-fs[j]<=4)
    elif kind=='pa':
        zs=[sum(zy(L)) for L in F]
        for j in range(9): m.Add(sum(zs[k] for k in range(9) if k!=j)-zs[j]<=7)
s=cp_model.CpSolver(); s.parameters.num_workers=1; s.parameters.max_time_in_seconds=600
t0=time.time(); it_=0
while True:
    it_+=1; r=s.Solve(m)
    if r==cp_model.INFEASIBLE: print(f"INFEASIBLE after {it_} rounds, {time.time()-t0:.0f}s",flush=True); break
    if r not in (cp_model.OPTIMAL,cp_model.FEASIBLE): print("unknown",s.StatusName(r),flush=True); break
    X=[blocks[i] for i in range(len(blocks)) if s.Value(x[i])]; Y=[blocks[i] for i in range(len(blocks)) if s.Value(y[i])]
    V=check(X,Y)
    if not V:
        print(f"round {it_}: CLEAN design B={len(X)} L={len(Y)}",flush=True)
        print(" blocks",[sorted(b) for b in X]); print(" lines",[sorted(b) for b in Y],flush=True)
        # block it and keep going to collect more candidates
        m.AddBoolOr([x[bidx[b]].Not() for b in X]+[y[bidx[b]].Not() for b in Y])
        continue
    for v in V[:200]: add(v)
    if it_%10==0: print(f"round {it_}: {len(V)} violations, {time.time()-t0:.0f}s",flush=True)
