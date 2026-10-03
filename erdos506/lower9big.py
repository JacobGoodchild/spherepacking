"""Abstract lower bound for Erdos #506 at n=9 (sphere model).
9 points on a sphere, not all on one circle; blocks = planes through >=3 points
(an exact cover of the 84 triples). Pole = projection centre (not a point).
#circles = B - L, L = blocks through the pole.
Necessary conditions used:
 C1 blocks partition the triples (3 points determine a plane).
 C2 blocks through the pole pairwise share <= 1 point.
 C3 projecting from the pole, the 9 points are not collinear, so by Csima-Sawyer
    there are >= ceil(6*9/13) = 5 ordinary lines: pairs lying in no pole-block.
 C4 inverting at any point q, the other 8 points are not collinear, so they have
    >= ceil(6*8/13) = 4 ordinary lines: blocks of size exactly 3 containing q.
"""
import itertools, sys, time
from ortools.sat.python import cp_model
n=9; TARGET=int(sys.argv[1]); SMAX=int(sys.argv[2]); C3=int(sys.argv[3]); C4=int(sys.argv[4]); T=float(sys.argv[5])
blocks=[frozenset(c) for s in range(3,SMAX+1) for c in itertools.combinations(range(n),s)]
m=cp_model.CpModel()
x=[m.NewBoolVar("") for _ in blocks]; y=[m.NewBoolVar("") for _ in blocks]
for t in itertools.combinations(range(n),3):
    m.AddExactlyOne([x[i] for i,b in enumerate(blocks) if set(t)<=b])
for i in range(len(blocks)): m.AddImplication(y[i],x[i])
for i,j in itertools.combinations(range(len(blocks)),2):
    if len(blocks[i]&blocks[j])>=2: m.AddBoolOr([y[i].Not(),y[j].Not()])
ordv=[]
for p in itertools.combinations(range(n),2):
    ov=m.NewBoolVar(""); ordv.append(ov)
    on=[y[i] for i,b in enumerate(blocks) if set(p)<=b]
    m.Add(sum(on)==0).OnlyEnforceIf(ov); m.Add(sum(on)>=1).OnlyEnforceIf(ov.Not())
m.Add(sum(ordv)>=C3)
for q in range(n):
    m.Add(sum(x[i] for i,b in enumerate(blocks) if q in b and len(b)==3)>=C4)
# Melchior (real, not all collinear): t2 >= 3 + sum_{k>=4} (k-3) t_k
# at the pole: lines are pole-blocks (size k) and ordinary pairs
m.Add(sum(ordv) >= 3 + sum((len(b)-3)*y[i] for i,b in enumerate(blocks) if len(b)>=4))
# de Bruijn-Erdos at the pole: #lines >= 9
m.Add(sum(ordv)+sum(y) >= n)
for q in range(n):
    thr=[(i,b) for i,b in enumerate(blocks) if q in b]
    # derived structure at q: block of size k through q -> line with k-1 points among the other 8
    t2=sum(x[i] for i,b in thr if len(b)==3)
    m.Add(t2 >= 3 + sum((len(b)-1-3)*x[i] for i,b in thr if len(b)-1>=4))
    m.Add(sum(x[i] for i,b in thr) >= n-1)
m.Add(sum(x)-sum(y)<=TARGET)
m.Add(sum(x[i] for i,b in enumerate(blocks) if len(b)>=6)>=1)
s=cp_model.CpSolver(); s.parameters.num_workers=int(sys.argv[6]) if len(sys.argv)>6 else 1; s.parameters.max_time_in_seconds=T
t0=time.time(); r=s.Solve(m)
print(f"target<= {TARGET}, blocks<= {SMAX}, C3>={C3}, C4>={C4}: {s.StatusName(r)} ({time.time()-t0:.1f}s)")
if r in (cp_model.OPTIMAL,cp_model.FEASIBLE):
    xb=[sorted(blocks[i]) for i in range(len(blocks)) if s.Value(x[i])]; yb=[sorted(blocks[i]) for i in range(len(blocks)) if s.Value(y[i])]
    print(" B",len(xb),"L",len(yb),"sizes",sorted(len(b) for b in xb)); print(" lines",yb); print(" blocks",xb)
