import itertools, sys, time
from math import sqrt
from ortools.sat.python import cp_model
from extend import blocks as fblocks
u=sqrt(15)
P8=[(u,-1),(-u,-1),(u/3,-1),(-u/3,-1),(u/2,1.5),(-u/2,1.5),(u/4,.25),(-u/4,.25)]
cube=[frozenset(b[2]) for b in fblocks(P8).values()]
print("cube blocks",len(cube),sorted(len(b) for b in cube))
n=9; SMAX=5; TARGET=int(sys.argv[1]) if len(sys.argv)>1 else 24
blocks=[frozenset(c) for s in range(3,SMAX+1) for c in itertools.combinations(range(n),s)]
idx={b:i for i,b in enumerate(blocks)}
m=cp_model.CpModel()
x=[m.NewBoolVar("") for _ in blocks]; y=[m.NewBoolVar("") for _ in blocks]
for t in itertools.combinations(range(n),3):
    m.AddExactlyOne([x[i] for i,b in enumerate(blocks) if set(t)<=b])
for b in cube:
    opts=[idx[b]]+([idx[b|{8}]] if len(b)+1<=SMAX else [])
    m.AddExactlyOne([x[i] for i in opts])
for i in range(len(blocks)): m.AddImplication(y[i],x[i])
for i,j in itertools.combinations(range(len(blocks)),2):
    if len(blocks[i]&blocks[j])>=2: m.AddBoolOr([y[i].Not(),y[j].Not()])
m.Add(sum(y)<=10)
m.Minimize(sum(x)-sum(y))
s=cp_model.CpSolver(); s.parameters.num_workers=1; s.parameters.max_time_in_seconds=600
r=s.Solve(m); print(s.StatusName(r), "min abstract circles (cube+1 designs) =", s.ObjectiveValue() if r in (cp_model.OPTIMAL,cp_model.FEASIBLE) else None)
if r in (cp_model.OPTIMAL,cp_model.FEASIBLE):
    xb=[sorted(blocks[i]) for i in range(len(blocks)) if s.Value(x[i])]
    yb=[sorted(blocks[i]) for i in range(len(blocks)) if s.Value(y[i])]
    print("B",len(xb),"L",len(yb)); print("extended",[b for b in xb if 8 in b and len(b)>3]); print("lines",yb)
