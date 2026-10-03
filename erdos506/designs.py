import itertools, sys, time
from ortools.sat.python import cp_model
n=9; TARGET=int(sys.argv[1]) if len(sys.argv)>1 else 24; SMAX=int(sys.argv[2]) if len(sys.argv)>2 else 5
blocks=[frozenset(c) for s in range(3,SMAX+1) for c in itertools.combinations(range(n),s)]
triples=list(itertools.combinations(range(n),3))
def model(extra=None, seed=0):
    m=cp_model.CpModel()
    x=[m.NewBoolVar(f"x{i}") for i in range(len(blocks))]
    y=[m.NewBoolVar(f"y{i}") for i in range(len(blocks))]
    for t in triples:
        m.AddExactlyOne([x[i] for i,b in enumerate(blocks) if set(t)<=b])
    for i in range(len(blocks)): m.AddImplication(y[i],x[i])
    # line blocks (through the pole) pairwise share <= 1 point
    for i,j in itertools.combinations(range(len(blocks)),2):
        if len(blocks[i]&blocks[j])>=2: m.AddBoolOr([y[i].Not(),y[j].Not()])
    m.Add(sum(x)-sum(y)<=TARGET)
    # symmetry breaking (weak): point 0..2 triple's block contains 0,1,2 and is largest-ish
    return m,x,y
if __name__=="__main__":
    m,x,y=model()
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=float(sys.argv[3]) if len(sys.argv)>3 else 120; s.parameters.num_workers=4
    t0=time.time(); r=s.Solve(m)
    print("status",s.StatusName(r),f"{time.time()-t0:.1f}s")
    if r in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        B=[sorted(blocks[i]) for i in range(len(blocks)) if s.Value(x[i])]
        L=[sorted(blocks[i]) for i in range(len(blocks)) if s.Value(y[i])]
        print("B",len(B),"L",len(L),"circles",len(B)-len(L))
        print("sizes",sorted(len(b) for b in B))
        print("blocks",B); print("lines",L)
