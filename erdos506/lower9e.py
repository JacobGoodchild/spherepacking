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
# No Fano plane in any derived structure (Fano is not realizable over R).
# z[q][abc]: the block through q,a,b also contains c (i.e. some chosen block contains {q,a,b,c}).
FANO=[(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]
import itertools as it
fanos=set()
for perm in it.permutations(range(7)):
    fanos.add(frozenset(frozenset(perm[i] for i in L) for L in FANO))
fanos=list(fanos)   # 30 labelled Fano planes on 7 symbols
zc={}
def z(q,trip):
    key=(q,frozenset(trip))
    if key not in zc:
        zc[key]=sum(x[i] for i,b in enumerate(blocks) if q in b and set(trip)<=b)
    return zc[key]
cnt=0
for q in range(n):
    others=[p for p in range(n) if p!=q]
    for S in it.combinations(others,7):
        for F in fanos:
            m.Add(sum(z(q,[S[i] for i in L]) for L in F) <= 6); cnt+=1
# also at the pole: the real lines cannot contain a Fano plane
yc={}
def zy(trip):
    key=frozenset(trip)
    if key not in yc: yc[key]=sum(y[i] for i,b in enumerate(blocks) if set(trip)<=b)
    return yc[key]
for S in it.combinations(range(n),7):
    for F in fanos:
        m.Add(sum(zy([S[i] for i in L]) for L in F) <= 6); cnt+=1
print("Fano exclusions",cnt,flush=True)
def labelled(conf,k):
    out=set()
    for perm in it.permutations(range(k)):
        out.add(frozenset(frozenset(perm[i] for i in L) for L in conf))
    return [list(map(sorted,c)) for c in out]
MK=[(0,1,3),(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,0),(6,7,1),(7,0,2)]
mks=labelled(MK,8)
c2=0
for q in range(n):       # derived structure at q lives on the 8 other points
    others=[p for p in range(n) if p!=q]
    for F in mks:
        m.Add(sum(z(q,[others[i] for i in L]) for L in F) <= 7); c2+=1
for S in it.combinations(range(n),8):   # real lines at the pole
    for F in mks:
        m.Add(sum(zy([S[i] for i in L]) for L in F) <= 7); c2+=1
print("Mobius-Kantor exclusions",c2,flush=True)
PAP=[(0,1,2),(3,4,5),(6,7,8),(0,4,8),(0,5,7),(1,3,8),(1,5,6),(2,3,7),(2,4,6)]
paps=labelled(PAP,9); c3=0
for F in paps:            # Pappus closure for real lines: 8 lines force the 9th
    zs=[zy(L) for L in F]
    for j in range(9):
        m.Add(sum(zs[k] for k in range(9) if k!=j) - zs[j] <= 7); c3+=1
print("Pappus closures",c3,flush=True)
# Miquel closure: 8 points on cube vertices; 5 concyclic faces force the 6th.
CUBEF=[(0,1,2,3),(4,5,6,7),(0,1,4,5),(2,3,6,7),(0,2,4,6),(1,3,5,7)]
def on_block(Q):
    return sum(x[i] for i,b in enumerate(blocks) if set(Q)<=b)
qc={}
def ob(Q):
    k=frozenset(Q)
    if k not in qc: qc[k]=on_block(Q)
    return qc[k]
cubes=set()
for perm in it.permutations(range(8)):
    cubes.add(frozenset(frozenset(perm[i] for i in f) for f in CUBEF))
cubes=[list(map(sorted,c)) for c in cubes]
c4=0
for S in it.combinations(range(n),8):
    for C in cubes:
        fs=[ob([S[i] for i in f]) for f in C]
        for j in range(6):
            m.Add(sum(fs[k] for k in range(6) if k!=j) - fs[j] <= 4); c4+=1
print("Miquel closures",c4,flush=True)
m.Add(sum(x)-sum(y)<=TARGET)
# symmetry breaking: label points by non-increasing number of big (>=4) blocks through them
deg=[sum(x[i] for i,b in enumerate(blocks) if q in b and len(b)>=4) for q in range(n)]
for q in range(n-1): m.Add(deg[q]>=deg[q+1])
LMIN=int(sys.argv[7]); LMAX=int(sys.argv[8])
m.Add(sum(y)>=LMIN); m.Add(sum(y)<=LMAX)
print("L range",LMIN,LMAX,flush=True)
s=cp_model.CpSolver(); s.parameters.num_workers=int(sys.argv[6]) if len(sys.argv)>6 else 1; s.parameters.max_time_in_seconds=T
t0=time.time(); r=s.Solve(m)
print(f"target<= {TARGET}, blocks<= {SMAX}, C3>={C3}, C4>={C4}: {s.StatusName(r)} ({time.time()-t0:.1f}s)")
if r in (cp_model.OPTIMAL,cp_model.FEASIBLE):
    xb=[sorted(blocks[i]) for i in range(len(blocks)) if s.Value(x[i])]; yb=[sorted(blocks[i]) for i in range(len(blocks)) if s.Value(y[i])]
    print(" B",len(xb),"L",len(yb),"sizes",sorted(len(b) for b in xb)); print(" lines",yb); print(" blocks",xb)
