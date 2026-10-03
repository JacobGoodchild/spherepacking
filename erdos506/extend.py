import itertools, numpy as np, sys, time
from math import sqrt
EPS=1e-7
def blocks(P):
    """Return dict key->(type, params, set of point idx). type 'C' (cx,cy,r) or 'L' (a,b,c) normalized."""
    B={}
    for i,j,k in itertools.combinations(range(len(P)),3):
        (x1,y1),(x2,y2),(x3,y3)=P[i],P[j],P[k]
        D=2*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
        if abs(D)<1e-9:
            a,b=y2-y1,x1-x2; n=np.hypot(a,b); a,b=a/n,b/n; c=-(a*x1+b*y1)
            if a<-1e-12 or (abs(a)<1e-12 and b<0): a,b,c=-a,-b,-c
            key=('L',round(a,6),round(b,6),round(c,6)); par=(a,b,c)
        else:
            s1,s2,s3=x1*x1+y1*y1,x2*x2+y2*y2,x3*x3+y3*y3
            cx=(s1*(y2-y3)+s2*(y3-y1)+s3*(y1-y2))/D; cy=(s1*(x3-x2)+s2*(x1-x3)+s3*(x2-x1))/D
            r=np.hypot(x1-cx,y1-cy); key=('C',round(cx,6),round(cy,6),round(r,6)); par=(cx,cy,r)
        if key not in B: B[key]=[key[0],par,set()]
        B[key][2].update((i,j,k))
    return B
def count(P):
    B=blocks(P); return sum(1 for b in B.values() if b[0]=='C'), B
def inter(o1,o2):
    t1,p1=o1; t2,p2=o2
    if t1=='L' and t2=='L':
        a1,b1,c1=p1; a2,b2,c2=p2; d=a1*b2-a2*b1
        if abs(d)<1e-12: return []
        return [((b1*c2-b2*c1)/d,(c1*a2-c2*a1)/d)]
    if t1=='C' and t2=='L': o1,o2=o2,o1; t1,p1,t2,p2=t2,p2,t1,p1
    if t1=='L':
        a,b,c=p1; cx,cy,r=p2; dist=a*cx+b*cy+c
        if abs(dist)>r+1e-12: return []
        fx,fy=cx-a*dist,cy-b*dist; h=sqrt(max(r*r-dist*dist,0))
        return [(fx+b*h,fy-a*h),(fx-b*h,fy+a*h)]
    (x1,y1,r1),(x2,y2,r2)=p1,p2; dx,dy=x2-x1,y2-y1; d=np.hypot(dx,dy)
    if d<1e-12 or d>r1+r2+1e-12 or d<abs(r1-r2)-1e-12: return []
    a=(r1*r1-r2*r2+d*d)/(2*d); h=sqrt(max(r1*r1-a*a,0)); mx,my=x1+a*dx/d,y1+a*dy/d
    return [(mx+h*dy/d,my-h*dx/d),(mx-h*dy/d,my+h*dx/d)]
def candidates(P,B):
    objs=[(b[0],b[1]) for b in B.values()]
    for i,j in itertools.combinations(range(len(P)),2):   # all lines through pairs (free)
        (x1,y1),(x2,y2)=P[i],P[j]; a,b=y2-y1,x1-x2; n=np.hypot(a,b); a,b=a/n,b/n
        objs.append(('L',(a,b,-(a*x1+b*y1))))
    pts={}
    for o1,o2 in itertools.combinations(objs,2):
        for q in inter(o1,o2):
            if abs(q[0])>1e4 or abs(q[1])>1e4: continue
            if min(np.hypot(q[0]-x,q[1]-y) for x,y in P)<1e-6: continue
            pts[(round(q[0],6),round(q[1],6))]=q
    return list(pts.values())
def beam(P0,target_n,width=30,keep=5):
    level=[(count(P0)[0],P0)]
    while len(level[0][1])<target_n:
        nxt={}
        for c,P in level:
            _,B=count(P)
            for q in candidates(P,B):
                P2=P+[q]; c2,B2=count(P2)
                if any(len(b[2])==len(P2) for b in B2.values()): continue
                key=tuple(sorted(len(b[2]) for b in B2.values()))+(c2,)
                if key not in nxt or c2<nxt[key][0]: nxt[key]=(c2,P2)
        level=sorted(nxt.values(),key=lambda t:t[0])[:width]
        print(f"  n={len(level[0][1])}: best {level[0][0]}  ({[t[0] for t in level[:8]]})",flush=True)
    return level
if __name__=="__main__":
    u=sqrt(15)
    seeds={
     '8pt_17':[(u,-1),(-u,-1),(u/3,-1),(-u/3,-1),(u/2,1.5),(-u/2,1.5),(u/4,.25),(-u/4,.25)],
     '8pt_18':[(1,1),(1,-1),(-1,1),(-1,-1),(2,2),(2,-2),(-2,2),(-2,-2)],
     '7pt_11':[(0,0),(175,0),(-225,300),(-225,-300),(-225,0),(63,84),(63,-84)],
     '6pt_8':[(0,0),(25,0),(-15,0),(-15,-20),(5,-10),(9,12)],
    }
    for name in sys.argv[1:]:
        P=seeds[name]; print(name,"start circles",count(P)[0],flush=True)
        res=beam(P,9)
        print(name,"=> best n=9:",res[0][0],res[0][1],flush=True)
