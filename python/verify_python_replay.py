#!/usr/bin/env python3
"""Rigorous fixed-point interval replay of the vector Bernstein certificate.

The exact integer chart polynomials are converted to common degree four
Bernstein control nets.  Each component is normalized by its exact maximum
absolute root control coefficient.  Endpoints are enclosed on a 2^-52 grid.
Dyadic de Casteljau subdivision is then performed with directed integer
floor/ceiling, so every true control coefficient remains enclosed.

A leaf is removed only when an IEEE-double (hence dyadic-rational) vector u
has a rigorously positive dot product with *every* interval control vector;
an explicit roundoff allowance of 64 machine eps times the absolute sum is
subtracted.  Thus an exhausted tree certifies that zero is outside the image
of all five polynomials on the cube.  Integer overflow is guarded explicitly.
"""
import argparse,json,math,time
import numpy as np
from pathlib import Path

from bernstein_tools import (DEG, T, to_bernstein, split_half,
                             midpoint_value, split_score)

S=1<<52
EPS=np.finfo(float).eps

def exact_root(raw):
    # 6*T is the exact integer conversion matrix power(x on [-1,1]) ->
    # degree-4 Bernstein.  After six axes all entries have denominator 6^6.
    t6=np.rint(6*T).astype(np.int64)
    assert np.max(abs(6*T-t6))==0
    roots=[]
    for poly in raw:
        a=np.zeros((5,)*6,dtype=np.int64)
        for m,c in poly:a[tuple(m)]+=int(c)
        for axis in range(6):
            a=np.moveaxis(np.tensordot(t6,a,axes=(1,axis)),0,axis)
        roots.append(a)
    mids=[];los=[];his=[]
    for a,poly in zip(roots,raw):
        den=int(np.max(np.abs(a)));assert den>0
        flat=a.ravel();lo=np.fromiter(((int(v)*S)//den for v in flat),dtype=np.int64,count=flat.size)
        hi=np.fromiter((-((-int(v)*S)//den) for v in flat),dtype=np.int64,count=flat.size)
        lo=lo.reshape(a.shape);hi=hi.reshape(a.shape)
        los.append(lo);his.append(hi);mids.append(to_bernstein(poly))
    return mids,los,his

# Integer numerator matrices for half subdivision, all row sums 16.
L=np.zeros((5,5),dtype=np.int64);R=np.zeros((5,5),dtype=np.int64)
for i in range(5):
    for k in range(i+1):L[i,k]=(1<<(4-i))*math.comb(i,k)
for i in range(5):
    # right_i = sum_{k=i}^4 C(4-i,k-i)/2^(4-i) b_k
    for k in range(i,5):R[i,k]=(1<<i)*math.comb(4-i,k-i)
assert np.all(L.sum(1)==16) and np.all(R.sum(1)==16)

def divlo(x):return x//16
def divhi(x):return -((-x)//16)
def split_enclosure(arr,axis,lower):
    x=np.moveaxis(arr,axis,-1)
    assert int(np.max(np.abs(x))) < (1<<58)
    l=np.tensordot(x,L.T,axes=(-1,0));r=np.tensordot(x,R.T,axes=(-1,0))
    op=divlo if lower else divhi
    return np.moveaxis(op(l),-1,axis),np.moveaxis(op(r),-1,axis)

def candidate_direction(polys,iterations=4):
    # Scalar Bernstein signs first.
    for i,p in enumerate(polys):
        if float(np.min(p))>0:
            u=np.zeros(5);u[i]=1;return float(np.min(p)),u
        if float(np.max(p))<0:
            u=np.zeros(5);u[i]=-1;return float(-np.max(p)),u
    y=np.asarray([midpoint_value(p) for p in polys]);best=(-np.inf,None)
    if np.linalg.norm(y)==0:return best
    for _ in range(iterations+1):
        combo=sum(z*p for z,p in zip(y,polys));idx=int(np.argmin(combo))
        mdot=float(combo.flat[idx]);margin=mdot/np.linalg.norm(y)
        if margin>best[0]:best=(margin,y.copy())
        ix=np.unravel_index(idx,combo.shape);v=np.asarray([p[ix] for p in polys]);d=v-y;dd=float(d@d)
        if dd==0:break
        gamma=float(np.clip(-(y@d)/dd,0,1))
        if gamma==0:break
        y=y+gamma*d
        if np.linalg.norm(y)==0:break
    return best

def exact_positive(u,los,his):
    terms=[]
    for z,lo,hi in zip(u,los,his):
        a=lo if z>=0 else hi
        terms.append(float(z)*(a.astype(np.float64)/S))
    val=sum(terms);ab=sum(abs(x) for x in terms)
    # Each integer/S conversion is exact (|integer|<2^53); allow far more
    # than the five products/four additions actually used.
    lower=val-64*EPS*ab-np.nextafter(0.,1.)
    return float(np.min(lower)),float(np.max(ab))

def run(source, maxnodes=1_000_000, maxdepth=100, report=1000):
    raw = json.loads(Path(source).read_text(encoding="utf-8"))["polynomials"]
    ps,lo,hi=exact_root(raw)
    stack=[(np.zeros(6,dtype=np.int16),ps,lo,hi,
            np.zeros(6),np.ones(6))]
    nodes=removed=deep=failed=0;worst=float('inf');t=time.time()
    while stack and nodes<maxnodes:
        ds,p,l,h,blo,bhi=stack.pop();nodes+=1;margin,u=candidate_direction(p)
        certified=False
        # This exact fixed-point scalar test is both cheaper and more robust
        # than consulting the floating scout.  It is decisive on the one
        # exceptionally deep boundary branch where repeatedly subdivided
        # float control nets lose useful componentwise range information.
        # The ordinary vector certificate closes all normal boxes below
        # depth 50.  Delay this more expensive full integer range scan until
        # depth 60; it is needed only for the exceptional float-wrapping
        # branch near a chart boundary.
        if int(ds.sum())>=60:
            scalar=[]
            for i in range(5):
                scalar.append(max(int(np.min(l[i])), -int(np.max(h[i])))/S)
            if max(scalar)>0:
                em=max(scalar);certified=True
        if margin>0 and not certified:
            em,ab=exact_positive(u,l,h)
            certified=em>0
            if not certified:
                margin,u=candidate_direction(p,iterations=32)
                if margin>0:
                    em,ab=exact_positive(u,l,h)
                    certified=em>0
        if certified:
            removed+=1;worst=min(worst,em);continue
        if margin>0:failed+=1
        if int(ds.sum())>=maxdepth:
            deep+=1
            print('UNRESOLVED_BOX depths',ds.tolist(),'lo',blo.tolist(),
                  'hi',bhi.tolist(),'float_margin',margin,flush=True)
            continue
        axis=int(np.argmax([split_score(p,i) for i in range(6)]))
        fp=[split_half(z,axis) for z in p]
        il=[split_enclosure(z,axis,True) for z in l]
        ih=[split_enclosure(z,axis,False) for z in h]
        nd=ds.copy();nd[axis]+=1
        mid=(blo[axis]+bhi[axis])/2
        lhi=bhi.copy();lhi[axis]=mid
        rlo=blo.copy();rlo[axis]=mid
        stack.append((nd.copy(),[z[1] for z in fp],[z[1] for z in il],
                      [z[1] for z in ih],rlo,bhi.copy()))
        stack.append((nd,[z[0] for z in fp],[z[0] for z in il],
                      [z[0] for z in ih],blo.copy(),lhi))
        if report and nodes%report==0:
            print('PROGRESS',nodes,'stack',len(stack),'removed',removed,'failed',failed,
                  'deep',deep,'depth',int(ds.sum()),'worst',worst,
                  'sec',round(time.time()-t,1),flush=True)
    print('FINAL nodes',nodes,'stack',len(stack),'removed',removed,'failed',failed,
          'deep',deep,'worst_exact_margin',worst,'seconds',time.time()-t)
    if not stack and deep==0:
        print('CERTIFIED: exact interval Bernstein separation excludes a common zero on [-1,1]^6')
        return True
    print('INCOMPLETE')
    return False

if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path,
                    default=here.parent / "data" / "final_chart_polynomials.json")
    ap.add_argument("--maxnodes", type=int, default=1_000_000)
    ap.add_argument("--maxdepth", type=int, default=100)
    ap.add_argument("--report", type=int, default=1000)
    args = ap.parse_args()
    raise SystemExit(0 if run(args.source, args.maxnodes, args.maxdepth, args.report) else 1)
