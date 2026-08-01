import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
FIGDIR = _os.path.abspath(_os.path.join(_HERE, '..', 'figures'))
_os.makedirs(FIGDIR, exist_ok=True)
def figpath(name): return _os.path.join(FIGDIR, name)
def datapath(name): return _os.path.join(_HERE, name)

from fractions import Fraction as F
from itertools import combinations
import numpy as np
from shapely.geometry import Polygon, LineString
from shapely.ops import polygonize, unary_union

# ---- refined palette ----
COL={'cornerT':'#5B8FD6','innerT':'#48B0A5','quad':'#74C078','pent':'#E58C6A','hex':'#D9534F'}
ACCENT='#555555'  # neutral accent
LAB={'cornerT':'corner triangle (6)','innerT':'inner triangle (6)','quad':'quadrilateral (3)',
     'pent':'pentagon (3)','hex':'hexagon (1)'}
SIDES_N={'cornerT':3,'innerT':3,'quad':4,'pent':5,'hex':6}

def cevian_lines(t):
    return {'A1':(F(1),-t,F(0)),'A2':(-t,F(1),F(0)),'B1':(t,1+t,-t),
            'B2':(F(1),1+t,F(-1)),'C1':(1+t,F(1),F(-1)),'C2':(1+t,t,-t)}
SIDELINES=[(F(0),F(1),F(0)),(F(1),F(0),F(0)),(F(1),F(1),F(-1))]
def inter(l1,l2):
    a1,b1,c1=l1;a2,b2,c2=l2;d=a1*b2-a2*b1
    if d==0:return None
    return ((-c1*b2+c2*b1)/d,(-a1*c2+a2*c1)/d)
def inside(p):X,Y=p;return X>=0 and Y>=0 and X+Y<=1
def on(l,p):a,b,c=l;return a*p[0]+b*p[1]+c==0
TRI=Polygon([(0,0),(1,0),(0,1)])

def arrangement(keys,t):
    L=cevian_lines(t); lines=[L[k] for k in keys]+SIDELINES
    nodes=set()
    for l1,l2 in combinations(lines,2):
        p=inter(l1,l2)
        if p and inside(p):nodes.add(p)
    nodes=list(nodes); edges=[]
    for l in lines:
        pts=[p for p in nodes if on(l,p)]
        pts.sort(key=(lambda p:(p[0],p[1])) if l[1]!=0 else (lambda p:(p[1],p[0])))
        for i in range(len(pts)-1):edges.append((pts[i],pts[i+1]))
    fl=[LineString([(float(a[0]),float(a[1])),(float(b[0]),float(b[1]))]) for a,b in edges]
    fs=[f for f in polygonize(unary_union(fl)) if f.representative_point().within(TRI)]
    nf=[(float(p[0]),float(p[1]),p) for p in nodes]
    def nearest(pt):return min(nf,key=lambda z:(z[0]-pt[0])**2+(z[1]-pt[1])**2)[2]
    def corners(face):
        cs=[nearest(c) for c in list(face.exterior.coords)[:-1]];cc=[]
        for p in cs:
            if not cc or cc[-1]!=p:cc.append(p)
        if len(cc)>1 and cc[0]==cc[-1]:cc.pop()
        n=len(cc);k=0
        for i in range(n):
            p0=cc[i-1];p1=cc[i];p2=cc[(i+1)%n]
            cr=(p1[0]-p0[0])*(p2[1]-p1[1])-(p1[1]-p0[1])*(p2[0]-p1[0])
            if cr!=0:k+=1
        return k
    return [(corners(f),F(f.area).limit_denominator(10**9)*2,f) for f in fs]

def classify(ng,ratio,t):
    if ng==6:return 'hex'
    if ng==5:return 'pent'
    if ng==4:return 'quad'
    corner=F(1)/((t+1)*(t*t+t+1))
    return 'cornerT' if ratio==corner else 'innerT'

def to_eq(pt):
    X,Y=float(pt[0]),float(pt[1])
    return (X+0.5*Y,(np.sqrt(3)/2)*Y)

# equilateral vertices in display coords
EA=to_eq((0,0)); EB=to_eq((1,0)); EC=to_eq((0,1))
print("common loaded; palette:",COL)
