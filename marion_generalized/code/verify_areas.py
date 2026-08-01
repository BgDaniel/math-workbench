import sympy as sp
t=sp.symbols('t')

# normalized barycentric vertices; signed area ratio of a triangle = det of the 3x3 matrix of rows
def bary(*pts):
    return sp.Matrix([list(p) for p in pts])
def sarea_tri(P,Q,R):  # signed area ratio of triangle PQR relative to ABC
    return bary(P,Q,R).det()
def norm(P):
    s=sum(P); return tuple(sp.simplify(c/s) for c in P)

A=(sp.Integer(1),0,0); B=(0,sp.Integer(1),0); C=(0,0,sp.Integer(1))
G=(sp.Rational(1,3),)*3

# six cevian lines as (coeff on x,y,z) linear forms =0
L={'L1':(0,1,-t),   # y=tz
   'L2':(0,-t,1),   # z=ty
   'L3':(-t,0,1),   # z=tx
   'L4':(1,0,-t),   # x=tz
   'L5':(1,-t,0),   # x=ty
   'L6':(-t,1,0)}   # y=tx
SIDE={'AB':(0,0,1),'BC':(1,0,0),'CA':(0,1,0)}  # z=0, x=0, y=0

def meet(l1,l2):
    a=sp.Matrix([l1,l2])
    ns=a.nullspace()[0]
    return norm(tuple(ns))

# ---- explicit corner triangle at A (side AB, lines L5, L1) ----
P1=A
P2=meet(SIDE['AB'],L['L5'])      # foot on AB
P3=meet(L['L5'],L['L1'])         # inner double point
corner=sp.simplify(sarea_tri(P1,P2,P3))
print("Corner triangle vertices:")
print("  A         =",P1)
print("  AB ∩ L5   =",P2)
print("  L5 ∩ L1   =",P3)
print("  signed area ratio =",corner,"   [closed form 1/((t+1)(t^2+t+1))]")
print("  matches:",sp.simplify(corner-1/((t+1)*(t**2+t+1)))==0)
print()

# ---- hexagon: its six vertices are the 'inner' L_i∩L_j; get them, fan from centroid ----
# hexagon bounded by all six cevians; vertices = consecutive cevian intersections
hexpairs=[('L1','L4'),('L4','L5'),('L5','L2'),('L2','L3'),('L3','L6'),('L6','L1')]
V=[meet(L[a],L[b]) for a,b in hexpairs]
# signed area via fan from centroid G
hexarea=sp.simplify(sum(sarea_tri(G,V[i],V[(i+1)%6]) for i in range(6)))
print("Hexagon vertices (L_i∩L_j):")
for (a,b),v in zip(hexpairs,V): print(f"  {a}∩{b} = {v}")
print("  signed area ratio =",hexarea,"  [closed form 2(t-1)^2/((t+2)(2t+1))]")
print("  matches:",sp.simplify(hexarea-2*(t-1)**2/((t+2)*(2*t+1)))==0)
print()

# ---- signed values at t=2 and t=1/2, with multiplicities; verify signed sum=1 ----
forms={'Hex':2*(t-1)**2/((t+2)*(2*t+1)),
       'Inner':t*(t-1)**2/((t+2)*(2*t+1)*(t**2+t+1)),
       'Corner':1/((t+1)*(t**2+t+1)),
       'Quad':2*(t-1)/((t+2)*(t**2+t+1)),
       'Pent':(t-1)*(t**2+3*t+1)/((t+1)*(2*t+1)*(t**2+t+1))}
mult={'Hex':1,'Inner':6,'Corner':6,'Quad':3,'Pent':3}
for tv in [sp.Integer(2),sp.Rational(1,2)]:
    print(f"--- signed formula values at t={tv} ---")
    tot=0
    for k,f in forms.items():
        val=sp.nsimplify(f.subs(t,tv)); tot+=mult[k]*val
        sign='+' if val>=0 else '−'
        print(f"  {k:7s} x{mult[k]}: {str(val):>10s}   (each {sign})   group total {sp.nsimplify(mult[k]*val)}")
    print(f"  SIGNED SUM = {sp.nsimplify(tot)}")
    print()
