import sympy as sp
from itertools import product
t=sp.symbols('t')
# coefficient covectors (a,b,c) for a x + b y + c z = 0
A=[('L1',(0,1,-t)),('L2',(0,-t,1))]
B=[('L3',(-t,0,1)),('L4',(1,0,-t))]
C=[('L5',(1,-t,0)),('L6',(-t,1,0))]
print("Transversal triples: concurrence determinant det[coeff rows]\n")
seen={}
for (na,va),(nb,vb),(nc,vc) in product(A,B,C):
    M=sp.Matrix([va,vb,vc]); d=sp.factor(M.det())
    print(f"  {na},{nb},{nc}:  det = {d}")
    seen.setdefault(str(d),[]).append((na,nb,nc))
print("\nDistinct determinant values and which triples give them:")
for d,tr in seen.items():
    print(f"  {d}   <- {len(tr)} triples")
print("\nRoots check: t^2+t+1 discriminant =",sp.discriminant(t**2+t+1,t),"(negative => no real roots); so positive real root only at t=1.")
# sample explicit 3x3 for the write-up: L1,L3,L5
M=sp.Matrix([(0,1,-t),(-t,0,1),(1,-t,0)])
print("\nSample L1,L3,L5 matrix:")
sp.pprint(M)
print("det =",sp.factor(M.det()))
