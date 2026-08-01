from fractions import Fraction as F
from collections import Counter
from common import *
t=F(2); order=['A1','A2','B1','B2','C1','C2']
NG={3:'T',4:'Q',5:'P',6:'H'}
prev=None
state=None
import shapely
# reuse arrangement; track parent splits per stage
def faces(keys): return arrangement(keys,t)
cur=faces([]); state=[(0,cur[0][2])]; 
lab=['+A1','+A2','+B1','+B2','+C1','+C2']
for si,key in enumerate(order,start=1):
    new=faces(order[:si]); groups={}
    for ng,ratio,f in new:
        rp=f.representative_point()
        for pid,pp in state:
            if pp.contains(rp): groups.setdefault(pid,[]).append((ng,f));break
    rules=Counter(); ns=[]; nid=max(p for p,_ in state)+1
    for pid,pp in state:
        ch=groups.get(pid,[])
        if len(ch)>=2:
            png=None
            # parent ngon: recompute from prev arrangement not stored; infer as sum-1 of children? 
            # Instead: parent polygon corner count
            from shapely.geometry import Polygon
            # count parent corners via its polygon
            # (approx via number of exterior coords-1 after dedup is unreliable; use children types)
            child_ngs=tuple(sorted(c[0] for c in ch))
            rules[child_ngs]+=1
            for ng,f in ch: ns.append((nid,f)); nid+=1
        else:
            ns.append((pid,ch[0][1] if ch else pp))
    state=ns
    census=Counter(); 
    for _,f in state: pass
    cen=Counter(ng for ng,_,_ in new)
    rulestr="; ".join(f"{cnt}×(→{'+'.join(NG[x] for x in ch)})" for ch,cnt in sorted(rules.items()))
    censusstr=" ".join(f"{cen[k]}{NG[k]}" for k in sorted(cen))
    print(f"{lab[si-1]:5s}: splits {rulestr:35s} | census {censusstr}")
