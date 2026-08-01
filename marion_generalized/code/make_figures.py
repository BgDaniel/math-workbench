"""
make_figures.py — generate every figure for the paper in one go.

Run:
    python make_figures.py            # from anywhere; writes PNGs into ../figures/

Each figure is a function; main() runs them in dependency order. The shared geometry
(arrangement, classification, palette, path helpers) lives in common.py.
"""
import json
from fractions import Fraction as F
from itertools import combinations
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly, FancyArrowPatch
from matplotlib.lines import Line2D
import sympy as sp
from common import *   # arrangement, classify, to_eq, COL, LAB, EA/EB/EC, SIDELINES,
                       # cevian_lines, inter, on, inside, figpath, datapath, FIGDIR

plt.rcParams.update({'font.family': 'DejaVu Sans'})

def extract_representatives():

    t2=F(2)
    # exact nodes at t=2
    L=cevian_lines(t2); alllines={**L}
    LINELABELS={'A1':L['A1'],'A2':L['A2'],'B1':L['B1'],'B2':L['B2'],'C1':L['C1'],'C2':L['C2'],
                'AB':SIDELINES[0],'CA':SIDELINES[1],'BC':SIDELINES[2]}
    lines=list(LINELABELS.values())
    nodes=set()
    for l1,l2 in combinations(lines,2):
        p=inter(l1,l2)
        if p and inside(p): nodes.add(p)
    nodes=list(nodes)
    def lines_through(p):
        return [name for name,l in LINELABELS.items() if on(l,p)]

    faces=arrangement(['A1','A2','B1','B2','C1','C2'],t2)
    def nearest_node(pt):
        return min(nodes,key=lambda z:(float(z[0])-pt[0])**2+(float(z[1])-pt[1])**2)

    reps={}
    for ng,ratio,f in faces:
        key=classify(ng,ratio,t2)
        if key in reps: continue
        coords=[nearest_node(c) for c in list(f.exterior.coords)[:-1]]
        # dedup consecutive
        cc=[]
        for p in coords:
            if not cc or cc[-1]!=p: cc.append(p)
        if len(cc)>1 and cc[0]==cc[-1]: cc.pop()
        desc=[tuple(sorted(lines_through(p))[:2]) for p in cc]
        reps[key]=(cc,desc,ratio)

    # symbolic lines in (X,Y) as function of t
    T=sp.symbols('t')
    def symline(name):
        a,b,c={'A1':(1,-T,0),'A2':(-T,1,0),'B1':(T,1+T,-T),'B2':(1,1+T,-1),
               'C1':(1+T,1,-1),'C2':(1+T,T,-T),'AB':(0,1,0),'CA':(1,0,0),'BC':(1,1,-1)}[name]
        return (sp.Integer(a) if not hasattr(a,'free_symbols') else a,
                sp.sympify(b),sp.sympify(c))
    def symmeet(n1,n2):
        a1,b1,c1=symline(n1); a2,b2,c2=symline(n2)
        d=a1*b2-a2*b1
        X=sp.simplify((-c1*b2+c2*b1)/d); Y=sp.simplify((-a1*c2+a2*c1)/d)
        return (X,Y)
    def signed_ratio(descs):
        V=[symmeet(*d) for d in descs]
        s=0
        n=len(V)
        for i in range(n):
            X1,Y1=V[i];X2,Y2=V[(i+1)%n]; s+=X1*Y2-X2*Y1
        return sp.simplify(s)  # = 2*signed area in (X,Y) = signed area ratio (since tri area=1/2)

    print("Representatives (vertex line-pairs), signed area ratio at t=2 and t=1/2:\n")
    order=['hex','innerT','cornerT','quad','pent']
    DESC={}
    for key in order:
        cc,desc,ratio=reps[key]
        sr=signed_ratio(desc)
        v2=sp.nsimplify(sr.subs(T,2)); vh=sp.nsimplify(sr.subs(T,sp.Rational(1,2)))
        # normalize orientation to + at t=2
        if v2<0: desc=desc[::-1]; sr=-sr; v2=-v2; vh=-vh
        DESC[key]=desc
        print(f"{key:8s} verts={desc}")
        print(f"         signed ratio(t) = {sr}")
        print(f"         t=2 -> {v2}    t=1/2 -> {vh}\n")

    json.dump({k:[list(d) for d in v] for k,v in DESC.items()}, open(datapath('representatives.json'),'w'))
    print("saved representatives.json")
    return DESC


def fig_family():

    plt.rcParams.update({'font.family':'DejaVu Sans'})

    def draw(ax,t,title):
        tr=F(t).limit_denominator(10**6)
        for ng,ratio,f in arrangement(['A1','A2','B1','B2','C1','C2'],tr):
            key=classify(ng,ratio,tr)
            pts=[to_eq(c) for c in list(f.exterior.coords)[:-1]]
            ax.add_patch(MplPoly(pts,closed=True,facecolor=COL[key],edgecolor='white',lw=1.1))
        ax.add_patch(MplPoly([EA,EB,EC],closed=True,fill=False,edgecolor='#222',lw=1.8))
        ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=12,color='#222',pad=6)

    def draw_medians(ax):
        A=np.array(EA);B=np.array(EB);C=np.array(EC)
        Mab=(A+B)/2;Mbc=(B+C)/2;Mca=(C+A)/2;G=(A+B+C)/3
        for tp in [(A,Mab,G),(Mab,B,G),(B,Mbc,G),(Mbc,C,G),(C,Mca,G),(Mca,A,G)]:
            ax.add_patch(MplPoly([tuple(p) for p in tp],closed=True,facecolor=COL['cornerT'],edgecolor='white',lw=1.1))
        ax.add_patch(MplPoly([tuple(A),tuple(B),tuple(C)],closed=True,fill=False,edgecolor='#222',lw=1.8))
        ax.set_aspect('equal');ax.axis('off');ax.set_title(r"$t=1$: medians",fontsize=12,pad=6)

    handles=[MplPoly([(0,0)],facecolor=COL[k],edgecolor='white',label=LAB[k]) for k in ['cornerT','innerT','quad','pent','hex']]

    # ---- panels ----
    fig,axs=plt.subplots(1,5,figsize=(16,3.7))
    fig.patch.set_facecolor('white')
    draw_medians(axs[0]); draw(axs[1],F(3,2),r"$t=3/2$"); draw(axs[2],F(2),r"$t=2$ (Marion)")
    draw(axs[3],F(4),r"$t=4$"); draw(axs[4],F(12),r"$t=12\ (\to\ \mathrm{sides})$")
    fig.legend(handles=handles,loc='lower center',ncol=5,frameon=False,fontsize=10.5,bbox_to_anchor=(0.5,-0.01))
    plt.tight_layout(rect=[0,0.06,1,1])
    plt.savefig(figpath('family_panels.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()

    # ---- labeled Marion ----
    fig,ax=plt.subplots(figsize=(6.3,5.9)); fig.patch.set_facecolor('white')
    draw(ax,F(2),r"Marion arrangement $t=2$:  19 cells, type $12/3/3/1$")
    fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,fontsize=10,bbox_to_anchor=(0.5,-0.03))
    plt.tight_layout(rect=[0,0.09,1,1])
    plt.savefig(figpath('marion_arrangement.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()
    print("arrangement figs done")


def fig_tree_and_areas():
    plt.rcParams.update({'font.family':'DejaVu Sans'})

    # ---------------- clean dendrogram ----------------
    t=F(2); order=['A1','A2','B1','B2','C1','C2']
    nodes={}; cur=arrangement([],t); root=cur[0]
    nodes[0]={'parent':None,'stage':0,'ng':root[0],'ratio':root[1],'poly':root[2]}
    state=[(0,root[2])]; nid=1
    for si,key in enumerate(order,start=1):
        new=arrangement(order[:si],t); groups={}
        for ng,ratio,f in new:
            rp=f.representative_point()
            for pid,pp in state:
                if pp.contains(rp): groups.setdefault(pid,[]).append((ng,ratio,f));break
        ns=[]
        for pid,pp in state:
            ch=groups.get(pid,[])
            if len(ch)<=1: ns.append((pid,ch[0][2] if ch else pp))
            else:
                for ng,ratio,f in ch:
                    nodes[nid]={'parent':pid,'stage':si,'ng':ng,'ratio':ratio,'poly':f}; ns.append((nid,f)); nid+=1
        state=ns
    leaves=[p for p,_ in state]
    for l in leaves: nodes[l]['leaf']=True
    children={}
    for i,nd in nodes.items():
        p=nd['parent']
        if p is not None: children.setdefault(p,[]).append(i)
    def lspan(i):
        if i in children: return np.mean([lspan(c) for c in children[i]])
        return nodes[i]['poly'].representative_point().x
    for p in children: children[p].sort(key=lspan)
    xpos={};cnt=[0]
    def dfs(i):
        if i not in children or not children[i]: xpos[i]=cnt[0];cnt[0]+=1;return xpos[i]
        xs=[dfs(c) for c in children[i]];xpos[i]=float(np.mean(xs));return xpos[i]
    dfs(0)
    def key_of(nd): return classify(nd['ng'],nd['ratio'],t)
    NGN={3:'3',4:'4',5:'5',6:'6'}
    MRK={3:'^',4:'s',5:'p',6:'h'}          # actual polygon glyphs
    MSZ_I={3:175,4:140,5:160,6:185}        # internal-node sizes per shape
    MSZ_L={3:430,4:340,5:380,6:440}        # leaf sizes per shape
    TMK={'cornerT':'^','innerT':'^','quad':'s','pent':'p','hex':'h'}
    Ybase=0.0
    def yof(stage): return 6-stage   # root top (6), stage grows downward
    fig,ax=plt.subplots(figsize=(13,7.2)); fig.patch.set_facecolor('white')
    # edges (elbow)
    for i,nd in nodes.items():
        p=nd['parent']
        if p is not None:
            x0,y0=xpos[p],yof(nodes[p]['stage']); x1,y1=xpos[i],yof(nd['stage'])
            ax.plot([x0,x0,x1,x1],[y0,(y0+y1)/2,(y0+y1)/2,y1],color='#c7ccd2',lw=1.3,zorder=1,solid_capstyle='round')
    # leaf stems down to baseline
    for l in leaves:
        x=xpos[l];y=yof(nodes[l]['stage'])
        if y>Ybase: ax.plot([x,x],[y,Ybase],color='#e2e5e9',lw=1.1,ls=(0,(2,2)),zorder=1)
    # internal nodes
    for i,nd in nodes.items():
        if nd.get('leaf'): continue
        x=xpos[i];y=yof(nd['stage'])
        ng=nd['ng']
        ax.scatter([x],[y],marker=MRK[ng],s=MSZ_I[ng],facecolor='white',edgecolor='#9aa2ab',lw=1.2,zorder=3)
    # leaves at baseline
    for l in leaves:
        k=key_of(nodes[l]);x=xpos[l]
        ng=nodes[l]['ng']
        ax.scatter([x],[Ybase],marker=MRK[ng],s=MSZ_L[ng],facecolor=COL[k],edgecolor='#2b2b2b',lw=1.0,zorder=3)
    # stage labels + split-rule annotations
    rules={1:'+A₁ : 1×(△→△+△)',2:'+A₂ : 1×(△→△+△)',3:'+B₁ : 1×(△→△+△), 2×(△→△+□)',
           4:'+B₂ : 1×(△→△+△), 2×(△→△+□)',5:'+C₁ : 3×(△→△+□), 2×(□→△+⬠)',
           6:'+C₂ : 2×(△→△+□), 2×(□→△+⬠), 1×(⬠→△+⬡)'}
    ax.text(-0.5,yof(0),'start  △',ha='right',va='center',fontsize=11,color='#333')
    for s in range(1,7):
        ax.text(cnt[0]-0.5+0.5, yof(s), '', fontsize=9)
        ax.annotate(rules[s],xy=(cnt[0]-1+0.2,yof(s)),xytext=(cnt[0]+0.2,yof(s)),
                    fontsize=9.2,color='#444',va='center',ha='left')
        ax.axhline(yof(s),color='#f4f5f6',lw=0.6,zorder=0)
    handles=[Line2D([0],[0],marker=TMK[k],color='w',markerfacecolor=COL[k],markeredgecolor='#2b2b2b',markersize=12,label=LAB[k]) for k in ['cornerT','innerT','quad','pent','hex']]
    ax.legend(handles=handles,loc='lower center',bbox_to_anchor=(0.42,-0.11),ncol=5,frameon=False,fontsize=10)
    ax.set_title("Cell genealogy of the Marion arrangement ($t=2$)\nEach row adds one cevian; every cut shaves a triangle. Node shape = polygon type (triangle / quadrilateral / pentagon / hexagon).",fontsize=12.5,pad=12)
    ax.set_xlim(-2.2,cnt[0]+7.0);ax.set_ylim(-0.7,6.6);ax.axis('off')
    plt.tight_layout()
    plt.savefig(figpath('genealogy_tree.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()
    print("tree v2 done")

    # ---------------- area plot reparametrized s=t/(t+2) ----------------
    def formulas(tt):
        return (2*(tt-1)**2/((tt+2)*(2*tt+1)),
                tt*(tt-1)**2/((tt+2)*(2*tt+1)*(tt**2+tt+1)),
                1/((tt+1)*(tt**2+tt+1)),
                2*(tt-1)/((tt+2)*(tt**2+tt+1)),
                (tt-1)*(tt**2+3*tt+1)/((tt+1)*(2*tt+1)*(tt**2+tt+1)))
    def true_area(tv):
        if tv==0: return (1.0,0,0,0,0)
        s=tv if tv>=1 else 1/tv
        return formulas(s)
    # axis coordinate s=t/(t+2): t=0->0, t=2->1/2, t=inf->1
    def S(tv): return tv/(tv+2)
    tgrid=np.concatenate([np.linspace(0,2,400),np.linspace(2,60,400)])
    sg=np.array([S(tv) for tv in tgrid])
    vals=np.array([true_area(tv) for tv in tgrid])
    fig,ax=plt.subplots(figsize=(9.4,5.4)); fig.patch.set_facecolor('white')
    names=['hexagon (×1)','inner triangle (×6)','corner triangle (×6)','quadrilateral (×3)','pentagon (×3)']
    cols=[COL['hex'],COL['innerT'],COL['cornerT'],COL['quad'],COL['pent']]
    orderplot=[0,4,3,2,1]
    for i in orderplot:
        ax.plot(sg,vals[:,i],color=cols[i],lw=2.3,label=names[i])
    # landmark ticks equally spaced: t=0,1,2,infty -> s=0,1/3,1/2,1 ; also add t=4 (s=2/3)
    ticks=[(0,'0'),(1,'1'),(2,'2'),(4,'4'),(10,'10')]
    xt=[S(tv) for tv,_ in ticks]+[1.0]
    xl=[lab for _,lab in ticks]+['∞']
    ax.set_xticks(xt); ax.set_xticklabels(xl)
    ax.axvline(S(1),color='#888',ls='--',lw=1.0); ax.text(S(1)+0.008,0.86,'medians',color='#888',fontsize=9,rotation=90,va='top')
    ax.axvline(S(2),color='#888',ls=':',lw=1.0); ax.text(S(2)+0.008,0.93,'Marion',color='#666',fontsize=9)
    ax.axvline(0,color='#bbb',lw=0.8); ax.axvline(1,color='#bbb',lw=0.8)
    ax.text(0.004,1.0,'sides',fontsize=8.5,color='#666',va='center')
    ax.text(0.995,1.0,'sides',fontsize=8.5,color='#666',va='center',ha='right')
    ax.set_xlabel('parameter  $t$   (axis $s=t/(t+2)$: the degenerations $t=0,\\infty$ and Marion $t=2$ equally spaced)',fontsize=10.5)
    ax.set_ylabel('true area ratio (cell / triangle)',fontsize=11.5)
    ax.set_title("Geometric cell areas across the whole range $t\\in[0,\\infty]$",fontsize=12.5)
    ax.legend(frameon=False,fontsize=10,loc='center left')
    ax.set_xlim(-0.01,1.01);ax.set_ylim(0,1.03);ax.grid(alpha=0.22)
    plt.tight_layout()
    plt.savefig(figpath('cell_areas.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()
    print("area v2 done")


def fig_orientation_and_signed():
    plt.rcParams.update({'font.family':'DejaVu Sans'})

    def ccw(poly):
        pts=list(poly.exterior.coords)[:-1]
        a=0
        for i in range(len(pts)):
            x1,y1=pts[i];x2,y2=pts[(i+1)%len(pts)];a+=x1*y2-x2*y1
        if a<0: pts=pts[::-1]
        return pts

    # ============ orientation figure ============
    fig,ax=plt.subplots(figsize=(7.0,6.5)); fig.patch.set_facecolor('white')
    tr=F(2)
    for ng,ratio,f in arrangement(['A1','A2','B1','B2','C1','C2'],tr):
        key=classify(ng,ratio,tr)
        pts=[to_eq(c) for c in list(f.exterior.coords)[:-1]]
        ax.add_patch(MplPoly(pts,closed=True,facecolor=COL[key],edgecolor='white',lw=1.1,alpha=0.9))
        # CCW arrows: map polygon vertices to eq, ensure ccw, draw arrowheads on 2 longest edges
        P=[to_eq(c) for c in ccw(f)]
        n=len(P)
        edges=sorted(range(n),key=lambda i:-(np.hypot(P[(i+1)%n][0]-P[i][0],P[(i+1)%n][1]-P[i][1])))[:2]
        cen=np.mean(P,axis=0)
        for i in edges:
            p=np.array(P[i]);q=np.array(P[(i+1)%n]); m=(p+q)/2
            d=(q-p);d=d/ (np.hypot(*d)+1e-9)
            m2=m+0.18*(m-cen)  # pull slightly outward? keep near edge
            a=FancyArrowPatch(tuple(m-0.055*d),tuple(m+0.055*d),arrowstyle='->',mutation_scale=6,color='#222',lw=0.8,zorder=5)
            ax.add_patch(a)
    ax.add_patch(MplPoly([EA,EB,EC],closed=True,fill=False,edgecolor='#222',lw=1.8))
    ax.set_aspect('equal');ax.axis('off')
    ax.set_title("Orientation (Umlaufrichtung) of the cells at $t=2$\nAll boundaries run counter-clockwise $\\Rightarrow$ all signed areas $+$, summing to $+1$",fontsize=12,pad=10)
    plt.tight_layout();plt.savefig(figpath('orientation.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()
    print("orient done")

    # ============ signed cancellation figure ============
    tt=sp.symbols('t')
    forms={'hexagon':2*(tt-1)**2/((tt+2)*(2*tt+1)),
           'inner tri.':tt*(tt-1)**2/((tt+2)*(2*tt+1)*(tt**2+tt+1)),
           'corner tri.':1/((tt+1)*(tt**2+tt+1)),
           'quadrilateral':2*(tt-1)/((tt+2)*(tt**2+tt+1)),
           'pentagon':(tt-1)*(tt**2+3*tt+1)/((tt+1)*(2*tt+1)*(tt**2+tt+1))}
    mult={'hexagon':1,'inner tri.':6,'corner tri.':6,'quadrilateral':3,'pentagon':3}
    cmap={'hexagon':COL['hex'],'inner tri.':COL['innerT'],'corner tri.':COL['cornerT'],'quadrilateral':COL['quad'],'pentagon':COL['pent']}
    labels=list(forms.keys())
    def totals(tv): return [float(mult[k]*forms[k].subs(tt,tv)) for k in labels]
    T2=totals(sp.Integer(2)); Th=totals(sp.Rational(1,2))

    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.2),gridspec_kw={'width_ratios':[1.35,1]}); fig.patch.set_facecolor('white')
    x=np.arange(len(labels)); w=0.38
    b1=ax1.bar(x-w/2,T2,w,color=[cmap[k] for k in labels],edgecolor='#333',lw=0.6,label='$t=2$')
    b2=ax1.bar(x+w/2,Th,w,color=[cmap[k] for k in labels],edgecolor='#333',lw=0.6,hatch='//',alpha=0.9,label='$t=1/2$ (mirror)')
    ax1.axhline(0,color='#333',lw=1.0)
    ax1.axhline(1,color='#555',lw=1.2,ls='--'); ax1.text(len(labels)-0.5,1.02,'signed sum $=1$',color='#555',fontsize=9,ha='right')
    ax1.set_xticks(x); ax1.set_xticklabels(labels,rotation=18,ha='right',fontsize=9.5)
    ax1.set_ylabel('signed group total (× multiplicity)',fontsize=10.5)
    ax1.set_title("Signed areas: group totals at $t=2$ and its mirror $t=1/2$",fontsize=11.5)
    for xi,v in zip(x-w/2,T2): ax1.text(xi,v+0.03,f'{v:.2f}',ha='center',fontsize=7.5,color='#333')
    for xi,v in zip(x+w/2,Th): ax1.text(xi,v+ (0.03 if v>=0 else -0.10),f'{v:+.2f}',ha='center',fontsize=7.5,color='#333')
    ax1.legend(frameon=False,fontsize=10,loc='upper left')
    ax1.grid(axis='y',alpha=0.2); ax1.set_ylim(-1.05,2.55)

    # right: corner construction small (t=2) vs large (t=1/2)
    def bary_to_eq(P):
        x_,y_,z_=[float(c) for c in P]; s=x_+y_+z_; X=y_/s; Y=z_/s
        return to_eq((X,Y))
    def corner_tri(tv):
        A=(1,0,0); F2=(tv/(tv+1),1/(tv+1),0); F3=(tv**2,tv,1)
        return [bary_to_eq(A),bary_to_eq(F2),bary_to_eq(F3)]
    ax2.add_patch(MplPoly([EA,EB,EC],closed=True,fill=False,edgecolor='#222',lw=1.8))
    # t=2 small (positive)
    c2=corner_tri(2.0)
    ax2.add_patch(MplPoly(c2,closed=True,facecolor=COL['cornerT'],edgecolor='#204060',lw=1.4,alpha=0.95))
    # t=1/2 large (positive but bigger)
    ch=corner_tri(0.5)
    ax2.add_patch(MplPoly(ch,closed=True,facecolor=COL['cornerT'],edgecolor='#204060',lw=1.4,alpha=0.35,hatch='//'))
    # arrows CCW on both
    for poly,txt,off in [(c2,'$t=2$:  $1/21$',(0.02,0.02)),(ch,'$t=1/2$:  $8/21$ (same labels, outer shape)',(0,0))]:
        P=poly; 
        # ensure ccw
        a=0
        for i in range(len(P)):
            a+=P[i][0]*P[(i+1)%len(P)][1]-P[(i+1)%len(P)][0]*P[i][1]
        if a<0: P=P[::-1]
        cen=np.mean(P,axis=0)
        for i in range(len(P)):
            p=np.array(P[i]);q=np.array(P[(i+1)%len(P)]);m=(p+q)/2;d=q-p;d=d/(np.hypot(*d)+1e-9)
            ax2.add_patch(FancyArrowPatch(tuple(m-0.05*d),tuple(m+0.05*d),arrowstyle='->',mutation_scale=6,color='#111',lw=0.9,zorder=6))
    ax2.text(0.5,-0.18,'corner-triangle construction\n(vertices $A,\\ AB\\cap L_5,\\ L_5\\cap L_1$)',ha='center',fontsize=9,transform=ax2.transAxes)
    ax2.text(0.5,1.02,'small at $t=2$  →  large ("outer") at $t=1/2$',ha='center',fontsize=10,transform=ax2.transAxes,color='#204060')
    ax2.set_aspect('equal');ax2.axis('off')
    plt.tight_layout();plt.savefig(figpath('signed_cancellation.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()
    print("signed done")


def fig_involution():
    plt.rcParams.update({'font.family':'DejaVu Sans'})

    A=np.array(EA);B=np.array(EB);C=np.array(EC)
    C1COL='#D9534F'  # first cevians (red)
    C2COL='#5B8FD6'  # second cevians (blue)

    def feet(t):
        f1=1/(t+1); f2=t/(t+1)
        # A->BC, B->CA, C->AB ; first=f1, second=f2
        first=[(A,B+f1*(C-B)),(B,C+f1*(A-C)),(C,A+f1*(B-A))]
        second=[(A,B+f2*(C-B)),(B,C+f2*(A-C)),(C,A+f2*(B-A))]
        return first,second

    def draw_config(ax,t,title,swapcolors=False):
        ax.add_patch(MplPoly([tuple(A),tuple(B),tuple(C)],closed=True,facecolor='#f7f7f4',edgecolor='#222',lw=1.8))
        first,second=feet(t)
        ca,cb=(C2COL,C1COL) if swapcolors else (C1COL,C2COL)
        for p,q in first:
            ax.plot([p[0],q[0]],[p[1],q[1]],color=ca,lw=2.4,solid_capstyle='round')
        for p,q in second:
            ax.plot([p[0],q[0]],[p[1],q[1]],color=cb,lw=2.4,solid_capstyle='round',ls=(0,(1,1)))
        ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=12,pad=8)

    def draw_medians(ax,title):
        ax.add_patch(MplPoly([tuple(A),tuple(B),tuple(C)],closed=True,facecolor='#f7f7f4',edgecolor='#222',lw=1.8))
        Mab=(A+B)/2;Mbc=(B+C)/2;Mca=(C+A)/2
        for p,q in [(A,Mbc),(B,Mca),(C,Mab)]:
            ax.plot([p[0],q[0]],[p[1],q[1]],color='#888',lw=2.6,solid_capstyle='round')
        ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=12,pad=8)

    # ---------- involution figure ----------
    fig,axs=plt.subplots(1,3,figsize=(13,4.4)); fig.patch.set_facecolor('white')
    draw_config(axs[0],F(3),r"$t=3$",swapcolors=False)
    draw_medians(axs[1],r"$t=1$ (fixed point): feet merge $\to$ medians")
    draw_config(axs[2],F(1,3),r"$t=1/3=1/t$: same lines, roles swapped",swapcolors=True)
    leg=[Line2D([0],[0],color=C1COL,lw=2.6,label='"first" cevians (feet at $\\frac{1}{t+1}$)'),
         Line2D([0],[0],color=C2COL,lw=2.6,ls=(0,(1,1)),label='"second" cevians (feet at $\\frac{t}{t+1}$)')]
    fig.legend(handles=leg,loc='lower center',ncol=2,frameon=False,fontsize=10.5,bbox_to_anchor=(0.5,-0.02))
    fig.suptitle(r"The configuration under $\iota:\,t\mapsto 1/t$ — the six lines are unchanged, only the pairing of cevians is swapped",
                 fontsize=12.5,y=1.02)
    plt.tight_layout(rect=[0,0.05,1,0.98])
    plt.savefig(figpath('involution_cevians.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()

    # ---------- invariant vs non-invariant cells ----------
    INVAR={'hex','innerT'}   # B(u)=0
    fig,ax=plt.subplots(figsize=(6.6,6.3)); fig.patch.set_facecolor('white')
    tr=F(2)
    for ng,ratio,f in arrangement(['A1','A2','B1','B2','C1','C2'],tr):
        key=classify(ng,ratio,tr)
        pts=[to_eq(c) for c in list(f.exterior.coords)[:-1]]
        if key in INVAR:
            ax.add_patch(MplPoly(pts,closed=True,facecolor=COL[key],edgecolor='white',lw=1.2))
        else:
            ax.add_patch(MplPoly(pts,closed=True,facecolor=COL[key],edgecolor='white',lw=1.2,alpha=0.45,hatch='///'))
    ax.add_patch(MplPoly([EA,EB,EC],closed=True,fill=False,edgecolor='#222',lw=1.8))
    ax.set_aspect('equal');ax.axis('off')
    ax.set_title(r"Invariant vs. non-invariant cells under $\iota:t\mapsto 1/t$  ($t=2$)",fontsize=12.5,pad=10)
    leg=[MplPoly([(0,0)],facecolor=COL['hex'],edgecolor='white',label='hexagon — invariant ($B=0$)'),
         MplPoly([(0,0)],facecolor=COL['innerT'],edgecolor='white',label='inner triangle — invariant ($B=0$)'),
         MplPoly([(0,0)],facecolor=COL['cornerT'],edgecolor='white',alpha=0.45,hatch='///',label='corner triangle — sign part ($B\\neq0$)'),
         MplPoly([(0,0)],facecolor=COL['quad'],edgecolor='white',alpha=0.45,hatch='///',label='quadrilateral — sign part ($B\\neq0$)'),
         MplPoly([(0,0)],facecolor=COL['pent'],edgecolor='white',alpha=0.45,hatch='///',label='pentagon — sign part ($B\\neq0$)')]
    fig.legend(handles=leg,loc='lower center',ncol=2,frameon=False,fontsize=9.5,bbox_to_anchor=(0.5,-0.05))
    plt.tight_layout(rect=[0,0.10,1,1])
    plt.savefig(figpath('involution_invariance.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()
    print("involution + invariance figs done")


def fig_involution_by_type(reps=None):
    if reps is None:
        reps=json.load(open(datapath('representatives.json')))

    def labelline(name,t):
        L=cevian_lines(t)
        D={**L,'AB':SIDELINES[0],'CA':SIDELINES[1],'BC':SIDELINES[2]}
        return D[name]
    def vert(pair,t): return inter(labelline(pair[0],t),labelline(pair[1],t))
    def poly_eq(desc,t): return [to_eq(vert(tuple(p),t)) for p in desc]
    def sarea_XY(desc,t):
        V=[vert(tuple(p),t) for p in desc]; s=F(0); n=len(V)
        for i in range(n):
            x1,y1=V[i];x2,y2=V[(i+1)%n];s+=x1*y2-x2*y1
        return s  # = signed area ratio

    order=['hex','innerT','cornerT','quad','pent']
    title={'hex':'hexagon','innerT':'inner triangle','cornerT':'corner triangle','quad':'quadrilateral','pent':'pentagon'}
    NEG='#8E2B22'; POS='#111111'

    fig,axes=plt.subplots(2,5,figsize=(15,6.4)); fig.patch.set_facecolor('white')
    for col,key in enumerate(order):
        desc=reps[key]
        for row,t in enumerate([F(2),F(1,2)]):
            ax=axes[row][col]
            ax.add_patch(MplPoly([EA,EB,EC],closed=True,fill=False,edgecolor='#333',lw=1.5))
            P=poly_eq(desc,t); sa=sarea_XY(desc,t)
            neg = sa<0
            fc=COL[key]
            ax.add_patch(MplPoly(P,closed=True,facecolor=fc,edgecolor='white',lw=1.2,
                                 alpha=0.55 if neg else 0.95, hatch='///' if neg else None))
            # orientation arrows following the stored vertex order
            Pa=np.array(P); cen=Pa.mean(axis=0); n=len(P)
            for i in range(n):
                p=Pa[i];q=Pa[(i+1)%n];m=(p+q)/2;d=q-p;L=np.hypot(*d)
                if L<1e-6: continue
                d=d/L
                ax.add_patch(FancyArrowPatch(tuple(m-0.045*d),tuple(m+0.045*d),arrowstyle='->',
                             mutation_scale=6,color=POS,lw=0.8,zorder=6))
            val=F(sa).limit_denominator(10**6)
            sign='' if val>=0 else '−'
            ax.text(0.5,-0.06,f"{'−' if neg else '+'}  area $= {'-' if val<0 else ''}\\frac{{{abs(val.numerator)}}}{{{abs(val.denominator)}}}$",
                    transform=ax.transAxes,ha='center',fontsize=11,color=NEG if neg else '#111')
            ax.set_aspect('equal'); ax.axis('off')
            ax.set_xlim(-0.08,1.08); ax.set_ylim(-0.18,0.98)
        axes[0][col].set_title(title[key],fontsize=12,pad=4,color=COL[key],fontweight='bold')
    # row labels
    axes[0][0].text(-0.13,0.5,'$t=2$',transform=axes[0][0].transAxes,rotation=90,va='center',ha='center',fontsize=13,fontweight='bold')
    axes[1][0].text(-0.13,0.5,'$t=\\frac{1}{2}$',transform=axes[1][0].transAxes,rotation=90,va='center',ha='center',fontsize=13,fontweight='bold')
    fig.suptitle("How each cell type moves under the involution $t\\mapsto 1/t$  (same label-construction, evaluated at $t=2$ and its mirror $t=\\frac{1}{2}$)",
                 fontsize=13,y=0.99)
    fig.text(0.5,0.005,"Hexagon and inner triangle are unchanged (invariant).  The corner triangle stays positively oriented but grows into a large \"outer\" triangle "
             "($\\frac{1}{21}\\!\\to\\!\\frac{8}{21}$).\nThe quadrilateral and pentagon reverse orientation (red, clockwise arrows) and become negative — the overcount of the corner triangles is cancelled exactly, signed total $=1$.",
             ha='center',fontsize=9.5,color='#444')
    plt.tight_layout(rect=[0.02,0.055,1,0.96])
    plt.savefig(figpath('involution_by_type.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()
    print("movement fig done")


def main():
    reps = extract_representatives()        # -> representatives.json
    fig_family()                            # family_panels, marion_arrangement
    fig_tree_and_areas()                    # genealogy_tree, cell_areas
    fig_orientation_and_signed()            # orientation, signed_cancellation
    fig_involution()                        # involution_cevians (+ involution_invariance)
    fig_involution_by_type(reps)            # involution_by_type
    print("all figures written to", FIGDIR)


if __name__ == "__main__":
    main()
