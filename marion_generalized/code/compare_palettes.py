import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly, Rectangle
from fractions import Fraction as F
from common import arrangement, classify, to_eq, EA, EB, EC, figpath

# three refined palettes: keys cornerT, innerT, quad, pent, hex (+accent)
PALETTES={
 'A — Jewel tones':      {'hex':'#C0392B','pent':'#0E8A6E','quad':'#C9992B','cornerT':'#2E5E8C','innerT':'#6C4C93','acc':'#188E8E'},
 'B — Muted editorial':  {'hex':'#B0413E','pent':'#5E8A55','quad':'#D0A24A','cornerT':'#41708E','innerT':'#836593','acc':'#4E8D7C'},
 'C — Vivid refined':    {'hex':'#E24A33','pent':'#2FB37A','quad':'#EDA419','cornerT':'#2E7FC1','innerT':'#8E5EC9','acc':'#1BB5B0'},
}
tval=F(2); faces=arrangement(['A1','A2','B1','B2','C1','C2'],tval)
typ_order=['hex','pent','quad','cornerT','innerT','acc']
typ_name={'hex':'hexagon','pent':'pentagon','quad':'quadrilateral','cornerT':'corner tri.','innerT':'inner tri.','acc':'accent'}

fig=plt.figure(figsize=(14.5,6.2)); fig.patch.set_facecolor('white')
for j,(name,pal) in enumerate(PALETTES.items()):
    ax=fig.add_axes([0.02+j*0.325, 0.28, 0.30, 0.62])
    for ng,ratio,f in faces:
        key=classify(ng,ratio,tval)
        pts=[to_eq(c) for c in list(f.exterior.coords)[:-1]]
        ax.add_patch(MplPoly(pts,closed=True,facecolor=pal[key],edgecolor='white',lw=1.2))
    ax.add_patch(MplPoly([EA,EB,EC],closed=True,fill=False,edgecolor='#2b2b2b',lw=1.8))
    ax.set_aspect('equal');ax.axis('off')
    ax.set_xlim(-0.05,1.05);ax.set_ylim(-0.03,0.95)
    ax.set_title(name,fontsize=13,pad=8,color='#222',fontweight='bold')
    # swatch strip
    axs=fig.add_axes([0.02+j*0.325, 0.06, 0.30, 0.17]); axs.axis('off')
    axs.set_xlim(0,6); axs.set_ylim(0,2)
    for i,k in enumerate(typ_order):
        axs.add_patch(Rectangle((i+0.08,0.75),0.84,0.9,facecolor=pal[k],edgecolor='#ccc',lw=0.6))
        axs.text(i+0.5,0.55,typ_name[k],ha='center',va='top',fontsize=6.6,color='#333')
        axs.text(i+0.5,1.2,pal[k],ha='center',va='center',fontsize=5.6,color='white',fontweight='bold')
fig.suptitle("Three refined palettes for the cell types — pick one (same six hues, more sophisticated tuning)",
             fontsize=13,y=0.99)
plt.savefig(figpath('palette_comparison.png'),dpi=170,bbox_inches='tight',facecolor='white');plt.close()
print("palette comparison done")
