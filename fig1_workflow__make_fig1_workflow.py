"""Figure 1 - benchmark workflow diagram (pure matplotlib drawing, no data needed).

Output: outputs/figures/fig1_workflow.png
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
fig,ax=plt.subplots(figsize=(15,8.2)); ax.set_xlim(0,150); ax.set_ylim(0,82); ax.axis("off")
DARK="#2f5d9e"; LIGHT="#d6e6c3"; PANEL={"pre":"#e6eff9","mod":"#f1f1f1","eval":"#fbf3e4","seq":"#efe8f6"}
def panel(x,y,w,h,title,key):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.3,rounding_size=1.2",fc=PANEL[key],ec="#555",lw=1.1,ls=(0,(4,3))))
    ax.text(x+w/2,y+h-2.6,title,ha="center",va="center",fontsize=12.5,weight="bold",color="#222")
def box(x,y,w,h,head,body=None,hc=DARK,bc=LIGHT,fs=9.6):
    if body is None:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.15,rounding_size=0.6",fc=hc,ec="none"))
        ax.text(x+w/2,y+h/2,head,ha="center",va="center",color="white",fontsize=fs,weight="bold",wrap=True); return
    hh=h*0.42
    ax.add_patch(FancyBboxPatch((x,y+h-hh),w,hh,boxstyle="round,pad=0.15,rounding_size=0.6",fc=hc,ec="none"))
    ax.text(x+w/2,y+h-hh/2,head,ha="center",va="center",color="white",fontsize=fs,weight="bold")
    ax.add_patch(FancyBboxPatch((x,y),w,h-hh-0.4,boxstyle="round,pad=0.15,rounding_size=0.6",fc=bc,ec="none"))
    ax.text(x+w/2,y+(h-hh-0.4)/2,body,ha="center",va="center",color="#222",fontsize=fs-1)
def arrow(x1,y1,x2,y2,c="#7fb069",lw=2.2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=16,color=c,lw=lw))
# Panel 1: data preprocessing
panel(1,1,38,79,"Data preprocessing","pre")
steps=[("Positive sites","Published hg38 site records\n(RMBase; +1,815 Nmix sites for Nm)"),
       ("Map to annotation","hg38 genome, GENCODE v47 genes / transcripts"),
       ("Redundancy removal","Exact duplicates; CD-HIT-EST or\nfixed-position identity filters (60–90%)"),
       ("Presumed negatives","Same gene / transcript and centre base,\nknown sites excluded, fixed seeds"),
       ("41-nt windows, 1:1 balanced","46,860 sequences across 6 modifications"),
       ("Gene-disjoint split","train / validation / held-out test")]
y=66.5
for i,(h,b) in enumerate(steps):
    box(4,y,32,9.2,h,b)
    if i<len(steps)-1: arrow(20,y-0.3,20,y-2.4)
    y-=12.2
# Panel 2: model construction
panel(42,1,40,79,"Model construction","mod")
mods=[("m6A","20,000"),("Nm","7,112"),("m5C","10,000"),("Ψ","4,656"),("m1A","1,810"),("m7G","1,482")]
for i,(m,n) in enumerate(mods):
    cx=45+(i%2)*18; cy=60-(i//2)*15.5
    box(cx,cy,16,12,f"{m} dataset",f"{n} sequences\n5 published methods")
ax.text(62,17,"30 local retrainings at 41 nt",ha="center",fontsize=11,weight="bold",color="#222")
ax.text(62,11.5,"Fit on training records only\nCheckpoint / threshold chosen on validation\nTest records never used for fitting",ha="center",va="center",fontsize=9.3,color="#333",linespacing=1.5)
# Panel 3: evaluation
panel(85,41,64,39,"Performance evaluation","eval")
ev=[("Held-out test","Acc, Pre, Sn, Sp, F1,\nMCC, AUROC, AUPRC\n+ 95% bootstrap CI"),
    ("Generalization gap","Full-data vs test\nscores; training −\ntest AUROC"),
    ("Error structure","Confusion matrices,\nSn − Sp bias,\nROC curves")]
for i,(h,b) in enumerate(ev): box(88+i*20.3,47,18.6,24,h,b,hc="#b5651d",bc="#f6dfc0",fs=9.4)
# Panel 4: sequence analysis
panel(85,1,64,37,"Sequence-pattern analysis","seq")
sq=[("Two-sample logos","Position-wise enrichment,\nBonferroni p < 0.05"),
    ("ENAC features","Window 2, ranked by\ninformation gain (top 30)"),
    ("GC composition","Flank GC content,\npositives vs negatives")]
for i,(h,b) in enumerate(sq): box(88+i*20.3,6,18.6,24,h,b,hc="#6a3d9a",bc="#e3d7f0",fs=9.4)
# arrows between panels
arrow(36.5,8.4,44.5,8.4,c="#444")     # split -> models
arrow(82.6,62,87.2,62,c="#444")        # models -> evaluation
arrow(82.6,20,87.2,20,c="#444")        # datasets -> sequence analysis
out = config.OUT_FIG / "fig1_workflow.png"
fig.savefig(out,dpi=170,bbox_inches="tight",facecolor="white")
print("written", out)
