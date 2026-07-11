import sys, warnings; warnings.filterwarnings("ignore"); sys.path.insert(0,"scripts")
import pandas as pd, numpy as np, re, cv2, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import RidgeCV
from embedding_io import read_embedding_table, image_key
D=2016; OUT=240; DST=np.array([[0,0],[D-1,0],[D-1,D-1],[0,D-1]],dtype=np.float32)
BLUE="#2e6e8e"; RED="#cc4b27"

ic=read_embedding_table("data/generatable/dimreduction_horn90/ic_scores.csv")
iccols=sorted([c for c in ic.columns if re.fullmatch(r"IC\d+",c)],key=lambda c:int(c[2:]))
ic["ik"]=ic.source_image_path.map(image_key)
HS=pd.read_csv("data/provided/human_disease_scores.csv"); HS["ik"]=HS.image_id.map(image_key)
HS=HS.dropna(subset=["human_score"]).groupby("ik",as_index=False).agg(
     human_score=("human_score","mean"),environment=("environment","first"))
HS=HS[HS.environment=="Nebraska2025"]
m=ic.groupby("ik",as_index=False).agg({**{c:"mean" for c in iccols},"genotype":"first"}).merge(HS,on="ik")
y=m.human_score.to_numpy(); g=m.genotype.to_numpy()
absr=m[iccols].corrwith(m.human_score,method="spearman").abs().sort_values(ascending=False)
above=absr.index.tolist()[:19]
A=np.logspace(-2,4,13); genos=pd.unique(g)
rng=np.random.RandomState(0); fo=dict(zip(genos,rng.randint(0,5,len(genos)))); folds=np.array([fo[x] for x in g])
def oof(cols):
    X=m[cols].to_numpy(float); o=np.full(len(y),np.nan)
    for f in range(5):
        te=folds==f; tr=~te; o[te]=RidgeCV(alphas=A).fit(X[tr],y[tr]).predict(X[te])
    return o
predA,predB=oof(above),oof(above+["IC66"])
r2=lambda o:1-np.nansum((y-o)**2)/np.nansum((y-y.mean())**2); r2A,r2B=r2(predA),r2(predB)
prr,spr=pearsonr(m.IC66,y)[0],spearmanr(m.IC66,y)[0]

# montage selected from the SAME Nebraska-scored images, by image-mean IC66
low_imgs=m.nsmallest(16,"IC66"); high_imgs=m.nlargest(16,"IC66")
low_iks=set(low_imgs.ik); high_iks=set(high_imgs.ik)
def crop_row(ik,low):
    sub=ic[ic.ik==ik]; return (sub.nsmallest(1,"IC66") if low else sub.nlargest(1,"IC66")).iloc[0]
lows=[crop_row(k,True)  for k in low_imgs.sort_values("IC66").ik]
highs=[crop_row(k,False) for k in high_imgs.sort_values("IC66",ascending=False).ik]
def reconstruct(row):
    src=cv2.imread(str(row["source_image_path"]))
    if src is None: return None
    c=np.array([[row[f"crop_corner_{i}_x"],row[f"crop_corner_{i}_y"]] for i in range(4)],dtype=np.float32)
    return cv2.cvtColor(cv2.resize(cv2.warpPerspective(src,cv2.getPerspectiveTransform(c,DST),(D,D)),(OUT,OUT)),cv2.COLOR_BGR2RGB)

plt.rcParams.update({"font.size":11})
fig=plt.figure(figsize=(6.5, 3.25))
outer=fig.add_gridspec(1,2,width_ratios=[2.95,1.0],wspace=0.10)
gm=outer[0,0].subgridspec(4,9,width_ratios=[1,1,1,1,0.3,1,1,1,1],wspace=0.04,hspace=0.04)
lax=[]; hax=[]
for idx in range(16):
    r4,c4=divmod(idx,4)
    for axl,data,c0,col in [(lax,lows,0,BLUE),(hax,highs,5,RED)]:
        ax=fig.add_subplot(gm[r4,c0+c4]); im=reconstruct(data[idx])
        if im is not None: ax.imshow(im)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_color(col); sp.set_linewidth(1.4)
        axl.append(ax)
gr=outer[0,1].subgridspec(2,1,hspace=0.30); ax1=fig.add_subplot(gr[0]); ax2=fig.add_subplot(gr[1])
for a in (ax1,ax2): a.set_box_aspect(1)
oth=m[~m.ik.isin(low_iks|high_iks)]
ax1.scatter(oth.IC66,oth.human_score,s=10,alpha=.3,color="#c2c2c2",edgecolor="none")
ax1.scatter(m[m.ik.isin(low_iks)].IC66,m[m.ik.isin(low_iks)].human_score,s=30,color=BLUE,edgecolor="k",lw=.4,zorder=3)
ax1.scatter(m[m.ik.isin(high_iks)].IC66,m[m.ik.isin(high_iks)].human_score,s=30,color=RED,edgecolor="k",lw=.4,zorder=3)
ax1.set_xlabel("IC66 score (image mean)"); ax1.set_ylabel("human disease score")
ax1.set_title("IC66 vs human disease scores",fontsize=11,fontweight="bold")
ax1.text(.96,.04,f"Pearson r²={prr**2:.2f}\nSpearman ρ={spr:.2f}",transform=ax1.transAxes,
         va="bottom",ha="right",fontsize=9,bbox=dict(boxstyle="round",fc="white",ec="grey",alpha=.85))
ax1.grid(alpha=.25)
lo,hi=np.nanmin(y),np.nanmax(y); ax2.plot([lo,hi],[lo,hi],ls="--",color="grey",lw=1)
ax2.scatter(y,predA,s=12,alpha=.45,color="#6a89a8",edgecolor="none",label=f"Top 19 ICs (R²={r2A:.2f})")
ax2.scatter(y,predB,s=12,alpha=.45,color=RED,edgecolor="none",label=f"Top 19 + IC66 (R²={r2B:.2f})")
ax2.set_xlabel("observed human disease score"); ax2.set_ylabel("predicted (out-of-fold)")
ax2.set_title("Effect of adding IC66",fontsize=11,fontweight="bold")
ax2.legend(fontsize=9,loc="upper left"); ax2.grid(alpha=.25)
fig.canvas.draw()
def ext(axs): p=[a.get_position() for a in axs]; return min(q.x0 for q in p),max(q.x1 for q in p),max(q.y1 for q in p)
for axs,txt,col in [(lax,"lower IC66",BLUE),(hax,"higher IC66",RED)]:
    x0,x1,yt=ext(axs); fig.text((x0+x1)/2,yt+0.015,txt,ha="center",va="bottom",fontsize=12,fontweight="bold",color=col)
    fig.add_artist(Line2D([x0,x1],[yt+0.010]*2,color=col,lw=1.4,transform=fig.transFigure))
x0=min(a.get_position().x0 for a in lax); yt=max(a.get_position().y1 for a in lax)
fig.text(x0-0.01,yt+0.05,"A",fontsize=17,fontweight="bold",ha="left",va="top")
ax1.text(-0.20,1.18,"B",transform=ax1.transAxes,fontsize=17,fontweight="bold",va="top")
ax2.text(-0.20,1.18,"C",transform=ax2.transAxes,fontsize=17,fontweight="bold",va="top")
fig.savefig("figures/ica_characterization/IC66_detail.png",bbox_inches="tight",dpi=300)
print(f"blue dots={len(low_iks)} red dots={len(high_iks)}  r2A={r2A:.3f} r2B={r2B:.3f}")
print("wrote IC66_detail.png")
