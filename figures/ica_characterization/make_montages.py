import sys, warnings; warnings.filterwarnings("ignore"); sys.path.insert(0,"scripts")
import pandas as pd, numpy as np, cv2, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
D=2016; OUT=240
DST=np.array([[0,0],[D-1,0],[D-1,D-1],[0,D-1]],dtype=np.float32)
ic=pd.read_csv("data/generatable/dimreduction_horn90/ic_scores.csv")

def reconstruct(row):
    src=cv2.imread(str(row["source_image_path"]))
    if src is None: return None
    c=np.array([[row[f"crop_corner_{i}_x"],row[f"crop_corner_{i}_y"]] for i in range(4)],dtype=np.float32)
    crop=cv2.warpPerspective(src,cv2.getPerspectiveTransform(c,DST),(D,D))
    return cv2.cvtColor(cv2.resize(crop,(OUT,OUT)),cv2.COLOR_BGR2RGB)
def pick(df,col,n,low):
    s=df.sort_values(col,ascending=low); out=[]; seen=set()
    for _,r in s.iterrows():
        if r["source_image_path"] in seen: continue
        seen.add(r["source_image_path"]); out.append(r)
        if len(out)==n: break
    return out

def montage(df,rows,title,outpath,N,S=1.5, max_width_in=6.5):
    ncol=2*N+1; wr=[1]*N+[0.35]+[1]*N
    imgcols=list(range(N))+list(range(N+1,2*N+1))
    # Cap figure width at journal single-column-ish max (6.5 in @ 300 dpi).
    S = min(S, (max_width_in - 0.9) / (2 * N + 0.35))
    fig_w=(2*N+0.35)*S+0.9; fig_h=len(rows)*S+1.25
    fig,axs=plt.subplots(len(rows),ncol,figsize=(fig_w,fig_h),gridspec_kw={"width_ratios":wr})
    if len(rows)==1: axs=axs[None,:]
    for ri,(icn,lab) in enumerate(rows):
        sel=pick(df,icn,N,True)+pick(df,icn,N,False)
        for k,r in enumerate(sel):
            ax=axs[ri,imgcols[k]]; im=reconstruct(r)
            if im is not None: ax.imshow(im)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values(): sp.set_color("k"); sp.set_linewidth(1.1)
        axs[ri,0].set_ylabel(lab,rotation=0,ha="right",va="center",fontsize=8,labelpad=5)
        sp=axs[ri,N]; sp.set_xticks([]); sp.set_yticks([])
        for s in sp.spines.values(): s.set_visible(False)
    fig.tight_layout(rect=[0,0,1,0.85])
    fig.canvas.draw()
    xL0=axs[0,imgcols[0]].get_position().x0; xL1=axs[0,imgcols[N-1]].get_position().x1
    xR0=axs[0,imgcols[N]].get_position().x0; xR1=axs[0,imgcols[2*N-1]].get_position().x1
    yline=max(axs[0,c].get_position().y1 for c in imgcols)+0.03
    for x0,x1,txt in [(xL0,xL1,"lower IC"),(xR0,xR1,"higher IC")]:
        fig.text((x0+x1)/2,yline+0.015,txt,ha="center",va="bottom",fontsize=12,fontweight="bold")
        fig.add_artist(Line2D([x0,x1],[yline,yline],color="k",lw=1.3,transform=fig.transFigure))
    fig.suptitle(title,fontsize=13,fontweight="bold",y=0.99)
    fig.savefig(outpath,bbox_inches="tight",dpi=300); plt.close(fig)
    print("wrote",outpath)

montage(ic,[("IC60","IC60\nleaf rotation\nR²=0.40\nH²=0.04"),
            ("IC37","IC37\nleaf rotation\nR²=0.39\nH²=0.08")],
        "ICs associated with leaf rotation",
        "figures/ica_characterization/montage_geometry.png",N=6)

ne=ic[ic.environment=="Nebraska2025"].copy()
montage(ne,[("IC1","IC1\ndisease\nr=-0.33\nH²=0.55"),
            ("IC57","IC57\ndisease\nr=-0.30\nH²=0.41"),
            ("IC20","IC20\ndisease\nr=-0.23\nH²=0.35"),
            ("IC76","IC76\ndisease\nr=+0.21\nH²=0.33")],
        "ICs most associated with human disease scores",
        "figures/ica_characterization/montage_disease.png",N=5)

# position figure: IC29 boxplot by crop segment; IC38 plain scatter (no color coding)
fig,ax=plt.subplots(1,2,figsize=(6.5, 2.7))
g0=ic.loc[ic.crop_index==0,"IC29"].dropna(); g1=ic.loc[ic.crop_index==1,"IC29"].dropna()
ax[0].boxplot([g0,g1],labels=["crop 1\n(proximal)","crop 2\n(distal)"],showfliers=False,
              widths=.55,patch_artist=True,boxprops=dict(facecolor="#9ecae1",edgecolor="k"),
              medianprops=dict(color="k",linewidth=1.5))
ax[0].set_ylabel("IC29 score")
ax[0].set_title("IC29 by crop segment\nSpearman r(crop_index)=+0.81  (H²≈0, disease |r|≈0)",
                fontsize=10,fontweight="bold")
ax[0].grid(alpha=.25,axis="y")
ax[1].scatter(ic["crop_center_y"],ic["IC38"],s=7,alpha=.3,color="#3a7ca5",edgecolor="none")
ax[1].set_xlabel("crop_center_y"); ax[1].set_ylabel("IC38 score")
ax[1].set_title("IC38 vs vertical crop position\nSpearman r=+0.43  (H²≈0, disease |r|≈0)",
                fontsize=10,fontweight="bold")
ax[1].grid(alpha=.25)
fig.suptitle("Position-artifact ICs track crop placement, not leaf condition",
             fontsize=11,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig("figures/ica_characterization/geometry_position.png",bbox_inches="tight",dpi=300); plt.close(fig)
print("wrote geometry_position.png")
