import sys, re, warnings; warnings.filterwarnings("ignore"); sys.path.insert(0,"scripts")
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from embedding_io import image_key

T=pd.read_csv("data/generatable/dimreduction_horn90/ic_consolidated_nebraska.csv",index_col=0)
iccols=sorted([c for c in T.index],key=lambda c:int(c[2:]))

# --- per-IC |Spearman r| with log1p(ExG), Nebraska, image-level (matches human line) ---
ic=pd.read_csv("data/generatable/dimreduction_horn90/ic_scores.csv")
ne=ic[ic.environment=="Nebraska2025"].copy()
key=None
for col in ("source_image_path","image_path"):
    if col in ne: 
        k=ne[col].map(image_key); key=k if key is None else key.fillna(k)
ne["image_key"]=key
icc=[c for c in ne.columns if re.fullmatch(r"IC\d+",c)]
img=ne.groupby("image_key").agg({**{c:"mean" for c in icc},"disease_exg":"mean"})
exg_absr=img[icc].corrwith(img["disease_exg"],method="spearman").abs()
T["exg_absr"]=exg_absr
print("ExG |r| (Nebraska, image-level): max=%.3f  #ICs p-ish strong(>0.15)=%d"%(
      exg_absr.max(),(exg_absr>0.15).sum()))
print("human |r| max=%.3f"%T.disease_absr.max())

h2=T["h2"].to_numpy(); dis=T["disease_absr"].to_numpy(); exg=T["exg_absr"].to_numpy()
geom=T["geom_R2"].to_numpy()

plt.rcParams.update({"font.size":11,"axes.titlesize":12,"axes.titleweight":"bold",
                     "figure.dpi":120,"savefig.dpi":300})
fig=plt.figure(figsize=(6.5, 5.0))
gs=fig.add_gridspec(2,2,hspace=.28,wspace=.26)
axAB=fig.add_subplot(gs[0,0]); axC=fig.add_subplot(gs[0,1]); axD=fig.add_subplot(gs[1,:])

# ---- A+B merged: dual y-axis ----
thr=np.linspace(0.75,0,300)
l1,=axAB.plot(thr,[(h2>=t).sum() for t in thr],color="#1b6ca8",lw=2.6,label="H² (heritability)")
axAB.set_xlim(0.75,0); axAB.set_xlabel("threshold  (H²  or  |Spearman r|)")
axAB.set_ylabel("# ICs ≥ H² threshold",color="#1b6ca8")
axAB.tick_params(axis="y",labelcolor="#1b6ca8"); axAB.grid(alpha=.22)
ax2=axAB.twinx()
l2,=ax2.plot(thr,[(dis>=t).sum() for t in thr],color="#b8480f",lw=2.4,label="|r| human disease score")
l3,=ax2.plot(thr,[(exg>=t).sum() for t in thr],color="#2e7d32",lw=2.4,ls="--",
             label="|r| log1p(ExG) disease")
ax2.set_ylabel("# ICs ≥ disease |r| threshold",color="#b8480f")
ax2.tick_params(axis="y",labelcolor="#b8480f")
axAB.set_title("A  ICs passing heritability / disease thresholds")
axAB.legend(handles=[l1,l2,l3],frameon=False,fontsize=9,loc="upper left")

# ---- C scatter ----
sc=axC.scatter(h2,dis,c=geom,cmap="viridis",s=46,edgecolor="k",linewidth=.3,vmin=0,vmax=.75)
cb=fig.colorbar(sc,ax=axC); cb.set_label("variance from\nextraction geometry (R²)")
for lab in ["IC1","IC19","IC29","IC57","IC26"]:
    axC.annotate(lab,(T.loc[lab,"h2"],T.loc[lab,"disease_absr"]),
                 xytext=(4,4),textcoords="offset points",fontsize=8)
axC.set_xlabel("broad-sense heritability (H²)"); axC.set_ylabel("|disease correlation| (human)")
axC.set_title("B  Heritability vs disease signal (per IC)"); axC.grid(alpha=.22)

# ---- D violins ----
data=[T["geno_pct"],T["spatial_pct"],T["device_pct"],T["resid_pct"]]
labels=["Genotype","Spatial\n(row+col+block)","Device","Residual"]
parts=axD.violinplot([d.to_numpy() for d in data],showmedians=True,widths=.85)
for pc in parts["bodies"]: pc.set_facecolor("#5aa469"); pc.set_alpha(.65); pc.set_edgecolor("k")
parts["cmedians"].set_color("k")
axD.set_xticks(range(1,5)); axD.set_xticklabels(labels)
axD.set_ylabel("proportion of variance explained (%)")
axD.set_title("C  Variance partition across all 90 ICs"); axD.grid(alpha=.22,axis="y")
for i,d in enumerate(data,1): axD.text(i,np.median(d)+2,f"med {np.median(d):.0f}%",ha="center",fontsize=8)

fig.suptitle("Characterization of 90 Horn-selected independent components — Nebraska 2025",
             fontsize=13,fontweight="bold",y=.995)
fig.savefig("figures/ica_characterization/ica_characterization.png",bbox_inches="tight")
T.round(4).to_csv("data/generatable/dimreduction_horn90/ic_consolidated_nebraska.csv")
print("wrote figure + updated consolidated table (added exg_absr)")
