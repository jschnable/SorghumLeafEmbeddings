import sys,warnings; warnings.filterwarnings("ignore"); sys.path.insert(0,"scripts")
import pandas as pd, numpy as np, re, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from embedding_io import read_embedding_table, image_key
from sklearn.linear_model import Ridge, RidgeCV

HS=pd.read_csv("data/provided/human_disease_scores.csv"); HS["ik"]=HS.image_id.map(image_key)
HS=HS.dropna(subset=["human_score"]).groupby("ik",as_index=False).agg(
     human_score=("human_score","mean"),environment=("environment","first"))
HS=HS[HS.environment=="Nebraska2025"]
ic=read_embedding_table("data/generatable/dimreduction_horn90/ic_scores.csv")
iccols=sorted([c for c in ic.columns if re.fullmatch(r"IC\d+",c)],key=lambda c:int(c[2:]))
ic["ik"]=ic.source_image_path.map(image_key)
img=ic.groupby("ik",as_index=False).agg({**{c:"mean" for c in iccols},"genotype":"first"})
m=img.merge(HS,on="ik",how="inner"); y=m.human_score.to_numpy(); g=m.genotype.to_numpy()
ranked=m[iccols].corrwith(m.human_score,method="spearman").abs().sort_values(ascending=False).index.tolist()
Xic=m[ranked].to_numpy(float)
emb=read_embedding_table("data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz")
embc=[c for c in emb.columns if c.startswith("embedding_mean_") or c.startswith("embedding_std_")]
emb["ik"]=emb.source_image_path.map(image_key)
Xemb=emb.groupby("ik",as_index=False)[embc].mean().set_index("ik").loc[m.ik].to_numpy(float)

A=np.logspace(-2,4,13); NREP=8; NF=5; ks=np.arange(1,91); genos=pd.unique(g)
# alpha per k (chosen once on full data), and for embeddings
alpha_k=[RidgeCV(alphas=A).fit(Xic[:,:k],y).alpha_ for k in ks]
alpha_e=RidgeCV(alphas=A).fit(Xemb,y).alpha_
def oof(X,a,folds):
    o=np.full(len(y),np.nan)
    for f in range(NF):
        te=folds==f; tr=~te
        if te.sum()<2: continue
        o[te]=Ridge(alpha=a).fit(X[tr],y[tr]).predict(X[te])
    mk=~np.isnan(o); return 1-((y[mk]-o[mk])**2).sum()/((y[mk]-y[mk].mean())**2).sum()
curve=np.zeros((NREP,len(ks))); embR=np.zeros(NREP)
for rep in range(NREP):
    rng=np.random.RandomState(rep); fo=dict(zip(genos,rng.randint(0,NF,len(genos))))
    folds=np.array([fo[x] for x in g])
    for j,k in enumerate(ks): curve[rep,j]=oof(Xic[:,:k],alpha_k[j],folds)
    embR[rep]=oof(Xemb,alpha_e,folds)
cm=curve.mean(0); cs=curve.std(0); em=embR.mean(); es=embR.std()
np.savez(sys.argv[1] if len(sys.argv)>1 else "/tmp/curve.npz",ks=ks,cm=cm,cs=cs,em=em,es=es)
for k in [1,8,20,40,90]: print(f"top {k:2d} ICs: R2={cm[k-1]:.3f}±{cs[k-1]:.3f}")
print(f"2048 embeddings: R2={em:.3f}±{es:.3f}")
fig,ax=plt.subplots(figsize=(6.5, 4.0))
ax.fill_between(ks,cm-cs,cm+cs,color="#2e6e8e",alpha=.2)
ax.plot(ks,cm,color="#2e6e8e",lw=2.3,label="top-k ICs (ranked by |disease r|)")
ax.axhspan(em-es,em+es,color="#b0772d",alpha=.18)
ax.axhline(em,color="#b0772d",lw=1.8,ls="--",label=f"2048 embeddings ({em:.2f} ± {es:.2f})")
ax.axvline(8,color="grey",ls=":",lw=1.2); ax.text(9,0.07,"|r|>0.15\ncutoff (8 ICs)",fontsize=8,color="grey")
ax.set_xlabel("number of ICs included (descending |disease correlation|)")
ax.set_ylabel("out-of-fold R²  (genotype-grouped 5-fold CV)")
ax.set_title("Disease prediction vs number of ICs (Nebraska human scores)",fontweight="bold")
ax.set_xlim(1,90); ax.set_ylim(0,max(em+es,(cm+cs).max())+0.05); ax.grid(alpha=.25); ax.legend(loc="lower right")
fig.tight_layout(); fig.savefig("figures/ica_characterization/disease_prediction_vs_n_ics.png",bbox_inches="tight",dpi=300)
print("wrote figures/ica_characterization/disease_prediction_vs_n_ics.png")
