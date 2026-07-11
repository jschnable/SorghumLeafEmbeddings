import sys, warnings, time; warnings.filterwarnings("ignore"); sys.path.insert(0,"scripts")
import numpy as np
from embedding_io import read_embedding_table, split_feature_metadata_columns
from sklearn.preprocessing import StandardScaler
from sklearn.utils.extmath import randomized_svd
t0=time.time()

df = read_embedding_table("data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz")
feat,_ = split_feature_metadata_columns(df)
X = df[feat].to_numpy(np.float64)
N, p = X.shape
Xs = StandardScaler().fit_transform(X)            # z-score each feature
# independence units (crops are NOT independent): unique source images
img = df["source_image_path"] if "source_image_path" in df else df["image_id"]
n_img = img.nunique(); n_geno = df["genotype"].nunique()
print(f"N(crops)={N}  p(features)={p}  n_images={n_img}  n_genotypes={n_geno}")

# full covariance eigenvalues (sorted desc)
C = np.cov(Xs, rowvar=False)
lam = np.sort(np.linalg.eigvalsh(C))[::-1]
lam = np.clip(lam, 1e-12, None)
print(f"[eig] top5={np.round(lam[:5],2)}  ... tail5={np.round(lam[-5:],4)}  (computed {time.time()-t0:.0f}s)")

# ---------- 1) SCREE / ELBOW (max distance to chord) ----------
def elbow(y):
    n=len(y); x=np.arange(n)
    x0,y0,x1,y1=x[0],y[0],x[-1],y[-1]
    d=np.abs((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)/np.hypot(y1-y0,x1-x0)
    return int(np.argmax(d))+1
k_elbow_full = elbow(lam)
k_elbow_200  = elbow(lam[:200])
k_kaiser     = int((lam>1).sum())
print(f"\n[1] SCREE/ELBOW: full-curve elbow k={k_elbow_full} | elbow(top200)={k_elbow_200} | Kaiser(eig>1)={k_kaiser}")

# ---------- 2) HORN PARALLEL ANALYSIS ----------
K=300; NP=15
rng=np.random.RandomState(0)
null=np.zeros((NP,K))
for b in range(NP):
    Xp=Xs.copy()
    for j in range(p): Xp[:,j]=Xp[rng.permutation(N),j]
    _,s,_=randomized_svd(Xp-Xp.mean(0),n_components=K,random_state=b)
    null[b]=(s**2)/(N-1)
null_p95=np.percentile(null,95,axis=0)
obs=lam[:K]
above=obs>null_p95
k_horn=int(np.argmax(~above)) if (~above).any() else K   # first crossing
print(f"[2] HORN parallel analysis: k={k_horn}  (obs eig falls below 95th-pct null at component {k_horn+1})")

# ---------- 3) MDL / BIC (Wax-Kailath, GIFT-style) ----------
def order_select(lam, N, p):
    aic=np.full(p,np.inf); mdl=np.full(p,np.inf)
    for k in range(p):
        tail=lam[k:]; m=len(tail)
        if m<1: break
        a=tail.mean(); g=np.exp(np.log(tail).mean())
        llf=N*m*np.log(a/g)            # >=0
        nfree=k*(2*p-k)
        aic[k]=2*llf+2*nfree
        mdl[k]=llf+0.5*nfree*np.log(N)
    return int(np.argmin(aic)), int(np.argmin(mdl))
k_aic, k_mdl = order_select(lam, N, p)
# effective-N versions (crops are correlated -> GIFT corrects N down)
k_aic_img, k_mdl_img = order_select(lam, n_img, p)
k_aic_g,   k_mdl_g   = order_select(lam, n_geno, p)
print(f"[3] MDL/AIC (Wax-Kailath):")
print(f"      N=crops({N}):    MDL k={k_mdl}   AIC k={k_aic}")
print(f"      N=images({n_img}): MDL k={k_mdl_img}   AIC k={k_aic_img}")
print(f"      N=genos({n_geno}):  MDL k={k_mdl_g}   AIC k={k_aic_g}")
print(f"\nelapsed {time.time()-t0:.0f}s")
