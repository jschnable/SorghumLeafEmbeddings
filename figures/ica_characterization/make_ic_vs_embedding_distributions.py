"""Per-feature heritability and disease-correlation distributions: 2048 raw
embedding dims vs 90 ICs (Nebraska). See ICA_methods_and_results.md section 7."""
import sys, warnings; warnings.filterwarnings("ignore"); sys.path.insert(0, "scripts")
import pandas as pd, numpy as np, re, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from embedding_io import read_embedding_table, image_key

Hic = pd.read_csv("data/generatable/ic_blues_horn90/heritability_Nebraska2025.csv")
Hemb = pd.read_csv("data/generatable/blues/nebraska_sam3_embeddings_2016crop/heritability_Nebraska2025.csv")
h2_ic, h2_emb = Hic.broad_sense_h2.dropna().to_numpy(), Hemb.broad_sense_h2.dropna().to_numpy()

HS = pd.read_csv("data/provided/human_disease_scores.csv"); HS["ik"] = HS.image_id.map(image_key)
HS = HS.dropna(subset=["human_score"]).groupby("ik", as_index=False).agg(
    human_score=("human_score", "mean"), environment=("environment", "first"))
HS = HS[HS.environment == "Nebraska2025"]
def absr(path, cols):
    t = read_embedding_table(path); key = None
    for c in ("source_image_path", "image_path"):
        if c in t.columns:
            k = t[c].map(image_key); key = k if key is None else key.fillna(k)
    img = t.assign(ik=key).groupby("ik", as_index=False)[cols].mean().merge(HS, on="ik")
    return img[cols].corrwith(img["human_score"], method="spearman").abs().dropna().to_numpy()
ic_p = "data/generatable/dimreduction_horn90/ic_scores.csv"
iccols = [c for c in read_embedding_table(ic_p).columns if re.fullmatch(r"IC\d+", c)]
emb_p = "data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz"
embcols = [c for c in read_embedding_table(emb_p).columns
           if c.startswith("embedding_mean_") or c.startswith("embedding_std_")]
r_ic, r_emb = absr(ic_p, iccols), absr(emb_p, embcols)

plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold"})
fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.0))
def panel(a, d_ic, d_emb, ylab, title):
    parts = a.violinplot([d_emb, d_ic], showmedians=True, widths=.85)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(["#b0772d", "#2e6e8e"][i]); pc.set_alpha(.6); pc.set_edgecolor("k")
    parts["cmedians"].set_color("k")
    a.set_xticks([1, 2]); a.set_xticklabels(["2048\nembeddings", "90\nICs"])
    a.set_ylabel(ylab); a.set_title(title); a.grid(alpha=.25, axis="y")
    for i, d in enumerate([d_emb, d_ic], 1): a.text(i, np.median(d), f"  med {np.median(d):.2f}", va="center", fontsize=9)
panel(ax[0], h2_ic, h2_emb, "broad-sense heritability (H²)", "Heritability: embeddings vs ICs")
panel(ax[1], r_ic, r_emb, "|Spearman r| with human disease", "Disease correlation: embeddings vs ICs")
fig.suptitle("Per-feature distributions: 2048 raw embedding dims vs 90 ICs (Nebraska)", fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("figures/ica_characterization/ic_vs_embedding_distributions.png", bbox_inches="tight", dpi=300)
print("wrote ic_vs_embedding_distributions.png")
