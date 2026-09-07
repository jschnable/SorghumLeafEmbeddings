"""Prepare the DINOv2 hotspot inventory and sample counts for retained PheWAS tests.

Uses sliding windows and their GWAS inputs. Sample counts are from the saved
paper PheWAS runs.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/provided/manuscript_tables'
OUT.mkdir(parents=True, exist_ok=True)
windows = pd.read_csv(ROOT / 'figures/main/Fig3_hotspots/dino2_hits_per_100kb.csv.gz')
windows = windows.loc[windows.n_distinct_hits >= 10].sort_values(['CHROM', 'window_start'])
regions = []
for row in windows.itertuples():
    if regions and regions[-1][0] == row.CHROM and row.window_start <= regions[-1][2]:
        regions[-1][2] = max(regions[-1][2], row.window_end)
    else:
        regions.append([row.CHROM, row.window_start, row.window_end])
hits = pd.read_csv(ROOT / 'data/generatable/gwas/embedding_ne_dino2_2016crop_with_cov/significant_markers.csv')
rows = []
for chrom, start, end in regions:
    subset = hits.loc[(hits.CHROM == chrom) & hits.POS.between(start, end)]
    lead = subset.loc[subset.p_value.idxmin()]
    rows.append(dict(chrom=chrom, start_bp=start, end_bp=end,
                     embeddings=subset.trait.nunique(), peak_marker=lead.MARKER,
                     peak_p_value=lead.p_value))
assert len(rows) == 8
pd.DataFrame(rows).to_csv(OUT / 'dino2_hotspots.csv', index=False)

files = sorted((ROOT / 'data/generatable/phwas').glob('*_phwas_results.csv'))
tests = pd.concat([pd.read_csv(path).assign(source=path.name) for path in files])
groups = tests.groupby(['trait', 'location'])
assert groups.n_observations.nunique().eq(1).all(), 'Sample counts differ across markers'
counts = groups.n_observations.first().reset_index()
# The two Michigan leaf-number traits do not pass the retained minimum n=30.
counts = counts.loc[counts.n_observations >= 30].sort_values(['trait', 'location'])
assert len(counts) == 121 and counts.trait.nunique() == 48
counts.to_csv(OUT / 'phewas_sample_counts.csv', index=False)
print(f'Wrote {len(rows)} DINOv2 hotspots and {len(counts)} PheWAS sample counts to {OUT}')
