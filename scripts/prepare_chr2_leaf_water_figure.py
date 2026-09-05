#!/usr/bin/env python3
"""Export verified raw phenotypes and current-model tests for the GDSL supplement."""
from pathlib import Path
import numpy as np
import pandas as pd
from panicle.data.loaders import load_genotype_file
from run_single_marker_test import marker_frame, find_marker_index

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "data/generatable/chr2_leaf_water_verification"
out = ROOT / "data/provided/chr2_leaf_water_figure"
out.mkdir(parents=True, exist_ok=True)
phenotypes = pd.read_csv(source / "phenotypes.csv")
markers = pd.read_csv(ROOT / "figures/chr2_gloss_peak/story_biomass_data.csv")[["genotype", "peak_dose"]]
covariates = pd.read_csv(ROOT / "data/provided/gwas_covariates_leaf_area_flowering_time.csv")
eligible = covariates.dropna(subset=["mask_pixels_blue", "days_to_flower_blue"]).genotype
data = phenotypes.merge(markers, on="genotype", validate="many_to_one")
data = data[data.genotype.isin(eligible) & data.peak_dose.notna()]
tests = pd.read_csv(source / "summary.csv")
for row in tests.query("analysis == 'water_fraction_current'").itertuples():
    group = data[data.env_id == row.group]
    assert len(group) == row.n_observations, (row.group, len(group), row.n_observations)
    assert int((group.peak_dose == 2).sum()) == row.n_alt_homozygote
data.to_csv(out / "phenotypes.csv", index=False)
tests.to_csv(out / "tests.csv", index=False)
print(data.groupby(["env_id", "peak_dose"]).size())

# Match the expression display to the current VCF and complete-covariate sample.
vcf = ROOT / "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
geno, ids, genome_map = load_genotype_file(vcf, file_format="vcf", precompute_alleles=False)
idx = find_marker_index(marker_frame(genome_map), "4:65447981:G:A")
dose = geno.subset_markers(np.array([idx])).to_numpy()[:, 0]
expr_gt = pd.DataFrame({"genotype": [str(i).replace(" ", "") for i in ids],
                        "4:65447981:G:A": pd.Series(dose).map({0: "0/0", 2: "1/1"})})
expr_gt = expr_gt[expr_gt.genotype.isin(eligible)]
expr_gt = expr_gt[expr_gt["4:65447981:G:A"].isin(["0/0", "1/1"])]
expr_gt.to_csv(out / "chr4_expression_genotypes.csv", index=False)
