#!/usr/bin/env python3
"""Rebuild chr2 leaf-water phenotypes and verify with the current single-marker CLI.

Uses untransformed water fraction and weights; retains the historical pooling
rule (equal-weight mean of each genotype's environment-level phenotypes).
Outputs are kept separately from the historical figure inputs.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import zipfile
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sensitivity-only", action="store_true", help="Run only the homozygote check after a completed primary run")
args = parser.parse_args()
ARCHIVE = ROOT / "data/externalsourcerequired/sorghum_trait_data_v2.2.zip"
OUT = ROOT / "data/generatable/chr2_leaf_water_verification"
OUT.mkdir(parents=True, exist_ok=True)
FRESH = "single_plant_leaf_fresh_weight_g"
DRY = "single_plant_leaf_dry_weight_g"
with zipfile.ZipFile(ARCHIVE) as archive:
    name = next(n for n in archive.namelist() if n.endswith("/observations.tsv"))
    raw = archive.read(name)
    with archive.open(name) as stream:
        obs = pd.read_csv(stream, sep="\t", low_memory=False)
obs.genotype = obs.genotype.astype(str).str.replace(" ", "", regex=False)
obs.value = pd.to_numeric(obs.value, errors="coerce")
selected = obs[obs.canonical_name.isin([FRESH, DRY]) & obs.env_id.isin(["MI2020", "MI2021"])].copy()
wide = selected.groupby(["env_id", "genotype", "canonical_name"]).value.mean().unstack("canonical_name")
missing_pairs = int(wide[[FRESH, DRY]].isna().any(axis=1).sum())
wide = wide.dropna(subset=[FRESH, DRY]).rename(columns={FRESH: "fresh_g", DRY: "dry_g"})
invalid = (wide.fresh_g <= 0) | (wide.dry_g < 0) | (wide.dry_g > wide.fresh_g)
if invalid.any():
    raise ValueError(f"Invalid fresh/dry pairs require review: {wide[invalid]}")
wide["water_fraction"] = (wide.fresh_g - wide.dry_g) / wide.fresh_g
pooled = wide.groupby("genotype").mean().reset_index().assign(env_id="MI2020+MI2021")
phenotypes = pd.concat([wide.reset_index(), pooled], ignore_index=True)
phenotypes.to_csv(OUT / "phenotypes.csv", index=False)
old = pd.read_csv(ROOT / "figures/chr2_gloss_peak/story_biomass_data.csv")
heterozygotes = set(old.loc[old.peak_dose == 1, "genotype"])
phenotypes[~phenotypes.genotype.isin(heterozygotes)].to_csv(OUT / "phenotypes_no_heterozygotes.csv", index=False)
comparison = pooled.merge(old, on="genotype")
qc = {
    "archive": str(ARCHIVE), "archive_sha256": hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
    "observations_sha256": hashlib.sha256(raw).hexdigest(),
    "missing_fresh_dry_pairs": missing_pairs, "invalid_pairs": int(invalid.sum()),
    "paired_genotypes_by_environment": wide.groupby(level=0).size().to_dict(),
    "historical_pooled_overlap": len(comparison),
    "historical_pooled_max_abs_difference": {
        new: float((comparison[new] - comparison[oldcol]).abs().max())
        for new, oldcol in [("fresh_g", "fresh"), ("dry_g", "dry"), ("water_fraction", "water_frac")]
    },
    "pooling": "Mean of genotype-level environment means; water fraction calculated within environment before pooling",
}
(OUT / "input_audit.json").write_text(json.dumps(qc, indent=2) + "\n")
print(json.dumps(qc, indent=2), flush=True)
for trait, suffix, extra in [
    ("water_fraction", "current", []),
    ("fresh_g", "current", []),
    ("dry_g", "current", []),
    ("water_fraction", "pcs_only", ["--no-covariates"]),
    ("water_fraction", "homozygotes", []),
]:
    if args.sensitivity_only and suffix != "homozygotes":
        continue
    phenotype_file = "phenotypes_no_heterozygotes.csv" if suffix == "homozygotes" else "phenotypes.csv"
    subprocess.run([
        sys.executable, "scripts/run_single_marker_test.py", str(OUT / phenotype_file),
        trait, "2:52490664:GGAGT:G", "--group-column", "env_id", "--cpu", "4",
        "--out-file", str(OUT / f"{trait}_{suffix}.csv"), *extra,
    ], cwd=ROOT, check=True)
rows = []
for filename in ["water_fraction_current", "fresh_g_current", "dry_g_current", "water_fraction_pcs_only", "water_fraction_homozygotes"]:
    frame = pd.read_csv(OUT / f"{filename}.csv")
    frame["analysis"] = filename
    frame["alt_homozygote_minus_ref"] = 2 * frame.effect_alt_allele
    frame["homozygote_contrast_se"] = 2 * frame.se
    rows.append(frame)
summary = pd.concat(rows, ignore_index=True)
summary.to_csv(OUT / "summary.csv", index=False)
print(summary[["analysis", "group", "n_observations", "n_ref_homozygote", "n_heterozygote", "n_alt_homozygote", "alt_homozygote_minus_ref", "p_value", "status"]].to_string(index=False))
