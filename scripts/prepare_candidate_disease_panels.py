#!/usr/bin/env python3
"""Prepare current-model disease inputs used by the JA and LysM supplements.

Run from the repository root. Raw expression panels are deliberately unchanged.
Only final plot inputs are exported to data/provided/candidate_disease_panels;
single-marker calculation outputs remain under data/generatable.
"""
from pathlib import Path
import argparse
import subprocess
import sys

import numpy as np
import pandas as pd
from panicle.data.loaders import load_genotype_file
from run_single_marker_test import find_marker_index, marker_frame

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', type=int, default=4)
    parser.add_argument('--vcf', type=Path, default=ROOT / 'data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz')
    args = parser.parse_args()
    output = ROOT / 'data/provided/candidate_disease_panels'
    work = ROOT / 'data/generatable/candidate_disease_panels'
    output.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    base = ROOT / 'figures/supplemental/lysm_hotspot'
    jobs = [
        ('4:4724594:G:C', 'human_score_blue', 'human_score_blue_nebraska.csv'),
        ('9:62301540:T:A', 'human_score_blue', 'human_score_blue_nebraska.csv'),
        ('9:1768703:G:T', 'human_score_blue', 'human_score_blue_nebraska.csv'),
        ('9:1768703:G:T', 'exg_logit_blue', 'exg_logit_blue_nebraska.csv'),
    ]
    tests = []
    for marker, trait, filename in jobs:
        result_path = work / f'{marker}_{trait}.csv'
        subprocess.run([
            sys.executable, str(ROOT / 'scripts/run_single_marker_test.py'),
            str(base / filename), trait, marker, '--genotype', str(args.vcf),
            '--cpu', str(args.cpu), '--out-file', str(result_path),
        ], cwd=ROOT, check=True)
        result = pd.read_csv(result_path)
        assert len(result) == 1 and result.iloc[0]['status'] == 'tested'
        tests.append(result)
    pd.concat(tests, ignore_index=True).assign(group='Nebraska2025').to_csv(output / 'tests.csv', index=False)

    geno, ids, genome_map = load_genotype_file(str(args.vcf), file_format='vcf', precompute_alleles=False)
    calls = pd.DataFrame({'genotype': [str(g).replace(' ', '') for g in ids]})
    cov = pd.read_csv(ROOT / 'data/provided/gwas_covariates_leaf_area_flowering_time.csv')
    cov.genotype = cov.genotype.str.replace(' ', '', regex=False)
    eligible = set(cov.dropna(subset=['mask_pixels_blue', 'days_to_flower_blue']).genotype)
    for marker in dict.fromkeys(j[0] for j in jobs):
        index = find_marker_index(marker_frame(genome_map), marker)
        dosage = geno.subset_markers(np.array([index])).to_numpy()[:, 0]
        calls[marker] = pd.Series(dosage).map({0: '0/0', 2: '1/1'})
    calls = calls[calls.genotype.isin(eligible)]
    # Verify that the plotted homozygote groups reproduce each test's sample counts.
    for (marker, trait, filename), result in zip(jobs, tests):
        phen = pd.read_csv(base / filename)
        phen.genotype = phen.genotype.str.replace(' ', '', regex=False)
        joined = calls.merge(phen.dropna(subset=[trait]), on='genotype', validate='one_to_one')
        for value, column in [('0/0', 'n_ref_homozygote'), ('1/1', 'n_alt_homozygote')]:
            assert int((joined[marker] == value).sum()) == int(result.iloc[0][column])
    calls.to_csv(output / 'genotypes.csv', index=False)


if __name__ == '__main__':
    main()
