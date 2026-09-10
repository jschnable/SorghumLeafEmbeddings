#!/usr/bin/env python3
"""Prepare the six homozygote-only association tests plotted in FigS14."""
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
import pandas as pd

def fit_lysm_tests(root, genotype=None, cpu=1):
    root = Path(root)
    sys.path.insert(0, str(root / 'scripts'))
    import run_phwas_panicle as phwas
    from marker_utils import find_marker_index, marker_frame
    from panicle.data.loaders import load_genotype_file
    vcf = Path(genotype) if genotype else root / 'data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz'
    geno, ids, mapping = load_genotype_file(str(vcf), file_format='vcf', precompute_alleles=False)
    ids = [str(g).replace(' ', '') for g in ids]
    id_rows = {g: i for i, g in enumerate(ids)}
    idx = find_marker_index(marker_frame(mapping), '9:1768703:G:T')
    marker = geno.subset_markers(np.array([idx]))
    marker_map = mapping.subset_markers(np.array([idx]))
    dosage = pd.Series(marker.to_numpy()[:, 0].astype(float), index=ids)
    cov_columns = ['mask_pixels_blue', 'days_to_flower_blue']
    cov = phwas.load_covariates(root / 'data/provided/gwas_covariates_leaf_area_flowering_time.csv', ids, cov_columns)
    phenotypes = pd.read_csv(root / 'figures/supplemental/FigS14_lysm_yield/phenotypes.csv')
    if phenotypes.duplicated(['env_id', 'genotype']).any():
        raise ValueError('Figure phenotypes must have one row per genotype and environment')
    args = SimpleNamespace(cache_dir=None, n_pcs=5, max_line=5000, cpu=int(cpu),
                           recompute_cache=False, lrt_solver='GEMMA', lrt_batch_size=2048,
                           screen_threshold=0.0005)  # original single-marker refinement gate
    rows = []
    for environment, frame in phenotypes.groupby('env_id', sort=True):
        frame = frame.set_index('genotype').sort_index()
        if not set(frame.index).issubset(id_rows):
            raise ValueError('Figure genotype missing from marker panel')
        dose = dosage.loc[frame.index]
        if not dose.isin([0, 2]).all() or not np.isfinite(frame.mass_g).all():
            raise ValueError('Expected finite homozygote-only figure phenotypes')
        if not np.isfinite(cov.loc[frame.index, cov_columns]).all().all():
            raise ValueError('Figure population has missing model covariates')
        expected = dose.map({0: 'GG', 2: 'TT'})
        if not expected.equals(frame.allele):
            raise ValueError('Figure allele labels do not match loaded marker calls')
        counts = phwas.marker_counts(dose.to_numpy(), np.ones(len(frame), dtype=bool))
        if len(frame) < 30 or min(counts['n_ref_homozygote'], counts['n_alt_homozygote']) < 3:
            raise ValueError('Figure environment fails original sample-count gates')
        indices = np.array([id_rows[g] for g in frame.index])
        result = phwas.run_group(['mass_g'], frame, geno, marker, marker_map, mapping,
                                 ids, indices, cov, cov_columns, args)['mass_g']
        rows.append(dict(group=environment, n_observations=len(frame), **counts,
                         p_value=phwas.result_value(result.pvalues),
                         effect_alt_allele=phwas.result_value(result.effects),
                         se=phwas.result_value(result.se)))
    return pd.DataFrame(rows)

if __name__ == '__main__':
    import argparse
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genotype', type=Path)
    parser.add_argument('--cpu', type=int, default=1)
    parser.add_argument('--out', type=Path, default=root / 'figures/supplemental/FigS14_lysm_yield/association_tests.csv')
    args = parser.parse_args()
    tests = fit_lysm_tests(root, args.genotype, args.cpu)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tests.to_csv(args.out, index=False)
    print(tests.to_string(index=False))
