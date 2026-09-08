#!/usr/bin/env python3
"""Fit the three Nebraska disease BLUEs with the manuscript's lme4 model.

Read mask areas once per original photograph, not once per embedding crop.
Use the shared BLUE implementation; record model terms and cohort at each step.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import calculate_blues as cb
from embedding_io import image_key

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, default=ROOT/'data/generatable/disease_gwas')
    parser.add_argument('--area-npz', type=Path, default=ROOT/'data/generatable/embeddings/dino2_all3_embeddings_2016crop_float32.npz')
    parser.add_argument('--winsor-order', choices=['plot', 'image'], default='plot')
    a = parser.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    with np.load(a.area_npz, allow_pickle=False) as z:
        meta = json.loads(z['metadata_json'].item())
    area = pd.DataFrame(meta['data'], columns=meta['columns'])
    area['image_key'] = area.source_image_path.map(image_key)
    area = area[['image_key', 'mask_pixels']].drop_duplicates()
    if area.image_key.duplicated().any():
        raise ValueError('Conflicting mask areas among crops from the same image')
    inputs = [ROOT/'data/provided/human_disease_scores.csv', ROOT/'data/provided/exg_ratings.csv']
    audit = {'winsor_order': a.winsor_order, 'winsor_quantiles': [0.01, 0.99],
             'logit_epsilon': 5e-5, 'mask_metadata_sha256': hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest(), 'mask_area_npz': str(a.area_npz.relative_to(ROOT)),
             'input_sha256': {}, 'traits': {}}
    for p in [*inputs, ROOT/'data/provided/field_image_metadata.csv', ROOT/'data/provided/image_ids_exclude.csv', ROOT/'data/provided/gwas_covariates_leaf_area_flowering_time.csv']:
        audit['input_sha256'][str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    all_blues = []
    for name, source, source_col, transform in [
        ('human_score', inputs[0], 'human_score', False),
        ('exg_raw', inputs[1], 'ExG_P20_disease_pct', False),
        ('exg_logit', inputs[1], 'ExG_P20_disease_pct', True),
    ]:
        args = argparse.Namespace(scores=source, metadata=ROOT/'data/provided/field_image_metadata.csv',
            environment='Nebraska2025', image_col='image_id', trait_regex='^'+source_col+'$',
            winsor_strength=0.01 if a.winsor_order == 'plot' else 0,
            spatial_cols='row,column,block', exclude_list=ROOT/'data/provided/image_ids_exclude.csv',
            metadata_optional=False, vc_cpu=1, progress_every=1, no_plot_averaging=False)
        data, _ = cb.load_data(args)
        data = data.drop(columns=['mask_pixels']).merge(area, on='image_key', how='inner', validate='one_to_one')
        data = data.rename(columns={source_col: name})
        data = data.dropna(subset=['mask_pixels', name])
        if transform:
            data[name] = cb.logit_transform(data[name].to_numpy())
        if a.winsor_order == 'image':
            data[name] = cb.winsorize(data[name].to_numpy(), .01)
        plots, fixed = cb.model_plot_means(data, [name], args)
        random = cb.model_random_effects(plots, args, genotype_random=False)
        assert fixed == ['mask_pixels_scaled']
        assert random == ['row', 'column', 'block', 'device']
        blues = cb.run_lme4_blues(plots, [name], random, fixed, args)
        assert not blues[name].isna().any()
        all_blues.append(blues)
        plots.to_csv(a.out_dir/f'{name}_plot_inputs.csv', index=False)
        audit['traits'][name] = {'images': len(data), 'plot_observations': len(plots),
            'blue_genotypes': len(blues), 'fixed_effects': ['genotype', *fixed], 'random_effects': random}
    combined = all_blues[0]
    for frame in all_blues[1:]:
        combined = combined.merge(frame, on='genotype', how='outer', validate='one_to_one')
    combined.to_csv(a.out_dir/'blues_Nebraska2025.csv', index=False)
    (a.out_dir/'blue_metadata.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps(audit['traits'], indent=2))


if __name__ == '__main__':
    main()
