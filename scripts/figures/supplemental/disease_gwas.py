#!/usr/bin/env python3
"""Three Manhattan/QQ pairs from the full PANICLE marker-level results.

Style sources (local manuscripts, not a new palette):
- This paper: FigS13 lysm_hotspot.R (9 pt sans, italic p, dashed cutoff,
  left/bottom axes, no grid, small rasterized points).
- sorghum-maize-mycobiome/src/figures/figure7/
  ColletGWASCombinedManualRelAbundanceManhattan.R (alternating black/darkgrey,
  chromosome centers and the label 'Chromosome').
- This repo's run_gwas_panicle.py (QQ plotting positions (rank - .5)/N).
All valid markers are plotted; no thinning or significance prefiltering.
"""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[3]
TRAITS = [('human_score', 'Human disease scores'), ('exg_raw', 'Raw ExG'), ('exg_logit', 'Logit ExG')]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gwas-dir', type=Path, default=ROOT/'data/generatable/disease_gwas/gwas')
    ap.add_argument('--out-dir', type=Path, default=ROOT/'figures/supplemental/FigS19_disease_gwas')
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    meta = json.loads((args.gwas_dir/'effective_tests.json').read_text())
    effective = meta['effective_tests']['Me']
    assert effective == 4446367, effective
    cutoff = .05 / effective
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                         'font.size': 9, 'axes.titlesize': 9, 'axes.labelsize': 9,
                         'xtick.labelsize': 9, 'ytick.labelsize': 9,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .5, 'xtick.major.width': .5, 'ytick.major.width': .5,
                         'mathtext.default': 'regular', 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(6.5, 5.7))
    grid = fig.add_gridspec(3, 2, width_ratios=[3.15, 1], left=.085, right=.985,
                           bottom=.09, top=.935, hspace=.8, wspace=.4)
    stats = {'effective_tests': effective, 'p_threshold': cutoff,
             'minus_log10_threshold': -math.log10(cutoff), 'traits': {}}
    for row, (trait, title) in enumerate(TRAITS):
        df = pd.read_csv(args.gwas_dir/'traits'/f'{trait}_marker_pvalues.csv',
                         usecols=['CHROM', 'POS', 'p_value'], dtype={'CHROM': np.int16, 'POS': np.int64})
        assert len(df) == 6422975
        valid = np.isfinite(df.p_value) & df.p_value.gt(0) & df.p_value.le(1)
        if not valid.all():
            raise ValueError(f'{trait}: invalid marker p-values')
        p = df.p_value.to_numpy()
        observed = -np.log10(p)
        ax = fig.add_subplot(grid[row, 0])
        ticks, labels, offset = [], [], 0
        for i, chrom in enumerate(sorted(df.CHROM.unique())):
            mask = df.CHROM.to_numpy() == chrom
            pos = df.POS.to_numpy()[mask]
            ax.scatter(pos + offset, observed[mask], s=.8, c=['black', 'darkgrey'][i % 2],
                       alpha=.8, linewidths=0, rasterized=True)
            ticks.append(offset + pos.max()/2)
            labels.append(str(chrom))
            offset += int(pos.max())
        ax.axhline(-math.log10(cutoff), color='black', ls='--', lw=.4)
        ax.set(xlim=(0, offset), ylim=(0, math.ceil(max(observed.max(), -math.log10(cutoff)) * 1.04)),
               xlabel='Chromosome', ylabel=r'$-\log_{10}(\mathit{p})$')
        ax.set_xticks(ticks, labels)
        ax.set_title(title, loc='left', pad=9)
        ax.text(-.12, 1.085, chr(97+row*2), transform=ax.transAxes, fontsize=12, weight='bold')
        qq = fig.add_subplot(grid[row, 1])
        ordered = np.sort(p)
        expected = -np.log10((np.arange(1, len(p)+1)-.5)/len(p))
        observed_qq = -np.log10(ordered)
        limit = math.ceil(max(expected.max(), observed_qq.max()))
        qq.scatter(expected, observed_qq, s=.8, color='black', alpha=.8, linewidths=0, rasterized=True)
        qq.plot([0, limit], [0, limit], color='black', lw=.5)
        qq.set(xlim=(0, limit), ylim=(0, limit),
               xlabel=r'Expected $-\log_{10}(\mathit{p})$', ylabel=r'Observed $-\log_{10}(\mathit{p})$')
        qq.set_box_aspect(1)
        qq.locator_params(axis='both', nbins=3)
        inflation = float(chi2.isf(np.median(p), 1)/chi2.ppf(.5, 1))
        qq.text(.05, .95, rf'$\lambda={inflation:.2f}$', transform=qq.transAxes, va='top', fontsize=8)
        qq.text(-.42, 1.085, chr(98+row*2), transform=qq.transAxes, fontsize=12, weight='bold')
        stats['traits'][trait] = {'markers': len(p), 'min_p': float(p.min()),
            'significant_markers': int((p <= cutoff).sum()), 'lambda_gc': inflation}
        print(trait, stats['traits'][trait], flush=True)
    fig.savefig(args.out_dir/'disease_gwas.png', dpi=300, facecolor='white')
    fig.savefig(args.out_dir/'disease_gwas.pdf', dpi=300, facecolor='white')
    (args.out_dir/'summary.json').write_text(json.dumps(stats, indent=2)+'\n')


if __name__ == '__main__':
    main()
