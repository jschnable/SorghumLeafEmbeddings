# Six-environment, homozygote-only LysM grain/panicle mass follow-up.
.script_file <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
.repo_root <- dirname(normalizePath(.script_file))
while (!file.exists(file.path(.repo_root, 'scripts', 'extract_embeddings.py'))) {
  .parent <- dirname(.repo_root)
  if (.parent == .repo_root) stop('Cannot locate repository root')
  .repo_root <- .parent
}
library(tidyverse)
library(paletteer)
library(cowplot)
# Refit the six homozygote-only tests while building the figure. Results remain
# in memory; no tests.csv or audit/metadata export is part of the deliverable.
inputs <- file.path(.repo_root, 'figures/supplemental/FigS14_lysm_yield')
df <- read_csv(file.path(inputs, 'phenotypes.csv'), show_col_types=FALSE)
library(reticulate)
py_run_string(r"(
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
)")
vcf <- Sys.getenv('LEAF_GENOTYPE_VCF', unset=file.path(.repo_root, 'data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz'))
tests <- py$fit_lysm_tests(.repo_root, vcf, as.integer(Sys.getenv('LEAF_CPU', '1')))
environments <- c('MI2020','MI2021','NE2020','NE2021','NE2023','NE2025')
titles <- c('Michigan 2020','Michigan 2021','Nebraska 2020','Nebraska 2021','Nebraska 2023','Nebraska 2025')
traits <- c('Panicle dry mass (g)','Panicle dry mass (g)','Seed mass (g)', 'Seed mass (g)','Single-panicle mass (g)','Seed mass (g)')
plots <- lapply(seq_along(environments), function(i) {
  dat <- filter(df, env_id == environments[i]) %>% mutate(allele=factor(allele, levels=c('GG','TT')))
  test <- filter(tests, group == environments[i])
  stopifnot(nrow(dat) == test$n_observations, test$n_heterozygote == 0)
  p_label <- if (test$p_value < .001) sprintf('p = %.2e', test$p_value) else sprintf('p = %.4f', test$p_value)
  counts <- table(dat$allele)
  ggplot(dat, aes(allele, mass_g, fill=allele)) +
    geom_boxplot(width=.5, outlier.size=.6, linewidth=.4) +
    scale_fill_manual(values=as.character(paletteer_d('colorBlindness::Brown2Blue10Steps')[c(3,2)])) +
    scale_x_discrete(labels=c(sprintf('GG\nn = %d', counts[['GG']]), sprintf('TT\nn = %d', counts[['TT']]))) +
    scale_y_continuous(limits=c(0,NA), expand=expansion(mult=c(0,.08))) +
    labs(title=titles[i], subtitle=p_label, x=NULL, y=traits[i]) +
    theme_classic(base_size=9) +
    theme(legend.position='none', plot.title=element_text(size=10,hjust=.5),
          plot.subtitle=element_text(size=9,hjust=.5), axis.text=element_text(color='black'))
})
figure <- plot_grid(plotlist=plots, ncol=3, labels=letters[1:6], label_size=13)
ggsave(file.path(.repo_root, 'figures/supplemental/FigS14_lysm_yield/chr9_1_panicle_wt.png'), figure,
       width=6.5, height=4.6, dpi=300, bg='white')
