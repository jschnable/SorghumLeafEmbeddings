# `data/externalsourcerequired/`

Large assets that are **not committed to GitHub** but are required to regenerate
the files in `data/generatable/`. Obtain each from its external source and place
it in the subfolder named below; the scripts default to these locations.

- **Raw images** — the full per-environment leaf image set. Distributed
  separately (too large for GitHub). Used by `scripts/extract_embeddings.py` to
  produce embeddings.
- `vcf/` — the sorghum marker **VCF or PLINK** files for GWAS. See
  `vcf/README.txt`. Used by `scripts/run_gwas_panicle.py`.
- `sam3_weights/` — the Hugging Face **SAM3** model files
  (`facebook/sam3`). See `sam3_weights/README.txt`.
- `dino2_weights/` — optional local **DINOv2** (`dinov2_vitl14_reg`) `.pth`
  checkpoint; if absent, `torch.hub` downloads the official weights. See
  `dino2_weights/README.txt`.

None of these should be committed to this repository.
