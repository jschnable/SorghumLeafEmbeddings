Place optional local DINOv2 weights here after downloading them yourself.
The paper companion repository does not redistribute DINOv2 weights.

Default expected use:
  python scripts/extract_embeddings.py ... --backend dino2 --dino2-weights placeholders/dino2_weights

If this directory contains a .pth checkpoint, the script loads
dinov2_vitl14_reg with pretrained=False and applies that checkpoint
with strict=True. Mismatched checkpoints fail instead of loading partially. If
no .pth file is present, torch.hub downloads/caches the
official dinov2_vitl14_reg weights.

Local .pth files must be backbone-only state dicts for dinov2_vitl14_reg. Full
training checkpoints with classifier heads, wrapper modules, or a different
DINOv2 architecture are expected to fail strict loading.

Current cleaned-pipeline model:
  dinov2_vitl14_reg

This is the large ViT-L/14 model with registers. It replaced earlier legacy
exploration with the much smaller dinov2_vits14_reg model.

Download/source:
  Official repository: https://github.com/facebookresearch/dinov2
  Use torch.hub model name dinov2_vitl14_reg for current analyses.
  The old tiny-model dinov2_vits14_reg path is not part of this pipeline.
