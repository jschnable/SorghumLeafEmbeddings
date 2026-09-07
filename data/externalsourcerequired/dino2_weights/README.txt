Place DINOv2 weights here after downloading them yourself.

Default expected use:
  python scripts/extract_embeddings.py ... --backend dino2 --dino2-weights data/externalsourcerequired/dino2_weights

If this directory contains a .pth checkpoint, the script loads
dinov2_vitl14_reg with pretrained=False and applies that checkpoint
with strict=True. If no .pth file is present, torch.hub downloads/caches the
official dinov2_vitl14_reg weights.

Current cleaned-pipeline model:
  dinov2_vitl14_reg

This is the large ViT-L/14 model with registers.

Download/source:
  Official repository: https://github.com/facebookresearch/dinov2
  Use torch.hub model name dinov2_vitl14_reg for current analyses.
  The old tiny-model dinov2_vits14_reg path is not part of this pipeline.
