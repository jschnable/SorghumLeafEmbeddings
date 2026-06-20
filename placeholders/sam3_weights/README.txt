Place SAM3 model weights/config files here after downloading them yourself.
The paper companion repository does not redistribute SAM3 weights.

Default expected use:
  python scripts/extract_embeddings.py ... --backend sam3 --sam3-weights placeholders/sam3_weights

Current cleaned-pipeline model:
  Hugging Face model id: facebook/sam3
  Loader: transformers.Sam3Model and transformers.Sam3Processor

Download source:
  Official Hugging Face model repository: https://huggingface.co/facebook/sam3

For offline use, download the full model repository into this directory, for example:
  huggingface-cli download facebook/sam3 --local-dir placeholders/sam3_weights

SAM 3.x model code changes over time. Use a transformers version that supports
Sam3Model and Sam3Processor, and keep the model repository files together from
one downloaded snapshot. Do not mix config/tokenizer files and weights from
different SAM3 snapshots.
