SAM3 checkpoint used for the manuscript embeddings

Default expected use:
  python scripts/extract_embeddings.py ... --backend sam3 --sam3-weights data/externalsourcerequired/sam3_weights

Model:
  Hugging Face model id: facebook/sam3
  Loader: transformers.Sam3Model and transformers.Sam3Processor

Exact revision:
  3c879f39826c281e95690f02c7821c4de09afae7
  https://huggingface.co/facebook/sam3/tree/3c879f39826c281e95690f02c7821c4de09afae7

Checkpoint:
  model.safetensors (3,439,938,512 bytes)
  SHA-256: 6d06f0a5f84e435071fe6603e61d0b4cc7b40e0d39d487cfd4d67d8cc11cc14a

Required companion files: config.json, processor_config.json, tokenizer.json,
special_tokens_map.json, merges.txt and vocab.json. All seven local files match
this Hugging Face revision. checksums.sha256 identifies each file exactly.

Download from the repository root using an authenticated Hugging Face account
with access to facebook/sam3:

  python -c "from huggingface_hub import snapshot_download; snapshot_download(repo_id='facebook/sam3', revision='3c879f39826c281e95690f02c7821c4de09afae7', local_dir='data/externalsourcerequired/sam3_weights', allow_patterns=['config.json', 'model.safetensors', 'processor_config.json', 'tokenizer.json', 'special_tokens_map.json', 'merges.txt', 'vocab.json'])"

Verify from this directory:
  sha256sum -c checksums.sha256

Use the Python dependency versions in requirements.txt with the documented loader.
