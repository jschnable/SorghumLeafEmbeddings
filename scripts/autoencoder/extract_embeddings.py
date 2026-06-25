"""
Extract latent embeddings from a trained autoencoder model.

Usage:
    python src/autoencoder/extract_embeddings.py --model_dir models/model_name/ --output_dir models/model_name/
"""

import os
import sys
import argparse
from pathlib import Path
import pandas as pd
import torch
import yaml
from tqdm import tqdm

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent))

from model import ConvolutionalAutoencoder
from data_utils import create_dataloaders_from_config


def extract_and_save_embeddings(
    model,
    dataloader,
    split_name,
    output_path,
    is_vae,
    device
):
    """
    Extract embeddings from a dataloader and save to CSV.

    Args:
        model: Trained autoencoder model
        dataloader: DataLoader for the split
        split_name: Name of the split (train/val/test)
        output_path: Path to save CSV file
        is_vae: Whether model is a VAE
        device: Device for computation
    """
    model.eval()

    all_embeddings = []
    all_paths = []
    all_genotypes = []

    print(f"Extracting embeddings for {split_name} split...")

    with torch.no_grad():
        for batch in tqdm(dataloader, desc=f"Processing {split_name}"):
            images = batch['image'].to(device)

            # Get latent representation
            if is_vae:
                mu, logvar = model.encode(images)
                # Use mean for deterministic embeddings
                z = mu
            else:
                z = model.encode(images)

            # Store results
            all_embeddings.append(z.cpu().numpy())
            all_paths.extend(batch['image_path'])
            all_genotypes.extend(batch['genotype'])

    # Concatenate all embeddings
    import numpy as np
    embeddings = np.concatenate(all_embeddings, axis=0)

    # Create DataFrame
    data = {
        'image_path': all_paths,
        'genotype': all_genotypes
    }

    # Add latent dimensions as columns
    n_dims = embeddings.shape[1]
    for i in range(n_dims):
        data[f'latent_dim_{i}'] = embeddings[:, i]

    df = pd.DataFrame(data)

    # Save to CSV
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} embeddings to: {output_path}")

    return df


def main():
    parser = argparse.ArgumentParser(description="Extract embeddings from trained autoencoder")
    parser.add_argument('--model_dir', type=str, required=True,
                       help='Directory containing trained model')
    parser.add_argument('--output_dir', type=str, default=None,
                       help='Directory to save embeddings (default: same as model_dir)')
    parser.add_argument('--splits', type=str, nargs='+', default=['train', 'val', 'test'],
                       choices=['train', 'val', 'test'],
                       help='Splits to extract embeddings for')
    parser.add_argument('--device', type=str, default='cuda' if torch.cuda.is_available() else 'cpu',
                       choices=['cuda', 'cpu'],
                       help='Device for computation')

    args = parser.parse_args()

    model_dir = Path(args.model_dir)
    output_dir = Path(args.output_dir) if args.output_dir else model_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device(args.device)
    print(f"Using device: {device}")

    # Load config
    config_path = model_dir / 'config.yaml'
    print(f"Loading config from: {config_path}")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Load model
    print("Initializing model...")
    model = ConvolutionalAutoencoder(config)

    # Load checkpoint
    checkpoint_path = model_dir / 'checkpoints' / 'best_model.pt'
    print(f"Loading checkpoint from: {checkpoint_path}")
    checkpoint = torch.load(checkpoint_path, map_location=device)

    if 'model_state_dict' in checkpoint:
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        model.load_state_dict(checkpoint)

    model = model.to(device)
    model.eval()

    is_vae = config['model'].get('is_vae', False)
    print(f"Model type: {'VAE' if is_vae else 'Standard AE'}")

    # Create dataloaders
    print("\nCreating dataloaders...")
    train_loader, val_loader, test_loader = create_dataloaders_from_config(config)

    dataloaders = {
        'train': train_loader,
        'val': val_loader,
        'test': test_loader
    }

    # Extract embeddings for requested splits
    print("\nExtracting embeddings...")
    for split in args.splits:
        if split not in dataloaders:
            print(f"Warning: Split '{split}' not available")
            continue

        dataloader = dataloaders[split]
        output_path = output_dir / f'embeddings_{split}.csv'

        extract_and_save_embeddings(
            model=model,
            dataloader=dataloader,
            split_name=split,
            output_path=output_path,
            is_vae=is_vae,
            device=device
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
