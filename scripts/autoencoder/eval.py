"""
Comprehensive Model Evaluation and Comparison Script for Autoencoders.

This script performs thorough evaluation and comparison of trained autoencoder models,
generating detailed visualizations, metrics, and interactive reports.

Usage:
    # Evaluate all models on all splits
    python src/autoencoder/eval.py --splits train val test --models_dir models/

    # Evaluate specific models
    python src/autoencoder/eval.py --splits test --model_paths models/model1/ models/model2/

    # Generate specific comparison types
    python src/autoencoder/eval.py --splits val --comparison_types reconstruction latent training

    # Custom output directory
    python src/autoencoder/eval.py --output_dir custom_eval/
"""

import os
import sys
import argparse
import json
import warnings
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
import traceback

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import numpy as np
import pandas as pd
from tqdm import tqdm
import yaml

# Visualization libraries
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

# Dimensionality reduction
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
try:
    import umap
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False
    warnings.warn("UMAP not available. Install with: pip install umap-learn")

# Clustering metrics
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score

# Image quality metrics
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import peak_signal_noise_ratio as psnr

# HTML generation
from jinja2 import Template
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent))

from model import ConvolutionalAutoencoder
from loss import create_loss_from_config, VAELoss, L1ReconstructionLoss, DiseaseWeightedL1Loss
from data_utils import create_dataloaders_from_config

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10


# ============================================================================
# Model Discovery and Loading
# ============================================================================

def discover_models(models_dir: str) -> List[Path]:
    """
    Find all trained model directories with saved checkpoints.

    Args:
        models_dir: Base directory containing model subdirectories

    Returns:
        List of paths to valid model directories
    """
    models_dir = Path(models_dir)
    model_dirs = []

    if not models_dir.exists():
        print(f"Warning: Models directory does not exist: {models_dir}")
        return model_dirs

    for path in models_dir.iterdir():
        if path.is_dir():
            # Check for required files
            checkpoint_path = path / 'checkpoints' / 'best_model.pt'
            config_path = path / 'config.yaml'

            if checkpoint_path.exists() and config_path.exists():
                model_dirs.append(path)
            elif any(path.glob('*.pt')):  # Fallback: any .pt file
                print(f"Warning: Model directory {path.name} has .pt files but missing standard structure")

    return sorted(model_dirs)


def load_model_info(model_dir: Path, device: torch.device) -> Optional[Dict[str, Any]]:
    """
    Load model, configuration, and training information.

    Args:
        model_dir: Path to model directory
        device: Device to load model onto

    Returns:
        Dictionary containing model information, or None if loading fails
    """
    try:
        # Load config
        config_path = model_dir / 'config.yaml'
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        # Initialize model
        model = ConvolutionalAutoencoder(config)

        # Load checkpoint
        checkpoint_path = model_dir / 'checkpoints' / 'best_model.pt'
        checkpoint = torch.load(checkpoint_path, map_location=device)

        # Load model state
        if 'model_state_dict' in checkpoint:
            model.load_state_dict(checkpoint['model_state_dict'])
        else:
            model.load_state_dict(checkpoint)

        model = model.to(device)
        model.eval()

        # Load training log
        training_log = None
        log_path = model_dir / 'logs' / 'training_log.csv'
        if log_path.exists():
            training_log = pd.read_csv(log_path)

        # Load training summary
        summary = None
        summary_path = model_dir / 'training_summary.json'
        if summary_path.exists():
            with open(summary_path, 'r') as f:
                summary = json.load(f)

        # Create shortened name for display
        name = model_dir.name
        # Extract key parameters for shorter name
        short_name = _create_short_name(name, config)

        return {
            'model': model,
            'config': config,
            'checkpoint': checkpoint,
            'training_log': training_log,
            'summary': summary,
            'name': name,
            'short_name': short_name,
            'path': model_dir
        }

    except Exception as e:
        print(f"Error loading model from {model_dir}: {e}")
        traceback.print_exc()
        return None


def _create_short_name(full_name: str, config: Dict[str, Any]) -> str:
    """
    Create a shortened, readable model name.

    Args:
        full_name: Full model directory name
        config: Model configuration

    Returns:
        Shortened name highlighting key parameters
    """
    # Extract key info
    model_type = "VAE" if config['model'].get('is_vae', False) else "AE"
    loss_type = config['loss'].get('type', 'l1')
    lr = config['training'].get('learning_rate', 0.001)
    bottleneck = config['model'].get('bottleneck_dim', 128)
    attention = config['model'].get('use_attention', False)

    parts = [model_type]

    if loss_type != 'l1':
        if loss_type == 'disease_weighted_l1':
            parts.append('DW')
        else:
            parts.append(loss_type.upper())

    parts.append(f"z{bottleneck}")
    parts.append(f"lr{lr}")

    if attention:
        parts.append("ATT")

    return "_".join(parts)


# ============================================================================
# Metric Computation Utilities
# ============================================================================

def compute_image_metrics(
    original: np.ndarray,
    reconstruction: np.ndarray,
    mask: Optional[np.ndarray] = None,
    data_range: float = 1.0
) -> Dict[str, float]:
    """
    Compute comprehensive image quality metrics.

    Args:
        original: Original image (H, W, C) in [0, 1]
        reconstruction: Reconstructed image (H, W, C) in [0, 1]
        mask: Optional binary mask (H, W) or (H, W, C)
        data_range: Data range for PSNR/SSIM computation

    Returns:
        Dictionary of metrics
    """
    metrics = {}

    # Apply mask if provided
    if mask is not None:
        # Ensure mask is 2D
        if mask.ndim == 3:
            mask = mask.mean(axis=2) > 0.5
        else:
            mask = mask > 0.5

        # Extract masked pixels
        orig_masked = original[mask]
        recon_masked = reconstruction[mask]
    else:
        orig_masked = original
        recon_masked = reconstruction

    # MSE (Mean Squared Error)
    mse = np.mean((orig_masked - recon_masked) ** 2)
    metrics['mse'] = float(mse)

    # MAE (Mean Absolute Error)
    mae = np.mean(np.abs(orig_masked - recon_masked))
    metrics['mae'] = float(mae)

    # RMSE (Root Mean Squared Error)
    rmse = np.sqrt(mse)
    metrics['rmse'] = float(rmse)

    # NRMSE (Normalized RMSE)
    data_range_actual = orig_masked.max() - orig_masked.min()
    if data_range_actual > 0:
        nrmse = rmse / data_range_actual
    else:
        nrmse = 0.0
    metrics['nrmse'] = float(nrmse)

    # SSIM (Structural Similarity Index)
    if mask is None:
        # Compute SSIM on full image
        ssim_value = ssim(
            original, reconstruction,
            data_range=data_range,
            channel_axis=2,
            win_size=7
        )
    else:
        # For masked regions, compute per-channel on masked image
        orig_masked_img = original * mask[..., None]
        recon_masked_img = reconstruction * mask[..., None]
        ssim_value = ssim(
            orig_masked_img, recon_masked_img,
            data_range=data_range,
            channel_axis=2,
            win_size=7
        )
    metrics['ssim'] = float(ssim_value)

    # PSNR (Peak Signal-to-Noise Ratio)
    if mse > 0:
        psnr_value = psnr(original, reconstruction, data_range=data_range)
    else:
        psnr_value = float('inf')
    metrics['psnr'] = float(psnr_value)

    # Per-channel MSE
    if original.ndim == 3:
        for c in range(original.shape[2]):
            if mask is not None:
                orig_ch = original[:, :, c][mask]
                recon_ch = reconstruction[:, :, c][mask]
            else:
                orig_ch = original[:, :, c]
                recon_ch = reconstruction[:, :, c]

            ch_mse = np.mean((orig_ch - recon_ch) ** 2)
            metrics[f'mse_ch{c}'] = float(ch_mse)

    return metrics


def apply_mask_from_path(image: np.ndarray, mask_path: str) -> np.ndarray:
    """
    Load and apply mask to an image.

    Args:
        image: Image array (H, W, C)
        mask_path: Path to mask file

    Returns:
        Masked image
    """
    import cv2

    if not os.path.exists(mask_path):
        return image

    mask = cv2.imread(mask_path, cv2.IMREAD_COLOR)
    if mask is None:
        return image

    mask = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)

    # Resize if needed
    if mask.shape[:2] != image.shape[:2]:
        mask = cv2.resize(mask, (image.shape[1], image.shape[0]),
                         interpolation=cv2.INTER_NEAREST)

    # Convert to binary mask
    mask_binary = (mask / 255.0).astype(np.float32)

    # Apply mask
    return image * mask_binary


# ============================================================================
# ModelEvaluator Class
# ============================================================================

class ModelEvaluator:
    """
    Evaluator for a single autoencoder model.

    Computes reconstruction metrics, extracts latent representations,
    and identifies best/worst samples.
    """

    def __init__(self, model_info: Dict[str, Any], device: torch.device):
        """
        Initialize evaluator.

        Args:
            model_info: Dictionary with model, config, etc.
            device: Device for computation
        """
        self.model_info = model_info
        self.model = model_info['model']
        self.config = model_info['config']
        self.device = device
        self.is_vae = self.config['model'].get('is_vae', False)

        # Initialize loss function
        self.loss_fn = create_loss_from_config(self.config)
        self.loss_fn = self.loss_fn.to(device)

        # Cache for computed results
        self.cache = {}

    def compute_reconstruction_metrics(
        self,
        dataloader: DataLoader,
        use_mask: bool = False,
        max_samples: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Compute reconstruction metrics for all samples in dataloader.

        Args:
            dataloader: DataLoader with samples to evaluate
            use_mask: Whether to compute masked metrics
            max_samples: Maximum number of samples to evaluate (None = all)

        Returns:
            DataFrame with per-sample metrics
        """
        self.model.eval()

        results = []
        sample_count = 0

        with torch.no_grad():
            for batch in tqdm(dataloader, desc="Computing metrics"):
                images = batch['image'].to(self.device)
                paths = batch['image_path']

                # Forward pass
                if self.is_vae:
                    reconstructions, mu, logvar = self.model(images)

                    # Compute KL divergence for each sample
                    kl_div = -0.5 * torch.sum(
                        1 + logvar - mu.pow(2) - logvar.exp(),
                        dim=1
                    ).cpu().numpy()
                else:
                    reconstructions = self.model(images)
                    kl_div = None

                # Convert to numpy
                orig_np = images.cpu().numpy().transpose(0, 2, 3, 1)  # (B, H, W, C)
                recon_np = reconstructions.cpu().numpy().transpose(0, 2, 3, 1)

                # Compute metrics for each sample
                for i in range(len(images)):
                    orig = orig_np[i]
                    recon = recon_np[i]

                    # Compute metrics (unmasked)
                    metrics = compute_image_metrics(orig, recon, mask=None)
                    metrics['image_path'] = paths[i]
                    metrics['masked'] = False

                    if self.is_vae and kl_div is not None:
                        metrics['kl_divergence'] = float(kl_div[i])

                    results.append(metrics)

                    # Compute masked metrics if requested
                    if use_mask:
                        # TODO: Load mask from dataset or mask_dir
                        # For now, skip masked metrics
                        pass

                    sample_count += 1
                    if max_samples and sample_count >= max_samples:
                        break

                if max_samples and sample_count >= max_samples:
                    break

        return pd.DataFrame(results)

    def extract_latent_representations(
        self,
        dataloader: DataLoader,
        deterministic: bool = True,
        max_samples: Optional[int] = None
    ) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        Extract latent representations for all samples.

        Args:
            dataloader: DataLoader with samples
            deterministic: For VAE, use mean instead of sampling
            max_samples: Maximum number of samples (None = all)

        Returns:
            Tuple of (latent_vectors, image_paths, genotypes)
        """
        self.model.eval()

        latents = []
        paths = []
        genotypes = []
        sample_count = 0

        with torch.no_grad():
            for batch in tqdm(dataloader, desc="Extracting latents"):
                images = batch['image'].to(self.device)

                # Get latent representation
                if self.is_vae:
                    mu, logvar = self.model.encode(images)
                    if deterministic:
                        z = mu
                    else:
                        z = self.model.reparameterize(mu, logvar)
                else:
                    z = self.model.encode(images)

                latents.append(z.cpu().numpy())
                paths.extend(batch['image_path'])
                genotypes.extend(batch['genotype'])

                sample_count += len(images)
                if max_samples and sample_count >= max_samples:
                    break

        latents = np.concatenate(latents, axis=0)

        if max_samples:
            latents = latents[:max_samples]
            paths = paths[:max_samples]
            genotypes = genotypes[:max_samples]

        return latents, paths, genotypes

    def save_latent_representations_to_csv(
        self,
        dataloader: DataLoader,
        output_path: str,
        deterministic: bool = True,
        max_samples: Optional[int] = None
    ):
        """
        Extract and save latent representations to CSV file.

        Args:
            dataloader: DataLoader with samples
            output_path: Path to save CSV file
            deterministic: For VAE, use mean instead of sampling
            max_samples: Maximum number of samples (None = all)
        """
        # Extract latent representations
        latents, paths, genotypes = self.extract_latent_representations(
            dataloader, deterministic=deterministic, max_samples=max_samples
        )

        # Create DataFrame
        data = {
            'image_path': paths,
            'genotype': genotypes
        }

        # Add latent dimensions as columns
        n_latent_dims = latents.shape[1]
        for i in range(n_latent_dims):
            data[f'latent_dim_{i}'] = latents[:, i]

        df = pd.DataFrame(data)

        # Save to CSV
        df.to_csv(output_path, index=False)
        print(f"Saved {len(df)} latent representations to: {output_path}")

        return df

    def get_best_worst_samples(
        self,
        dataloader: DataLoader,
        n: int = 5,
        metric: str = 'mse'
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Find best and worst reconstruction samples.

        Args:
            dataloader: DataLoader with samples
            n: Number of best/worst samples
            metric: Metric to use for ranking ('mse', 'ssim', etc.)

        Returns:
            Tuple of (best_samples, worst_samples) lists
        """
        # Compute metrics for all samples
        metrics_df = self.compute_reconstruction_metrics(dataloader, use_mask=False)

        # Sort by metric
        ascending = metric not in ['ssim', 'psnr']  # Higher is better for SSIM/PSNR
        sorted_df = metrics_df.sort_values(by=metric, ascending=ascending)

        # Get best and worst
        best_df = sorted_df.head(n)
        worst_df = sorted_df.tail(n)

        # Load images and reconstructions for these samples
        best_samples = self._load_samples_with_reconstructions(best_df, dataloader)
        worst_samples = self._load_samples_with_reconstructions(worst_df, dataloader)

        return best_samples, worst_samples

    def _load_samples_with_reconstructions(
        self,
        metrics_df: pd.DataFrame,
        dataloader: DataLoader
    ) -> List[Dict]:
        """
        Load images and reconstructions for specified samples.

        Args:
            metrics_df: DataFrame with image_path and metrics
            dataloader: DataLoader to get images from

        Returns:
            List of dictionaries with images, reconstructions, and metrics
        """
        samples = []
        paths_to_load = set(metrics_df['image_path'].tolist())

        self.model.eval()
        with torch.no_grad():
            for batch in dataloader:
                batch_paths = batch['image_path']

                # Check if any paths in this batch are needed
                indices = [i for i, p in enumerate(batch_paths) if p in paths_to_load]

                if not indices:
                    continue

                images = batch['image'][indices].to(self.device)

                # Get reconstructions
                if self.is_vae:
                    reconstructions, mu, logvar = self.model(images)
                else:
                    reconstructions = self.model(images)

                # Store samples
                for i, idx in enumerate(indices):
                    path = batch_paths[idx]

                    # Get metrics for this sample
                    sample_metrics = metrics_df[metrics_df['image_path'] == path].iloc[0].to_dict()

                    sample = {
                        'image': images[i].cpu().numpy(),
                        'reconstruction': reconstructions[i].cpu().numpy(),
                        'path': path,
                        'metrics': sample_metrics
                    }

                    samples.append(sample)
                    paths_to_load.remove(path)

                if not paths_to_load:
                    break

        return samples


# ============================================================================
# Multi-Model Comparator Class
# ============================================================================

class MultiModelComparator:
    """
    Compares multiple trained autoencoder models.

    Generates comparative visualizations and metrics across models.
    """

    def __init__(
        self,
        model_evaluators: Dict[str, ModelEvaluator],
        output_dir: Path,
        device: torch.device
    ):
        """
        Initialize comparator.

        Args:
            model_evaluators: Dict mapping model names to ModelEvaluator instances
            output_dir: Directory to save outputs
            device: Device for computation
        """
        self.evaluators = model_evaluators
        self.output_dir = output_dir
        self.device = device

        # Create output subdirectories
        self.dirs = {
            'reconstruction': output_dir / 'reconstruction_quality',
            'latent': output_dir / 'latent_space',
            'training': output_dir / 'training_analysis',
        }

        for d in self.dirs.values():
            d.mkdir(parents=True, exist_ok=True)

        # Subdirectories
        (self.dirs['reconstruction'] / 'random_samples').mkdir(exist_ok=True)
        (self.dirs['reconstruction'] / 'best_worst').mkdir(exist_ok=True)
        (self.dirs['reconstruction'] / 'error_maps').mkdir(exist_ok=True)
        (self.dirs['reconstruction'] / 'distributions').mkdir(exist_ok=True)

        (self.dirs['latent'] / 'tsne').mkdir(exist_ok=True)
        (self.dirs['latent'] / 'pca').mkdir(exist_ok=True)
        (self.dirs['latent'] / 'comparison').mkdir(exist_ok=True)
        if UMAP_AVAILABLE:
            (self.dirs['latent'] / 'umap').mkdir(exist_ok=True)

        (self.dirs['training'] / 'individual_curves').mkdir(exist_ok=True)

    def compare_reconstruction_quality(
        self,
        model_dataloaders: Dict[str, Dict[str, DataLoader]],
        splits: List[str],
        n_random_samples: int = 12,
        n_best_worst: int = 5
    ):
        """
        Generate reconstruction quality comparisons.

        Args:
            model_dataloaders: Dict mapping model names to dict of {split: DataLoader}
            splits: List of splits to evaluate
            n_random_samples: Number of random samples to visualize
            n_best_worst: Number of best/worst samples
        """
        print("\n" + "="*80)
        print("Reconstruction Quality Evaluation")
        print("="*80)

        # Collect metrics for all models and splits
        all_metrics = {}

        for split in splits:
            print(f"\nEvaluating split: {split}")

            for model_name, evaluator in self.evaluators.items():
                # Get this model's dataloader for this split
                model_loaders = model_dataloaders.get(model_name)
                if model_loaders is None:
                    print(f"  Warning: No dataloaders for model '{model_name}'")
                    continue

                dataloader = model_loaders.get(split)
                if dataloader is None:
                    print(f"  Warning: Model '{model_name}' has no dataloader for split '{split}'")
                    continue

                print(f"  Model: {model_name}")

                # Compute metrics
                metrics_df = evaluator.compute_reconstruction_metrics(
                    dataloader, use_mask=False
                )

                all_metrics[(model_name, split, False)] = metrics_df

                # Generate random samples visualization
                self._visualize_random_samples(
                    evaluator, dataloader, model_name, split, n_random_samples
                )

                # Get best/worst samples
                best_samples, worst_samples = evaluator.get_best_worst_samples(
                    dataloader, n=n_best_worst, metric='mse'
                )

                # Visualize best/worst
                self._visualize_best_worst(
                    best_samples, worst_samples, model_name, split
                )

        # Save metrics summary
        self._save_metrics_summary(all_metrics)

        # Generate distribution comparisons
        self._plot_metric_distributions(all_metrics, splits)

    def _visualize_random_samples(
        self,
        evaluator: ModelEvaluator,
        dataloader: DataLoader,
        model_name: str,
        split: str,
        n_samples: int = 12
    ):
        """
        Visualize random reconstruction samples.

        Args:
            evaluator: ModelEvaluator instance
            dataloader: DataLoader
            model_name: Name of model
            split: Data split name
            n_samples: Number of samples
        """
        evaluator.model.eval()

        # Get random batch
        for batch in dataloader:
            images = batch['image'][:n_samples].to(evaluator.device)
            break

        with torch.no_grad():
            if evaluator.is_vae:
                reconstructions, _, _ = evaluator.model(images)
            else:
                reconstructions = evaluator.model(images)

        # Convert to numpy
        images_np = images.cpu().numpy().transpose(0, 2, 3, 1)
        recons_np = reconstructions.cpu().numpy().transpose(0, 2, 3, 1)

        # Create visualization
        fig = plt.figure(figsize=(20, 6))
        gs = GridSpec(3, n_samples, figure=fig, hspace=0.3, wspace=0.1)

        for i in range(n_samples):
            # Original
            ax = fig.add_subplot(gs[0, i])
            ax.imshow(np.clip(images_np[i], 0, 1))
            ax.axis('off')
            if i == 0:
                ax.set_title('Original', fontsize=10, loc='left')

            # Reconstruction
            ax = fig.add_subplot(gs[1, i])
            ax.imshow(np.clip(recons_np[i], 0, 1))
            ax.axis('off')
            if i == 0:
                ax.set_title('Reconstruction', fontsize=10, loc='left')

            # Error map
            ax = fig.add_subplot(gs[2, i])
            error = np.abs(images_np[i] - recons_np[i]).mean(axis=2)
            im = ax.imshow(error, cmap='hot', vmin=0, vmax=0.3)
            ax.axis('off')
            if i == 0:
                ax.set_title('Abs Error', fontsize=10, loc='left')

            # Compute metrics
            metrics = compute_image_metrics(images_np[i], recons_np[i])
            ax.text(
                0.5, -0.1, f"SSIM: {metrics['ssim']:.3f}\nPSNR: {metrics['psnr']:.1f}",
                transform=ax.transAxes, ha='center', fontsize=7
            )

        # Add colorbar
        cbar_ax = fig.add_axes([0.92, 0.11, 0.01, 0.22])
        fig.colorbar(im, cax=cbar_ax)

        plt.suptitle(f"{evaluator.model_info['short_name']} - {split.upper()} Split",
                    fontsize=12, y=0.98)

        # Save
        output_path = (self.dirs['reconstruction'] / 'random_samples' /
                      f"{model_name}_samples_{split}.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def _visualize_best_worst(
        self,
        best_samples: List[Dict],
        worst_samples: List[Dict],
        model_name: str,
        split: str
    ):
        """
        Visualize best and worst reconstruction samples.

        Args:
            best_samples: List of best samples with images and metrics
            worst_samples: List of worst samples
            model_name: Model name
            split: Split name
        """
        for sample_type, samples in [('best', best_samples), ('worst', worst_samples)]:
            n = len(samples)
            if n == 0:
                continue

            fig, axes = plt.subplots(3, n, figsize=(3*n, 9))
            if n == 1:
                axes = axes[:, np.newaxis]

            for i, sample in enumerate(samples):
                # Get images (C, H, W) -> (H, W, C)
                img = sample['image'].transpose(1, 2, 0)
                recon = sample['reconstruction'].transpose(1, 2, 0)
                error = np.abs(img - recon).mean(axis=2)

                # Original
                axes[0, i].imshow(np.clip(img, 0, 1))
                axes[0, i].axis('off')
                axes[0, i].set_title('Original', fontsize=10)

                # Reconstruction
                axes[1, i].imshow(np.clip(recon, 0, 1))
                axes[1, i].axis('off')
                axes[1, i].set_title('Reconstruction', fontsize=10)

                # Error map
                im = axes[2, i].imshow(error, cmap='hot', vmin=0, vmax=0.3)
                axes[2, i].axis('off')
                axes[2, i].set_title('Error', fontsize=10)

                # Add metrics text
                metrics = sample['metrics']
                metrics_text = (
                    f"MSE: {metrics['mse']:.6f}\n"
                    f"SSIM: {metrics['ssim']:.4f}\n"
                    f"PSNR: {metrics['psnr']:.2f}"
                )
                axes[2, i].text(
                    0.5, -0.15, metrics_text,
                    transform=axes[2, i].transAxes,
                    ha='center', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
                )

            plt.suptitle(f"{model_name} - {sample_type.upper()} {n} ({split})",
                        fontsize=12)
            plt.tight_layout()

            # Save
            output_path = (self.dirs['reconstruction'] / 'best_worst' /
                          f"{model_name}_{sample_type}{n}_{split}.png")
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()

    def _save_metrics_summary(self, all_metrics: Dict[Tuple, pd.DataFrame]):
        """
        Save comprehensive metrics summary CSV.

        Args:
            all_metrics: Dict mapping (model_name, split, masked) to metrics DataFrame
        """
        summary_rows = []

        for (model_name, split, masked), df in all_metrics.items():
            # Compute summary statistics for each metric
            metric_cols = [c for c in df.columns if c not in ['image_path', 'masked']]

            for metric in metric_cols:
                values = df[metric].dropna()
                if len(values) == 0:
                    continue

                row = {
                    'model_name': model_name,
                    'split': split,
                    'masked': masked,
                    'metric': metric,
                    'mean': values.mean(),
                    'std': values.std(),
                    'median': values.median(),
                    'min': values.min(),
                    'max': values.max(),
                    'q25': values.quantile(0.25),
                    'q75': values.quantile(0.75)
                }
                summary_rows.append(row)

        summary_df = pd.DataFrame(summary_rows)

        # Save
        output_path = self.dirs['reconstruction'] / 'metrics_summary.csv'
        summary_df.to_csv(output_path, index=False)
        print(f"\nSaved metrics summary to: {output_path}")

    def _plot_metric_distributions(
        self,
        all_metrics: Dict[Tuple, pd.DataFrame],
        splits: List[str]
    ):
        """
        Plot metric distributions across models.

        Args:
            all_metrics: Metrics data
            splits: List of splits
        """
        # Metrics to plot
        metrics_to_plot = ['mse', 'mae', 'ssim', 'psnr']

        for split in splits:
            # Filter metrics for this split
            split_data = {
                model: df for (model, s, masked), df in all_metrics.items()
                if s == split and not masked
            }

            if not split_data:
                continue

            # Create plots
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            axes = axes.flatten()

            for idx, metric in enumerate(metrics_to_plot):
                ax = axes[idx]

                # Prepare data for violin plot
                data_list = []
                labels = []

                for model_name, df in split_data.items():
                    if metric in df.columns:
                        data_list.append(df[metric].dropna().values)
                        labels.append(model_name)

                if not data_list:
                    continue

                # Violin plot
                parts = ax.violinplot(
                    data_list, positions=range(len(data_list)),
                    showmeans=True, showmedians=True
                )

                ax.set_xticks(range(len(labels)))
                ax.set_xticklabels(labels, rotation=45, ha='right')
                ax.set_ylabel(metric.upper())
                ax.set_title(f'{metric.upper()} Distribution')
                ax.grid(axis='y', alpha=0.3)

            plt.suptitle(f'Metric Distributions - {split.upper()} Split', fontsize=14)
            plt.tight_layout()

            # Save
            output_path = (self.dirs['reconstruction'] / 'distributions' /
                          f'distributions_{split}.png')
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()

    def compare_latent_spaces(
        self,
        model_dataloaders: Dict[str, Dict[str, DataLoader]],
        splits: List[str],
        methods: List[str] = ['tsne', 'umap', 'pca'],
        max_samples: int = 2000
    ):
        """
        Compare latent space representations across models.

        Args:
            model_dataloaders: Dict mapping model names to dict of {split: DataLoader}
            splits: List of splits to analyze
            methods: Dimensionality reduction methods to use
            max_samples: Maximum samples per split
        """
        print("\n" + "="*80)
        print("Latent Space Analysis")
        print("="*80)

        # Extract latent representations for all models and splits
        all_latents = {}

        for split in splits:
            print(f"\nProcessing split: {split}")

            for model_name, evaluator in self.evaluators.items():
                # Check if embeddings already exist in the model directory
                model_dir = evaluator.model_info['path']
                existing_embeddings_path = model_dir / f'embeddings_{split}.csv'

                if existing_embeddings_path.exists():
                    print(f"  Loading existing embeddings for: {model_name} ({split})")
                    print(f"    From: {existing_embeddings_path}")

                    # Load embeddings from CSV
                    df = pd.read_csv(existing_embeddings_path)

                    # Extract data
                    paths = df['image_path'].tolist()
                    genotypes = df['genotype'].tolist()

                    # Extract latent dimensions
                    latent_cols = [col for col in df.columns if col.startswith('latent_dim_')]
                    latents = df[latent_cols].values

                    all_latents[(model_name, split)] = {
                        'latents': latents,
                        'paths': paths,
                        'genotypes': genotypes
                    }

                    print(f"    Loaded {len(df)} embeddings from existing file")

                else:
                    # Embeddings don't exist, extract them
                    print(f"  Extracting latents for: {model_name} ({split})")

                    # Get this model's dataloader for this split
                    model_loaders = model_dataloaders.get(model_name)
                    if model_loaders is None:
                        print(f"  Warning: No dataloaders for model '{model_name}'")
                        continue

                    dataloader = model_loaders.get(split)
                    if dataloader is None:
                        print(f"  Warning: Model '{model_name}' has no dataloader for split '{split}'")
                        continue

                    latents, paths, genotypes = evaluator.extract_latent_representations(
                        dataloader, deterministic=True, max_samples=max_samples
                    )

                    all_latents[(model_name, split)] = {
                        'latents': latents,
                        'paths': paths,
                        'genotypes': genotypes
                    }

                # Save latent representations to evaluation output directory
                csv_filename = f"{model_name}_latents_{split}.csv"
                csv_path = self.dirs['latent'] / csv_filename

                # Create DataFrame and save
                data = {
                    'image_path': paths,
                    'genotype': genotypes
                }
                for i in range(latents.shape[1]):
                    data[f'latent_dim_{i}'] = latents[:, i]

                df = pd.DataFrame(data)
                df.to_csv(csv_path, index=False)
                print(f"    Saved {len(df)} latent representations to: {csv_filename}")

        # Apply dimensionality reduction and visualize
        for method in methods:
            if method == 'umap' and not UMAP_AVAILABLE:
                print(f"Skipping {method} (not installed)")
                continue

            print(f"\nApplying {method.upper()} dimensionality reduction...")
            self._apply_dim_reduction_and_plot(all_latents, method, splits)

        # Compute clustering metrics
        self._compute_clustering_metrics(all_latents, splits)

    def _apply_dim_reduction_and_plot(
        self,
        all_latents: Dict,
        method: str,
        splits: List[str]
    ):
        """
        Apply dimensionality reduction and create visualizations.

        Args:
            all_latents: Dictionary of latent representations
            method: Reduction method ('tsne', 'umap', 'pca')
            splits: List of splits
        """
        for split in splits:
            # Individual model plots
            for model_name, evaluator in self.evaluators.items():
                key = (model_name, split)
                if key not in all_latents:
                    continue

                data = all_latents[key]
                latents = data['latents']
                genotypes = data['genotypes']

                # Apply dimensionality reduction
                if method == 'tsne':
                    reducer = TSNE(n_components=2, random_state=42, perplexity=30)
                elif method == 'umap':
                    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15)
                elif method == 'pca':
                    reducer = PCA(n_components=2, random_state=42)
                else:
                    continue

                embedding = reducer.fit_transform(latents)

                # Create plot
                fig, ax = plt.subplots(figsize=(10, 8))

                # Check if this is a single-model evaluation
                is_single_model = len(self.evaluators) == 1

                if is_single_model:
                    # Single model: use single color for all points
                    ax.scatter(
                        embedding[:, 0], embedding[:, 1],
                        c='steelblue',
                        alpha=0.6,
                        s=20,
                        edgecolors='none'
                    )
                else:
                    # Multiple models: color by genotype for comparison
                    unique_genotypes = list(set(genotypes))
                    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_genotypes)))
                    genotype_to_color = dict(zip(unique_genotypes, colors))

                    for genotype in unique_genotypes:
                        mask = np.array(genotypes) == genotype
                        ax.scatter(
                            embedding[mask, 0], embedding[mask, 1],
                            c=[genotype_to_color[genotype]], label=genotype,
                            alpha=0.6, s=20
                        )

                    # Only show legend if not too many genotypes
                    if len(unique_genotypes) <= 20:
                        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

                ax.set_xlabel(f'{method.upper()} 1')
                ax.set_ylabel(f'{method.upper()} 2')
                ax.set_title(f'{evaluator.model_info["short_name"]} - {split.upper()} - {method.upper()}')

                plt.tight_layout()

                # Save
                output_path = (self.dirs['latent'] / method /
                              f"{model_name}_{method}_{split}.png")
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()

            # Cross-model comparison plot
            self._plot_cross_model_latent_comparison(all_latents, method, split)

    def _plot_cross_model_latent_comparison(
        self,
        all_latents: Dict,
        method: str,
        split: str
    ):
        """
        Plot latent spaces of all models on same plot.

        Args:
            all_latents: Latent data
            method: Reduction method
            split: Data split
        """
        fig, ax = plt.subplots(figsize=(12, 10))

        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
        colors = plt.cm.tab10(np.linspace(0, 1, len(self.evaluators)))

        for idx, (model_name, evaluator) in enumerate(self.evaluators.items()):
            key = (model_name, split)
            if key not in all_latents:
                continue

            data = all_latents[key]
            latents = data['latents']

            # Apply dimensionality reduction
            if method == 'tsne':
                reducer = TSNE(n_components=2, random_state=42, perplexity=30)
            elif method == 'umap':
                reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15)
            elif method == 'pca':
                reducer = PCA(n_components=2, random_state=42)
            else:
                continue

            embedding = reducer.fit_transform(latents)

            # Plot - color all points by model
            marker = markers[idx % len(markers)]
            color = colors[idx]

            ax.scatter(
                embedding[:, 0], embedding[:, 1],
                c=[color] * len(embedding),  # Color all points by model color
                marker=marker,
                label=evaluator.model_info['short_name'],
                alpha=0.5, s=30,
                edgecolors='none'
            )

        ax.set_xlabel(f'{method.upper()} 1')
        ax.set_ylabel(f'{method.upper()} 2')
        ax.set_title(f'All Models Comparison - {split.upper()} - {method.upper()}')
        ax.legend(loc='best')
        ax.grid(alpha=0.3)

        plt.tight_layout()

        # Save
        output_path = (self.dirs['latent'] / 'comparison' /
                      f"all_models_{method}_{split}.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def _compute_clustering_metrics(
        self,
        all_latents: Dict,
        splits: List[str]
    ):
        """
        Compute and save clustering quality metrics.

        Args:
            all_latents: Latent representations
            splits: Data splits
        """
        results = []

        for split in splits:
            for model_name, evaluator in self.evaluators.items():
                key = (model_name, split)
                if key not in all_latents:
                    continue

                data = all_latents[key]
                latents = data['latents']
                genotypes = data['genotypes']

                # Encode genotypes as integers
                unique_genotypes = list(set(genotypes))
                genotype_to_id = {g: i for i, g in enumerate(unique_genotypes)}
                labels = np.array([genotype_to_id[g] for g in genotypes])

                # Compute metrics if we have enough samples and labels
                if len(latents) < 10 or len(unique_genotypes) < 2:
                    continue

                try:
                    silhouette = silhouette_score(latents, labels)
                    davies_bouldin = davies_bouldin_score(latents, labels)
                    calinski = calinski_harabasz_score(latents, labels)

                    results.append({
                        'model_name': model_name,
                        'split': split,
                        'method': 'raw_latent',
                        'n_samples': len(latents),
                        'n_clusters': len(unique_genotypes),
                        'silhouette_score': silhouette,
                        'davies_bouldin_score': davies_bouldin,
                        'calinski_harabasz_score': calinski
                    })
                except Exception as e:
                    print(f"Warning: Could not compute clustering metrics for {model_name}/{split}: {e}")

        if results:
            df = pd.DataFrame(results)
            output_path = self.dirs['latent'] / 'clustering_metrics.csv'
            df.to_csv(output_path, index=False)
            print(f"\nSaved clustering metrics to: {output_path}")

    def compare_training_curves(self):
        """
        Compare training curves across all models.
        """
        print("\n" + "="*80)
        print("Training Analysis")
        print("="*80)

        # Collect training logs
        training_data = {}
        for model_name, evaluator in self.evaluators.items():
            if evaluator.model_info['training_log'] is not None:
                training_data[model_name] = evaluator.model_info['training_log']

        if not training_data:
            print("No training logs available")
            return

        # Plot overlay
        self._plot_training_curves_overlay(training_data)

        # Plot individual curves
        for model_name, log in training_data.items():
            self._plot_individual_training_curve(model_name, log)

        # Compare training summaries
        self._create_training_summary_comparison()

    def _plot_training_curves_overlay(self, training_data: Dict[str, pd.DataFrame]):
        """
        Plot training curves for all models overlaid.

        Args:
            training_data: Dict mapping model names to training log DataFrames
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        colors = plt.cm.tab10(np.linspace(0, 1, len(training_data)))

        for idx, (model_name, log) in enumerate(training_data.items()):
            color = colors[idx]
            label = self.evaluators[model_name].model_info['short_name']

            # Training loss
            if 'train_loss' in log.columns:
                axes[0].plot(log['epoch'], log['train_loss'], '-',
                           color=color, label=label, alpha=0.7)

            # Validation loss
            if 'val_loss' in log.columns:
                axes[1].plot(log['epoch'], log['val_loss'], '-',
                           color=color, label=label, alpha=0.7)

        axes[0].set_xlabel('Epoch')
        axes[0].set_ylabel('Training Loss')
        axes[0].set_title('Training Loss Comparison')
        axes[0].legend()
        axes[0].grid(alpha=0.3)

        axes[1].set_xlabel('Epoch')
        axes[1].set_ylabel('Validation Loss')
        axes[1].set_title('Validation Loss Comparison')
        axes[1].legend()
        axes[1].grid(alpha=0.3)

        plt.tight_layout()

        # Save
        output_path = self.dirs['training'] / 'loss_curves_comparison.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Saved training curves comparison to: {output_path}")

    def _plot_individual_training_curve(self, model_name: str, log: pd.DataFrame):
        """
        Plot detailed training curve for a single model.

        Args:
            model_name: Model name
            log: Training log DataFrame
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Loss curves
        ax = axes[0, 0]
        if 'train_loss' in log.columns and 'val_loss' in log.columns:
            ax.plot(log['epoch'], log['train_loss'], label='Train', alpha=0.7)
            ax.plot(log['epoch'], log['val_loss'], label='Val', alpha=0.7)
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Loss')
            ax.set_title('Loss Curves')
            ax.legend()
            ax.grid(alpha=0.3)

        # Learning rate
        ax = axes[0, 1]
        if 'learning_rate' in log.columns:
            ax.plot(log['epoch'], log['learning_rate'], color='orange')
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Learning Rate')
            ax.set_title('Learning Rate Schedule')
            ax.set_yscale('log')
            ax.grid(alpha=0.3)

        # Loss components (for VAE)
        ax = axes[1, 0]
        if 'recon_loss' in log.columns and 'kl_loss' in log.columns:
            ax.plot(log['epoch'], log['recon_loss'], label='Reconstruction', alpha=0.7)
            ax.plot(log['epoch'], log['kl_loss'], label='KL Divergence', alpha=0.7)
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Loss Component')
            ax.set_title('VAE Loss Components')
            ax.legend()
            ax.grid(alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'N/A (not VAE)', ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')

        # Train-Val gap
        ax = axes[1, 1]
        if 'train_loss' in log.columns and 'val_loss' in log.columns:
            gap = log['val_loss'] - log['train_loss']
            ax.plot(log['epoch'], gap, color='red', alpha=0.7)
            ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Val Loss - Train Loss')
            ax.set_title('Train-Val Gap (Overfitting Metric)')
            ax.grid(alpha=0.3)

        plt.suptitle(f"Training Analysis: {self.evaluators[model_name].model_info['short_name']}",
                    fontsize=14)
        plt.tight_layout()

        # Save
        output_path = (self.dirs['training'] / 'individual_curves' /
                      f"{model_name}_detailed.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def _create_training_summary_comparison(self):
        """
        Create CSV comparing training summaries.
        """
        summary_rows = []

        for model_name, evaluator in self.evaluators.items():
            summary = evaluator.model_info.get('summary')
            if summary is None:
                continue

            row = {
                'model_name': model_name,
                'short_name': evaluator.model_info['short_name'],
                **summary
            }
            summary_rows.append(row)

        if summary_rows:
            df = pd.DataFrame(summary_rows)
            output_path = self.dirs['training'] / 'training_summary_comparison.csv'
            df.to_csv(output_path, index=False)
            print(f"Saved training summary to: {output_path}")

    def generate_html_report(self):
        """
        Generate interactive HTML report.
        """
        print("\n" + "="*80)
        print("Generating HTML Report")
        print("="*80)

        # TODO: Implement comprehensive HTML report using Jinja2 and Plotly
        print("HTML report generation not yet implemented")
        pass


# ============================================================================
# Main Execution
# ============================================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Evaluate and compare trained autoencoder models'
    )

    parser.add_argument(
        '--models_dir',
        type=str,
        default='models/',
        help='Directory containing trained models'
    )

    parser.add_argument(
        '--model_paths',
        type=str,
        nargs='+',
        help='Specific model directories to evaluate'
    )

    parser.add_argument(
        '--splits',
        type=str,
        nargs='+',
        default=['val'],
        choices=['train', 'val', 'test'],
        help='Data splits to evaluate on'
    )

    parser.add_argument(
        '--comparison_types',
        type=str,
        nargs='+',
        default=['reconstruction', 'latent', 'training'],
        choices=['reconstruction', 'latent', 'training', 'all'],
        help='Types of comparisons to generate'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        help='Custom output directory (default: models/evaluations/comparison_TIMESTAMP/)'
    )

    parser.add_argument(
        '--n_random_samples',
        type=int,
        default=12,
        help='Number of random samples to visualize'
    )

    parser.add_argument(
        '--n_best_worst',
        type=int,
        default=5,
        help='Number of best/worst samples to show'
    )

    parser.add_argument(
        '--max_samples_latent',
        type=int,
        default=2000,
        help='Maximum samples for latent space analysis'
    )

    parser.add_argument(
        '--dim_reduction_methods',
        type=str,
        nargs='+',
        default=['tsne', 'pca'],
        choices=['tsne', 'umap', 'pca'],
        help='Dimensionality reduction methods to use'
    )

    parser.add_argument(
        '--device',
        type=str,
        default='cuda' if torch.cuda.is_available() else 'cpu',
        choices=['cuda', 'cpu'],
        help='Device for computation'
    )

    parser.add_argument(
        '--config',
        type=str,
        help='Path to data config file (for creating dataloaders)'
    )

    return parser.parse_args()


def main():
    """Main evaluation workflow."""
    args = parse_args()

    print("="*80)
    print("AUTOENCODER MODEL EVALUATION AND COMPARISON")
    print("="*80)

    # Set device
    device = torch.device(args.device)
    print(f"\nUsing device: {device}")

    # Discover models
    if args.model_paths:
        model_dirs = [Path(p) for p in args.model_paths]
    else:
        model_dirs = discover_models(args.models_dir)

    if not model_dirs:
        print(f"\nError: No models found!")
        if args.model_paths:
            print(f"Searched paths: {args.model_paths}")
        else:
            print(f"Searched directory: {args.models_dir}")
        sys.exit(1)

    print(f"\nFound {len(model_dirs)} model(s):")
    for d in model_dirs:
        print(f"  - {d.name}")

    # Load models
    print("\n" + "="*80)
    print("Loading Models")
    print("="*80)

    model_evaluators = {}
    for model_dir in model_dirs:
        print(f"\nLoading: {model_dir.name}")
        model_info = load_model_info(model_dir, device)

        if model_info is not None:
            evaluator = ModelEvaluator(model_info, device)
            model_evaluators[model_info['name']] = evaluator
            print(f"  ✓ Loaded successfully as '{model_info['short_name']}'")
        else:
            print(f"  ✗ Failed to load")

    if not model_evaluators:
        print("\nError: No models could be loaded!")
        sys.exit(1)

    print(f"\nSuccessfully loaded {len(model_evaluators)} model(s)")

    # Create output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path("models") / "evaluations" / f"comparison_{timestamp}"

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput directory: {output_dir}")

    # Create dataloaders for each model using its own config
    print("\n" + "="*80)
    print("Creating DataLoaders")
    print("="*80)

    model_dataloaders = {}  # Dict: model_name -> dataloaders dict

    for model_name, evaluator in model_evaluators.items():
        print(f"\nCreating dataloaders for: {model_name}")

        # Use model's own config
        config = evaluator.config

        # Override with custom config if provided
        if args.config:
            print(f"  Using custom config from: {args.config}")
            with open(args.config, 'r') as f:
                config = yaml.safe_load(f)

        try:
            train_loader, val_loader, test_loader = create_dataloaders_from_config(config)
            dataloaders = {
                'train': train_loader,
                'val': val_loader,
                'test': test_loader
            }
            model_dataloaders[model_name] = dataloaders
            print(f"  ✓ Created dataloaders for splits: {list(dataloaders.keys())}")
        except Exception as e:
            print(f"  ✗ Error creating dataloaders: {e}")
            traceback.print_exc()
            print(f"  Skipping model {model_name}")
            continue

    if not model_dataloaders:
        print("\nError: No dataloaders could be created!")
        sys.exit(1)

    # Create comparator
    comparator = MultiModelComparator(model_evaluators, output_dir, device)

    # Determine which comparisons to run
    comparison_types = args.comparison_types
    if 'all' in comparison_types:
        comparison_types = ['reconstruction', 'latent', 'training']

    # Run comparisons
    if 'reconstruction' in comparison_types:
        comparator.compare_reconstruction_quality(
            model_dataloaders,
            args.splits,
            n_random_samples=args.n_random_samples,
            n_best_worst=args.n_best_worst
        )

    if 'latent' in comparison_types:
        comparator.compare_latent_spaces(
            model_dataloaders,
            args.splits,
            methods=args.dim_reduction_methods,
            max_samples=args.max_samples_latent
        )

    if 'training' in comparison_types:
        comparator.compare_training_curves()

    # Generate HTML report
    # comparator.generate_html_report()

    print("\n" + "="*80)
    print("EVALUATION COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_dir}")
    print("\nGenerated outputs:")
    print(f"  - Reconstruction quality: {output_dir / 'reconstruction_quality'}")
    print(f"  - Latent space analysis: {output_dir / 'latent_space'}")
    print(f"  - Training analysis: {output_dir / 'training_analysis'}")


if __name__ == "__main__":
    main()

