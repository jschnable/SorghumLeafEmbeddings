"""
Loss functions for leaf image autoencoders and VAEs.

This module provides:
1. Standard L1 reconstruction loss
2. Disease-weighted L1 reconstruction loss (emphasizes non-green pixels)
3. KL divergence loss for VAE latent space regularization
4. Combined losses for VAE training

Mathematical Formulations:
--------------------------
1. L1 Loss: L = (1/N) * Σ|x - x̂|

2. Disease-Weighted L1: L = (1/N) * Σ[w(x) * |x - x̂|]
   where w(x) = 1 + α * distance_from_green(x)

3. KL Divergence: L_KL = -0.5 * Σ(1 + log(σ²) - μ² - σ²)
   where μ, σ are mean and std of learned latent distribution

4. VAE Loss: L_total = L_reconstruction + β * L_KL
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple, Optional, Dict, Any
import yaml


class L1ReconstructionLoss(nn.Module):
    """
    Standard L1 (Mean Absolute Error) reconstruction loss.

    Computes the element-wise absolute difference between original and
    reconstructed images, then averages across all elements.

    Formula: L = (1/N) * Σ|x - x̂|
    """

    def __init__(self):
        """Initialize L1 reconstruction loss."""
        super().__init__()

    def forward(
        self,
        reconstructed: torch.Tensor,
        original: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute L1 reconstruction loss.

        Args:
            reconstructed: Reconstructed images, shape (B, C, H, W)
            original: Original images, shape (B, C, H, W)

        Returns:
            Scalar loss value
        """
        return F.l1_loss(reconstructed, original, reduction='mean')


class DiseaseWeightedL1Loss(nn.Module):
    """
    L1 reconstruction loss weighted by deviation from green color.

    This loss emphasizes pixels that deviate from green, which typically
    indicates disease symptoms in leaf images. Healthy leaf tissue appears
    green, while diseased areas often show yellowing, browning, or necrosis.

    The weighting function increases the penalty for pixels that are less green,
    helping the autoencoder learn to better reconstruct disease features.

    Formula: L = (1/N) * Σ[w(x) * |x - x̂|]
    where w(x) = 1 + α * distance_from_green(x)

    Distance from green is computed in RGB space as:
    - greenness = G channel value
    - non_greenness = (R + B) / 2 - G
    - distance = max(0, non_greenness)
    """

    def __init__(
        self,
        alpha: float = 1.0,
        colorspace: str = 'RGB',
        epsilon: float = 1e-6
    ):
        """
        Initialize disease-weighted L1 loss.

        Args:
            alpha: Strength of disease weighting. Higher values increase the
                   penalty for non-green pixels. Default: 1.0
            colorspace: Color space of input images ('RGB' or 'LAB').
                       Default: 'RGB'
            epsilon: Small constant for numerical stability. Default: 1e-6
        """
        super().__init__()
        self.alpha = alpha
        self.colorspace = colorspace.upper()
        self.epsilon = epsilon

        if self.colorspace not in ['RGB', 'LAB']:
            raise ValueError(f"Unsupported colorspace: {colorspace}. Use 'RGB' or 'LAB'.")

    def _compute_greenness_weight_rgb(self, image: torch.Tensor) -> torch.Tensor:
        """
        Compute disease weighting based on deviation from green in RGB space.

        Green pixels have high G channel and lower R, B channels.
        Non-green pixels (disease symptoms) have higher R and/or B relative to G.

        Args:
            image: Input image tensor, shape (B, 3, H, W)

        Returns:
            Weight map, shape (B, 1, H, W), values >= 1.0
        """
        # Separate channels (assuming RGB order)
        r = image[:, 0:1, :, :]  # Red channel
        g = image[:, 1:2, :, :]  # Green channel
        b = image[:, 2:3, :, :]  # Blue channel

        # Compute non-greenness: average of R and B minus G
        # Positive values indicate deviation from green
        non_greenness = ((r + b) / 2.0) - g

        # Clamp to positive values (only penalize non-green, not extra-green)
        distance_from_green = torch.clamp(non_greenness, min=0.0)

        # Compute weight: 1.0 for green pixels, higher for non-green
        # Weight = 1 + alpha * distance_from_green
        weight = 1.0 + self.alpha * distance_from_green

        return weight

    def _compute_greenness_weight_lab(self, image: torch.Tensor) -> torch.Tensor:
        """
        Compute disease weighting based on deviation from green in LAB space.

        In LAB color space:
        - L: Lightness (0-100)
        - a: green-red axis (negative = green, positive = red)
        - b: blue-yellow axis (negative = blue, positive = yellow)

        Green corresponds to negative 'a' values.
        Disease symptoms (yellowing, browning) have positive 'a' and/or 'b'.

        Args:
            image: Input image tensor in LAB space, shape (B, 3, H, W)

        Returns:
            Weight map, shape (B, 1, H, W), values >= 1.0
        """
        # Assuming image is normalized to [0, 1]
        # L channel: image[:, 0:1, :, :]
        a = image[:, 1:2, :, :]  # a channel (green-red)
        b = image[:, 2:3, :, :]  # b channel (blue-yellow)

        # In normalized [0, 1] space, 0.5 is neutral
        # Values < 0.5 in 'a' channel indicate green
        # Values > 0.5 indicate red (disease)
        a_deviation = torch.clamp(a - 0.5, min=0.0)

        # Yellowing is indicated by positive b values
        b_deviation = torch.clamp(b - 0.5, min=0.0)

        # Combine deviations (Euclidean distance from green)
        distance_from_green = torch.sqrt(a_deviation**2 + b_deviation**2 + self.epsilon)

        # Compute weight
        weight = 1.0 + self.alpha * distance_from_green

        return weight

    def forward(
        self,
        reconstructed: torch.Tensor,
        original: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute disease-weighted L1 reconstruction loss.

        Args:
            reconstructed: Reconstructed images, shape (B, C, H, W)
            original: Original images, shape (B, C, H, W)

        Returns:
            Scalar loss value
        """
        # Compute base L1 loss (element-wise absolute difference)
        l1_loss = torch.abs(reconstructed - original)

        # Compute disease weights based on original image
        if self.colorspace == 'RGB':
            weights = self._compute_greenness_weight_rgb(original)
        else:  # LAB
            weights = self._compute_greenness_weight_lab(original)

        # Apply weights to L1 loss
        # Broadcast weights across channels if needed
        if weights.shape[1] == 1 and l1_loss.shape[1] > 1:
            weights = weights.expand_as(l1_loss)

        weighted_loss = weights * l1_loss

        # Return mean loss
        return weighted_loss.mean()


class KLDivergenceLoss(nn.Module):
    """
    KL divergence loss for VAE latent space regularization.

    Computes KL divergence between learned latent distribution q(z|x) and
    standard normal prior p(z) = N(0, I).

    For a diagonal Gaussian q(z|x) = N(μ, σ²I):
    KL(q||p) = -0.5 * Σ(1 + log(σ²) - μ² - σ²)

    This encourages the latent space to be continuous and well-structured,
    enabling smooth interpolation and sampling.
    """

    def __init__(self):
        """Initialize KL divergence loss."""
        super().__init__()

    def forward(
        self,
        mu: torch.Tensor,
        logvar: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute KL divergence loss.

        Args:
            mu: Mean of latent distribution, shape (B, latent_dim)
            logvar: Log variance of latent distribution, shape (B, latent_dim)
                   Using log(σ²) instead of σ² for numerical stability

        Returns:
            Scalar loss value (mean across batch and latent dimensions)
        """
        # KL divergence formula: -0.5 * Σ(1 + log(σ²) - μ² - σ²)
        # = -0.5 * Σ(1 + logvar - μ² - exp(logvar))
        kl_divergence = -0.5 * torch.sum(
            1 + logvar - mu.pow(2) - logvar.exp(),
            dim=1  # Sum over latent dimensions
        )

        # Return mean across batch
        return kl_divergence.mean()


class VAELoss(nn.Module):
    """
    Combined loss for Variational Autoencoder (VAE) training.

    Combines reconstruction loss with KL divergence loss using a beta weighting:
    L_total = L_reconstruction + β * L_KL

    The beta parameter controls the trade-off between reconstruction quality
    and latent space regularization (beta-VAE framework). Higher beta encourages
    more disentangled representations but may reduce reconstruction quality.

    Beta values:
    - β = 1.0: Standard VAE
    - β > 1.0: Beta-VAE with stronger regularization
    - β < 1.0: Emphasizes reconstruction over regularization
    """

    def __init__(
        self,
        reconstruction_loss: nn.Module,
        beta: float = 1.0
    ):
        """
        Initialize VAE loss.

        Args:
            reconstruction_loss: Reconstruction loss module (L1 or DiseaseWeightedL1)
            beta: Weight for KL divergence term. Default: 1.0 (standard VAE)
        """
        super().__init__()
        self.reconstruction_loss = reconstruction_loss
        self.kl_loss = KLDivergenceLoss()
        self.beta = beta

    def forward(
        self,
        reconstructed: torch.Tensor,
        original: torch.Tensor,
        mu: torch.Tensor,
        logvar: torch.Tensor,
        return_components: bool = False
    ) -> torch.Tensor | Dict[str, torch.Tensor]:
        """
        Compute VAE loss.

        Args:
            reconstructed: Reconstructed images, shape (B, C, H, W)
            original: Original images, shape (B, C, H, W)
            mu: Mean of latent distribution, shape (B, latent_dim)
            logvar: Log variance of latent distribution, shape (B, latent_dim)
            return_components: If True, return dict with loss components.
                              Default: False

        Returns:
            If return_components is False:
                Scalar total loss value
            If return_components is True:
                Dictionary with keys:
                    'total': Total loss
                    'reconstruction': Reconstruction loss
                    'kl': KL divergence loss
                    'beta': Beta value used
        """
        # Compute reconstruction loss
        recon_loss = self.reconstruction_loss(reconstructed, original)

        # Compute KL divergence loss
        kl_loss = self.kl_loss(mu, logvar)

        # Compute total loss
        total_loss = recon_loss + self.beta * kl_loss

        if return_components:
            return {
                'total': total_loss,
                'reconstruction': recon_loss,
                'kl': kl_loss,
                'beta': torch.tensor(self.beta)
            }

        return total_loss


def create_loss_from_config(
    config: Dict[str, Any] | str
) -> nn.Module:
    """
    Factory function to create loss function from configuration.

    Reads configuration and returns appropriate loss function for training.

    Args:
        config: Either a configuration dictionary or path to YAML config file

    Returns:
        Loss function module (nn.Module)

    Raises:
        ValueError: If configuration is invalid or loss type is unsupported

    Example config structure:
        loss:
            type: 'vae'  # Options: 'l1', 'disease_weighted_l1', 'vae'
            reconstruction_type: 'disease_weighted_l1'  # For VAE only
            disease_weight_alpha: 2.0  # For disease-weighted loss
            beta: 1.0  # For VAE
            colorspace: 'RGB'  # 'RGB' or 'LAB'
    """
    # Load config from file if path is provided
    if isinstance(config, str):
        with open(config, 'r') as f:
            config = yaml.safe_load(f)

    # Extract loss configuration
    loss_config = config.get('loss', {})
    loss_type = loss_config.get('type', 'l1').lower()
    colorspace = config.get('image', {}).get('colorspace', 'RGB')

    # Create reconstruction loss
    if loss_type == 'l1':
        return L1ReconstructionLoss()

    elif loss_type == 'disease_weighted_l1':
        alpha = loss_config.get('disease_weight_alpha', 1.0)
        return DiseaseWeightedL1Loss(
            alpha=alpha,
            colorspace=colorspace
        )

    elif loss_type == 'vae':
        # Determine reconstruction loss type
        recon_type = loss_config.get('reconstruction_type', 'l1').lower()

        if recon_type == 'l1':
            recon_loss = L1ReconstructionLoss()
        elif recon_type == 'disease_weighted_l1':
            alpha = loss_config.get('disease_weight_alpha', 1.0)
            recon_loss = DiseaseWeightedL1Loss(
                alpha=alpha,
                colorspace=colorspace
            )
        else:
            raise ValueError(
                f"Unsupported reconstruction type for VAE: {recon_type}. "
                "Use 'l1' or 'disease_weighted_l1'."
            )

        # Get beta parameter
        beta = loss_config.get('beta', 1.0)

        # Create VAE loss
        return VAELoss(
            reconstruction_loss=recon_loss,
            beta=beta
        )

    else:
        raise ValueError(
            f"Unsupported loss type: {loss_type}. "
            "Use 'l1', 'disease_weighted_l1', or 'vae'."
        )


# ============================================================================
# Example Usage
# ============================================================================

if __name__ == "__main__":
    """
    Example demonstrating how to use the loss functions.
    """
    import torch

    # Create dummy data
    batch_size = 8
    channels = 3
    height, width = 256, 256
    latent_dim = 128

    # Simulated original and reconstructed images
    original = torch.rand(batch_size, channels, height, width)
    reconstructed = torch.rand(batch_size, channels, height, width)

    # Simulated VAE latent parameters
    mu = torch.randn(batch_size, latent_dim)
    logvar = torch.randn(batch_size, latent_dim)

    print("=" * 80)
    print("Loss Function Examples")
    print("=" * 80)

    # Example 1: Standard L1 Loss
    print("\n1. Standard L1 Reconstruction Loss")
    print("-" * 80)
    l1_loss = L1ReconstructionLoss()
    loss_value = l1_loss(reconstructed, original)
    print(f"Loss value: {loss_value.item():.6f}")

    # Example 2: Disease-Weighted L1 Loss (RGB)
    print("\n2. Disease-Weighted L1 Loss (RGB, alpha=2.0)")
    print("-" * 80)
    disease_loss_rgb = DiseaseWeightedL1Loss(alpha=2.0, colorspace='RGB')
    loss_value = disease_loss_rgb(reconstructed, original)
    print(f"Loss value: {loss_value.item():.6f}")

    # Example 3: Disease-Weighted L1 Loss (LAB)
    print("\n3. Disease-Weighted L1 Loss (LAB, alpha=1.5)")
    print("-" * 80)
    disease_loss_lab = DiseaseWeightedL1Loss(alpha=1.5, colorspace='LAB')
    loss_value = disease_loss_lab(reconstructed, original)
    print(f"Loss value: {loss_value.item():.6f}")

    # Example 4: KL Divergence Loss
    print("\n4. KL Divergence Loss")
    print("-" * 80)
    kl_loss = KLDivergenceLoss()
    loss_value = kl_loss(mu, logvar)
    print(f"Loss value: {loss_value.item():.6f}")

    # Example 5: VAE Loss with L1 reconstruction
    print("\n5. VAE Loss (L1 reconstruction, beta=1.0)")
    print("-" * 80)
    vae_loss_l1 = VAELoss(
        reconstruction_loss=L1ReconstructionLoss(),
        beta=1.0
    )
    loss_dict = vae_loss_l1(
        reconstructed, original, mu, logvar,
        return_components=True
    )
    print(f"Total loss: {loss_dict['total'].item():.6f}")
    print(f"  Reconstruction: {loss_dict['reconstruction'].item():.6f}")
    print(f"  KL divergence: {loss_dict['kl'].item():.6f}")
    print(f"  Beta: {loss_dict['beta'].item():.2f}")

    # Example 6: VAE Loss with disease-weighted reconstruction
    print("\n6. VAE Loss (Disease-weighted reconstruction, beta=1.5)")
    print("-" * 80)
    vae_loss_disease = VAELoss(
        reconstruction_loss=DiseaseWeightedL1Loss(alpha=2.0, colorspace='RGB'),
        beta=1.5
    )
    loss_dict = vae_loss_disease(
        reconstructed, original, mu, logvar,
        return_components=True
    )
    print(f"Total loss: {loss_dict['total'].item():.6f}")
    print(f"  Reconstruction: {loss_dict['reconstruction'].item():.6f}")
    print(f"  KL divergence: {loss_dict['kl'].item():.6f}")
    print(f"  Beta: {loss_dict['beta'].item():.2f}")

    # Example 7: Create loss from config
    print("\n7. Create Loss from Configuration")
    print("-" * 80)

    # Example config
    example_config = {
        'loss': {
            'type': 'vae',
            'reconstruction_type': 'disease_weighted_l1',
            'disease_weight_alpha': 2.0,
            'beta': 1.0
        },
        'image': {
            'colorspace': 'RGB'
        }
    }

    loss_fn = create_loss_from_config(example_config)
    print(f"Created loss function: {loss_fn.__class__.__name__}")

    if isinstance(loss_fn, VAELoss):
        loss_dict = loss_fn(
            reconstructed, original, mu, logvar,
            return_components=True
        )
        print(f"Total loss: {loss_dict['total'].item():.6f}")

    print("\n" + "=" * 80)
    print("Examples completed successfully!")
    print("=" * 80)
