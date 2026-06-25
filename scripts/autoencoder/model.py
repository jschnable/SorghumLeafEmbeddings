"""
Convolutional Autoencoder (CAE) based on PDSE-Lite architecture.

This module implements a configurable convolutional autoencoder for leaf image
analysis, with optional attention mechanisms and variational autoencoder (VAE)
support.

Architecture Overview:
---------------------
Encoder:
    Conv2D #1: 16 filters → MaxPool2D (×½) → 128×128×16
    Conv2D #2: 8 filters → MaxPool2D (×½) → 64×64×8
    Conv2D #3: 8 filters → MaxPool2D (×½) → 32×32×8
    Bottleneck: 8 filters → 32×32×8

Decoder:
    UpSample2D #1 (×2) → Conv2D #5: 8 filters → 64×64×8
    UpSample2D #2 (×2) → Conv2D #6: 8 filters → 128×128×8
    UpSample2D #3 (×2) → Conv2D #7: 3 filters → 256×256×3

Reference: PDSE-Lite (Plant Disease Severity Estimation)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple, Optional, Dict, Any, List
import yaml
import math


class ChannelAttention(nn.Module):
    """
    Channel Attention Module (Squeeze-and-Excitation).

    Adaptively recalibrates channel-wise feature responses by explicitly
    modeling interdependencies between channels.

    Args:
        channels: Number of input channels
        reduction: Reduction ratio for bottleneck layer (default: 16)
    """

    def __init__(self, channels: int, reduction: int = 16):
        super().__init__()
        self.avg_pool = nn.AdaptiveAvgPool2d(1)
        self.max_pool = nn.AdaptiveMaxPool2d(1)

        # Shared MLP
        self.fc = nn.Sequential(
            nn.Linear(channels, channels // reduction, bias=False),
            nn.ReLU(inplace=True),
            nn.Linear(channels // reduction, channels, bias=False)
        )
        self.sigmoid = nn.Sigmoid()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Apply channel attention.

        Args:
            x: Input tensor (B, C, H, W)

        Returns:
            Attention-weighted tensor (B, C, H, W)
        """
        b, c, _, _ = x.size()

        # Average pooling branch
        avg_out = self.fc(self.avg_pool(x).view(b, c))

        # Max pooling branch
        max_out = self.fc(self.max_pool(x).view(b, c))

        # Combine and apply sigmoid
        attention = self.sigmoid(avg_out + max_out).view(b, c, 1, 1)

        return x * attention.expand_as(x)


class SpatialAttention(nn.Module):
    """
    Spatial Attention Module.

    Focuses on 'where' is informative by computing attention across spatial
    dimensions using both average and max pooling.

    Args:
        kernel_size: Size of convolutional kernel (default: 7)
    """

    def __init__(self, kernel_size: int = 7):
        super().__init__()
        padding = kernel_size // 2
        self.conv = nn.Conv2d(2, 1, kernel_size, padding=padding, bias=False)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Apply spatial attention.

        Args:
            x: Input tensor (B, C, H, W)

        Returns:
            Attention-weighted tensor (B, C, H, W)
        """
        # Compute channel-wise statistics
        avg_out = torch.mean(x, dim=1, keepdim=True)
        max_out, _ = torch.max(x, dim=1, keepdim=True)

        # Concatenate and apply convolution
        attention = torch.cat([avg_out, max_out], dim=1)
        attention = self.sigmoid(self.conv(attention))

        return x * attention


class AttentionBlock(nn.Module):
    """
    Combined Channel and Spatial Attention Block.

    Args:
        channels: Number of input channels
        attention_type: Type of attention ('channel', 'spatial', or 'both')
        reduction: Channel reduction ratio for channel attention
    """

    def __init__(
        self,
        channels: int,
        attention_type: str = 'both',
        reduction: int = 16
    ):
        super().__init__()
        self.attention_type = attention_type.lower()

        if self.attention_type in ['channel', 'both']:
            self.channel_attention = ChannelAttention(channels, reduction)

        if self.attention_type in ['spatial', 'both']:
            self.spatial_attention = SpatialAttention()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Apply attention mechanism(s).

        Args:
            x: Input tensor (B, C, H, W)

        Returns:
            Attention-weighted tensor (B, C, H, W)
        """
        if self.attention_type in ['channel', 'both']:
            x = self.channel_attention(x)

        if self.attention_type in ['spatial', 'both']:
            x = self.spatial_attention(x)

        return x


class EncoderBlock(nn.Module):
    """
    Encoder block consisting of Conv2D + ReLU + MaxPool + optional Attention.

    Args:
        in_channels: Number of input channels
        out_channels: Number of output channels
        kernel_size: Size of convolutional kernel
        use_attention: Whether to use attention mechanism
        attention_type: Type of attention if enabled
        use_pooling: Whether to apply max pooling (default: True)
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int = 3,
        use_attention: bool = False,
        attention_type: str = 'both',
        use_pooling: bool = True
    ):
        super().__init__()
        padding = kernel_size // 2
        self.use_pooling = use_pooling

        self.conv = nn.Conv2d(
            in_channels, out_channels,
            kernel_size=kernel_size,
            padding=padding,
            bias=True
        )
        self.relu = nn.ReLU(inplace=True)

        if use_pooling:
            self.pool = nn.MaxPool2d(kernel_size=2, stride=2)

        self.attention = None
        if use_attention:
            self.attention = AttentionBlock(out_channels, attention_type)

        # Initialize weights
        self._init_weights()

    def _init_weights(self):
        """Initialize convolutional weights using Kaiming initialization."""
        nn.init.kaiming_normal_(self.conv.weight, mode='fan_out', nonlinearity='relu')
        if self.conv.bias is not None:
            nn.init.constant_(self.conv.bias, 0)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass through encoder block.

        Args:
            x: Input tensor

        Returns:
            Encoded features
        """
        x = self.conv(x)
        x = self.relu(x)

        if self.use_pooling:
            x = self.pool(x)

        if self.attention is not None:
            x = self.attention(x)

        return x


class DecoderBlock(nn.Module):
    """
    Decoder block consisting of UpSample + Conv2D + ReLU + optional Attention.

    Args:
        in_channels: Number of input channels
        out_channels: Number of output channels
        kernel_size: Size of convolutional kernel
        use_attention: Whether to use attention mechanism
        attention_type: Type of attention if enabled
        use_upsampling: Whether to apply upsampling (default: True)
        final_layer: Whether this is the final layer (no ReLU)
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int = 3,
        use_attention: bool = False,
        attention_type: str = 'both',
        use_upsampling: bool = True,
        final_layer: bool = False
    ):
        super().__init__()
        padding = kernel_size // 2
        self.use_upsampling = use_upsampling
        self.final_layer = final_layer

        if use_upsampling:
            self.upsample = nn.Upsample(scale_factor=2, mode='nearest')

        self.conv = nn.Conv2d(
            in_channels, out_channels,
            kernel_size=kernel_size,
            padding=padding,
            bias=True
        )

        if not final_layer:
            self.relu = nn.ReLU(inplace=True)

        self.attention = None
        if use_attention and not final_layer:
            self.attention = AttentionBlock(out_channels, attention_type)

        # Initialize weights
        self._init_weights()

    def _init_weights(self):
        """Initialize convolutional weights using Kaiming initialization."""
        if not self.final_layer:
            nn.init.kaiming_normal_(self.conv.weight, mode='fan_out', nonlinearity='relu')
        else:
            nn.init.xavier_normal_(self.conv.weight)
        if self.conv.bias is not None:
            nn.init.constant_(self.conv.bias, 0)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass through decoder block.

        Args:
            x: Input tensor

        Returns:
            Decoded features
        """
        if self.use_upsampling:
            x = self.upsample(x)

        x = self.conv(x)

        if not self.final_layer:
            x = self.relu(x)

            if self.attention is not None:
                x = self.attention(x)

        return x


class ConvolutionalAutoencoder(nn.Module):
    """
    Convolutional Autoencoder based on PDSE-Lite architecture.

    This model implements a configurable CAE with optional attention mechanisms
    and variational autoencoder support for leaf image analysis.

    Architecture matches PDSE-Lite CAE block:
    - Encoder: Progressive downsampling with increasing depth
    - Bottleneck: Compressed latent representation
    - Decoder: Progressive upsampling with reconstruction

    Args:
        config: Configuration dictionary or path to YAML config file
    """

    def __init__(self, config: Dict[str, Any] | str):
        super().__init__()

        # Load config from file if path is provided
        if isinstance(config, str):
            with open(config, 'r') as f:
                config = yaml.safe_load(f)

        # Extract model configuration
        model_config = config.get('model', {})
        image_config = config.get('image', {})

        # Model parameters
        self.input_size = tuple(model_config.get('input_size', [256, 256]))
        self.input_channels = image_config.get('input_channels',
                                               image_config.get('channels', 3))
        self.bottleneck_dim = model_config.get('bottleneck_dim', 128)
        self.use_attention = model_config.get('use_attention', False)
        self.attention_type = model_config.get('attention_type', 'both')
        self.is_variational = model_config.get('is_vae', False)

        # Validate configuration
        self._validate_config()

        # Calculate filter sizes and spatial dimensions
        self.filters = self._calculate_filter_sizes()

        # Calculate spatial dimensions after encoder (before flattening)
        # After 3 pooling layers: input_size / 8
        self.encoded_h = self.input_size[0] // 8
        self.encoded_w = self.input_size[1] // 8
        self.encoded_channels = self.filters['enc3']
        self.encoded_dim = self.encoded_h * self.encoded_w * self.encoded_channels

        # Build encoder
        self.encoder = self._build_encoder()

        # Build bottleneck with flattening (different for VAE vs standard AE)
        if self.is_variational:
            # VAE: separate fully connected layers for mean and log-variance
            self.fc_mu = nn.Linear(self.encoded_dim, self.bottleneck_dim)
            self.fc_logvar = nn.Linear(self.encoded_dim, self.bottleneck_dim)
            # Decoder input from bottleneck
            self.fc_decode = nn.Linear(self.bottleneck_dim, self.encoded_dim)

            # Initialize VAE-specific layers to prevent posterior collapse
            nn.init.xavier_normal_(self.fc_mu.weight)
            nn.init.constant_(self.fc_mu.bias, 0)
            nn.init.xavier_normal_(self.fc_logvar.weight)
            nn.init.constant_(self.fc_logvar.bias, 0)
            nn.init.xavier_normal_(self.fc_decode.weight)
            nn.init.constant_(self.fc_decode.bias, 0)
        else:
            # Standard AE: single bottleneck
            self.fc_encode = nn.Linear(self.encoded_dim, self.bottleneck_dim)
            self.fc_decode = nn.Linear(self.bottleneck_dim, self.encoded_dim)

            # Initialize standard AE layers
            nn.init.xavier_normal_(self.fc_encode.weight)
            nn.init.constant_(self.fc_encode.bias, 0)
            nn.init.xavier_normal_(self.fc_decode.weight)
            nn.init.constant_(self.fc_decode.bias, 0)

        # Build decoder
        self.decoder = self._build_decoder()

    def _validate_config(self):
        """Validate configuration parameters."""
        # Check input size
        h, w = self.input_size
        if h != w:
            print(f"Warning: Non-square input size {self.input_size}. "
                  "Model may not work correctly.")

        # Check input size is divisible by 8 (3 pooling layers)
        if h % 8 != 0 or w % 8 != 0:
            print(f"Warning: Input size {self.input_size} not divisible by 8. "
                  "May cause issues with pooling/upsampling.")

        # Validate bottleneck dimension
        if not isinstance(self.bottleneck_dim, int) or self.bottleneck_dim <= 0:
            raise ValueError(f"bottleneck_dim must be positive integer, got {self.bottleneck_dim}")

    def _calculate_filter_sizes(self) -> Dict[str, int]:
        """
        Calculate number of filters for each layer.

        Based on PDSE-Lite architecture, scaled for different input sizes.

        Returns:
            Dictionary with filter counts for each layer
        """
        # PDSE-Lite original: 256x256 input
        # Encoder: 16 -> 8 -> 8
        # Decoder: 8 -> 8 -> 3 (output channels)

        # For now, use fixed architecture as in paper
        # Could scale based on input size if needed
        return {
            'enc1': 16,
            'enc2': 8,
            'enc3': 8,
            'dec1': 8,
            'dec2': 8,
            'dec3': self.input_channels
        }

    def _build_encoder(self) -> nn.ModuleList:
        """
        Build encoder network.

        Returns:
            ModuleList containing encoder blocks
        """
        encoder_blocks = nn.ModuleList([
            # Encoder Block 1: input_channels -> 16, downsample to 128x128
            EncoderBlock(
                self.input_channels, self.filters['enc1'],
                kernel_size=3,
                use_attention=self.use_attention,
                attention_type=self.attention_type,
                use_pooling=True
            ),

            # Encoder Block 2: 16 -> 8, downsample to 64x64
            EncoderBlock(
                self.filters['enc1'], self.filters['enc2'],
                kernel_size=3,
                use_attention=self.use_attention,
                attention_type=self.attention_type,
                use_pooling=True
            ),

            # Encoder Block 3: 8 -> 8, downsample to 32x32
            EncoderBlock(
                self.filters['enc2'], self.filters['enc3'],
                kernel_size=3,
                use_attention=self.use_attention,
                attention_type=self.attention_type,
                use_pooling=True
            )
        ])

        return encoder_blocks

    def _build_decoder(self) -> nn.ModuleList:
        """
        Build decoder network.

        Returns:
            ModuleList containing decoder blocks
        """
        decoder_blocks = nn.ModuleList([
            # Decoder Block 1: encoded_channels -> 8, upsample to 64x64
            DecoderBlock(
                self.encoded_channels, self.filters['dec1'],
                kernel_size=3,
                use_attention=self.use_attention,
                attention_type=self.attention_type,
                use_upsampling=True,
                final_layer=False
            ),

            # Decoder Block 2: 8 -> 8, upsample to 128x128
            DecoderBlock(
                self.filters['dec1'], self.filters['dec2'],
                kernel_size=3,
                use_attention=self.use_attention,
                attention_type=self.attention_type,
                use_upsampling=True,
                final_layer=False
            ),

            # Decoder Block 3: 8 -> input_channels, upsample to 256x256
            DecoderBlock(
                self.filters['dec2'], self.filters['dec3'],
                kernel_size=3,
                use_attention=False,  # No attention on final layer
                use_upsampling=True,
                final_layer=True
            )
        ])

        return decoder_blocks

    def encode(self, x: torch.Tensor) -> torch.Tensor | Tuple[torch.Tensor, torch.Tensor]:
        """
        Encode input to latent representation.

        Args:
            x: Input tensor (B, C, H, W)

        Returns:
            If variational: (mu, logvar) both (B, latent_dim)
            Otherwise: latent representation (B, latent_dim)
        """
        # Pass through encoder blocks
        for encoder_block in self.encoder:
            x = encoder_block(x)

        # Flatten spatial dimensions
        batch_size = x.size(0)
        x = x.view(batch_size, -1)  # (B, encoded_dim)

        # Bottleneck encoding
        if self.is_variational:
            mu = self.fc_mu(x)  # (B, bottleneck_dim)
            logvar = self.fc_logvar(x)  # (B, bottleneck_dim)
            return mu, logvar
        else:
            z = self.fc_encode(x)  # (B, bottleneck_dim)
            # NO activation function on latent code - it should be unbounded
            # to allow full representational capacity
            return z

    def reparameterize(
        self,
        mu: torch.Tensor,
        logvar: torch.Tensor
    ) -> torch.Tensor:
        """
        Reparameterization trick for VAE.

        z = μ + σ * ε, where ε ~ N(0, I)

        Args:
            mu: Mean of latent distribution (B, bottleneck_dim)
            logvar: Log variance of latent distribution (B, bottleneck_dim)

        Returns:
            Sampled latent vector (B, bottleneck_dim)
        """
        # Clamp logvar to prevent numerical instability
        # Prevents std from becoming too small (underflow) or too large (overflow)
        logvar = torch.clamp(logvar, min=-10, max=10)

        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """
        Decode latent representation to reconstruction.

        Args:
            z: Latent representation (B, bottleneck_dim)

        Returns:
            Reconstructed output (B, C, H, W)
        """
        # Decode from bottleneck to spatial representation
        batch_size = z.size(0)
        z = self.fc_decode(z)  # (B, encoded_dim)
        z = F.relu(z)

        # Reshape to spatial dimensions
        z = z.view(batch_size, self.encoded_channels, self.encoded_h, self.encoded_w)

        # Pass through decoder blocks
        for decoder_block in self.decoder:
            z = decoder_block(z)

        return z

    def forward(
        self,
        x: torch.Tensor
    ) -> torch.Tensor | Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Forward pass through the autoencoder.

        Args:
            x: Input tensor (B, C, H, W)

        Returns:
            If variational: (reconstruction, mu, logvar)
            Otherwise: reconstruction
        """
        if self.is_variational:
            # VAE forward pass
            mu, logvar = self.encode(x)
            z = self.reparameterize(mu, logvar)
            reconstruction = self.decode(z)
            return reconstruction, mu, logvar
        else:
            # Standard AE forward pass
            z = self.encode(x)
            reconstruction = self.decode(z)
            return reconstruction

    def get_latent_representation(
        self,
        x: torch.Tensor,
        deterministic: bool = True
    ) -> torch.Tensor:
        """
        Get latent representation for an input.

        Args:
            x: Input tensor (B, C, H, W)
            deterministic: If True (and VAE), return mean instead of sampling

        Returns:
            Latent representation (B, C, H, W)
        """
        with torch.no_grad():
            if self.is_variational:
                mu, logvar = self.encode(x)
                if deterministic:
                    return mu
                else:
                    return self.reparameterize(mu, logvar)
            else:
                return self.encode(x)

    def summary(self, input_shape: Optional[Tuple[int, int, int, int]] = None):
        """
        Print model architecture summary.

        Args:
            input_shape: Optional input shape (B, C, H, W).
                        If None, uses config default.
        """
        if input_shape is None:
            input_shape = (1, self.input_channels, *self.input_size)

        # Get device of model parameters
        device = next(self.parameters()).device

        print("=" * 80)
        print("Convolutional Autoencoder Architecture")
        print("=" * 80)
        print(f"Input shape: {input_shape}")
        print(f"Variational: {self.is_variational}")
        print(f"Attention: {self.use_attention} ({self.attention_type})")
        print(f"Bottleneck dim: {self.bottleneck_dim}")
        print(f"Encoded spatial dims: {self.encoded_h}×{self.encoded_w}×{self.encoded_channels}")
        print(f"Encoded flat dim: {self.encoded_dim}")
        print()

        # Create dummy input on same device as model
        x = torch.randn(input_shape, device=device)

        print("ENCODER:")
        print("-" * 80)
        for i, encoder_block in enumerate(self.encoder):
            x = encoder_block(x)
            print(f"  Block {i+1}: {tuple(x.shape)} "
                  f"({self.filters[f'enc{i+1}']} filters)")

        print(f"  Flatten:  ({input_shape[0]}, {self.encoded_dim})")

        print()
        if self.is_variational:
            print("BOTTLENECK (VAE):")
            print("-" * 80)
            mu, logvar = self.encode(torch.randn(input_shape, device=device))
            print(f"  fc_mu:     {tuple(mu.shape)}")
            print(f"  fc_logvar: {tuple(logvar.shape)}")
            z = self.reparameterize(mu, logvar)
            print(f"  z (sampled): {tuple(z.shape)}")
        else:
            print("BOTTLENECK:")
            print("-" * 80)
            z = self.encode(torch.randn(input_shape, device=device))
            print(f"  fc_encode: {tuple(z.shape)}")

        print()
        print("DECODER:")
        print("-" * 80)
        print(f"  fc_decode: ({input_shape[0]}, {self.encoded_dim})")
        print(f"  Reshape:   ({input_shape[0]}, {self.encoded_channels}, {self.encoded_h}, {self.encoded_w})")

        # Decode from latent
        if self.is_variational:
            x = self.decode(self.reparameterize(mu, logvar))
        else:
            x = self.decode(z)

        # Show decoder block shapes by re-running
        z_reshaped = torch.randn(input_shape[0], self.encoded_channels,
                                 self.encoded_h, self.encoded_w, device=device)
        for i, decoder_block in enumerate(self.decoder):
            z_reshaped = decoder_block(z_reshaped)
            print(f"  Block {i+1}: {tuple(z_reshaped.shape)} "
                  f"({self.filters[f'dec{i+1}']} filters)")

        print()
        print("=" * 80)
        total_params = sum(p.numel() for p in self.parameters())
        trainable_params = sum(p.numel() for p in self.parameters() if p.requires_grad)
        print(f"Total parameters: {total_params:,}")
        print(f"Trainable parameters: {trainable_params:,}")
        print("=" * 80)


def create_model_from_config(config: Dict[str, Any] | str) -> ConvolutionalAutoencoder:
    """
    Factory function to create model from configuration.

    Args:
        config: Either a configuration dictionary or path to YAML config file

    Returns:
        ConvolutionalAutoencoder instance
    """
    return ConvolutionalAutoencoder(config)


# ============================================================================
# Example Usage and Testing
# ============================================================================

if __name__ == "__main__":
    """
    Test the convolutional autoencoder with various configurations.
    """
    print("=" * 80)
    print("Convolutional Autoencoder Testing")
    print("=" * 80)

    # Test 1: Standard Autoencoder
    print("\n1. Standard Autoencoder (No Attention)")
    print("-" * 80)

    config_standard = {
        'model': {
            'input_size': [256, 256],
            'bottleneck_dim': 128,
            'use_attention': False,
            'is_vae': False
        },
        'image': {
            'input_channels': 3
        }
    }

    model = ConvolutionalAutoencoder(config_standard)
    model.summary()

    # Test forward pass
    dummy_input = torch.randn(2, 3, 256, 256)
    output = model(dummy_input)
    print(f"\nForward pass test:")
    print(f"  Input shape:  {tuple(dummy_input.shape)}")
    print(f"  Output shape: {tuple(output.shape)}")
    assert output.shape == dummy_input.shape, "Output shape mismatch!"
    print("  ✓ Shape test passed")

    # Test 2: Autoencoder with Attention
    print("\n" + "=" * 80)
    print("2. Autoencoder with Channel Attention")
    print("-" * 80)

    config_attention = {
        'model': {
            'input_size': [256, 256],
            'bottleneck_dim': 256,
            'use_attention': True,
            'attention_type': 'channel',
            'is_vae': False
        },
        'image': {
            'input_channels': 3
        }
    }

    model_attention = ConvolutionalAutoencoder(config_attention)
    model_attention.summary()

    output = model_attention(dummy_input)
    print(f"\nForward pass test:")
    print(f"  Input shape:  {tuple(dummy_input.shape)}")
    print(f"  Output shape: {tuple(output.shape)}")
    assert output.shape == dummy_input.shape, "Output shape mismatch!"
    print("  ✓ Shape test passed")

    # Test 3: Variational Autoencoder
    print("\n" + "=" * 80)
    print("3. Variational Autoencoder (VAE)")
    print("-" * 80)

    config_vae = {
        'model': {
            'input_size': [256, 256],
            'bottleneck_dim': 128,
            'use_attention': False,
            'is_vae': True
        },
        'image': {
            'input_channels': 3
        }
    }

    model_vae = ConvolutionalAutoencoder(config_vae)
    model_vae.summary()

    reconstruction, mu, logvar = model_vae(dummy_input)
    print(f"\nForward pass test:")
    print(f"  Input shape:         {tuple(dummy_input.shape)}")
    print(f"  Reconstruction shape: {tuple(reconstruction.shape)}")
    print(f"  mu shape:            {tuple(mu.shape)}")
    print(f"  logvar shape:        {tuple(logvar.shape)}")
    assert reconstruction.shape == dummy_input.shape, "Reconstruction shape mismatch!"
    assert mu.shape == (2, 128), "Mu shape mismatch!"
    assert logvar.shape == (2, 128), "Logvar shape mismatch!"
    print("  ✓ All shape tests passed")

    # Test latent representation
    latent = model_vae.get_latent_representation(dummy_input, deterministic=True)
    print(f"\nLatent representation test:")
    print(f"  Latent shape (deterministic): {tuple(latent.shape)}")
    latent_sampled = model_vae.get_latent_representation(dummy_input, deterministic=False)
    print(f"  Latent shape (sampled):       {tuple(latent_sampled.shape)}")
    print("  ✓ Latent extraction test passed")

    # Test 4: VAE with Both Attentions
    print("\n" + "=" * 80)
    print("4. VAE with Channel + Spatial Attention")
    print("-" * 80)

    config_vae_attention = {
        'model': {
            'input_size': [256, 256],
            'bottleneck_dim': 256,
            'use_attention': True,
            'attention_type': 'both',
            'is_vae': True
        },
        'image': {
            'input_channels': 3
        }
    }

    model_vae_attention = ConvolutionalAutoencoder(config_vae_attention)
    model_vae_attention.summary()

    reconstruction, mu, logvar = model_vae_attention(dummy_input)
    print(f"\nForward pass test:")
    print(f"  Input shape:         {tuple(dummy_input.shape)}")
    print(f"  Reconstruction shape: {tuple(reconstruction.shape)}")
    assert reconstruction.shape == dummy_input.shape, "Reconstruction shape mismatch!"
    print("  ✓ Shape test passed")

    # Parameter count comparison
    print("\n" + "=" * 80)
    print("Parameter Comparison")
    print("=" * 80)
    print(f"Standard AE:              {sum(p.numel() for p in model.parameters()):,} params")
    print(f"AE + Channel Attention:   {sum(p.numel() for p in model_attention.parameters()):,} params")
    print(f"VAE:                      {sum(p.numel() for p in model_vae.parameters()):,} params")
    print(f"VAE + Both Attentions:    {sum(p.numel() for p in model_vae_attention.parameters()):,} params")

    print("\n" + "=" * 80)
    print("All tests passed successfully!")
    print("=" * 80)
