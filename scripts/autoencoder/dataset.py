"""
PyTorch Dataset class for preprocessed leaf images.

This module provides a Dataset class for loading preprocessed leaf images
with genotype-based splitting to prevent data leakage.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import cv2
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
from PIL import Image


class LeafImageDataset(Dataset):
    """
    Dataset for preprocessed leaf images with metadata.

    Images are loaded from preprocessed directories and resized to target dimensions.
    Each sample returns both the image tensor and associated metadata from CSV.

    Attributes:
        image_paths: List of paths to image files
        metadata_df: DataFrame containing image metadata
        target_size: Tuple of (height, width) for resizing
        colorspace: Colorspace of images ('RGB', 'LAB', etc.)
        use_masked: Whether using masked images
        normalize_colorspace: Whether images are color-normalized
        mask_dir: Directory containing masks (used when normalize_colorspace=False and use_masked=True)
        cache_images: Whether to cache loaded images in memory
        image_cache: Dictionary storing cached images if enabled
    """

    def __init__(
        self,
        image_paths: List[str],
        field_index_df: pd.DataFrame,
        metadata_df: Optional[pd.DataFrame] = None,
        target_size: Tuple[int, int] = (256, 256),
        colorspace: str = 'RGB',
        use_masked: bool = False,
        normalize_colorspace: bool = True,
        mask_dir: Optional[str] = None,
        cache_images: bool = False,
    ):
        """
        Initialize the dataset.

        Args:
            image_paths: List of paths to image files
            field_index_df: DataFrame with plotNumber to genotype mapping (required columns: plotNumber, genotype)
            metadata_df: Optional DataFrame with additional metadata
            target_size: Tuple of (height, width) for resizing images
            colorspace: Colorspace of the images
            use_masked: Whether images are masked
            normalize_colorspace: Whether images are color-normalized
            mask_dir: Directory containing masks (used when normalize_colorspace=False and use_masked=True)
            cache_images: Whether to cache loaded images in memory
        """
        self.image_paths = image_paths
        self.field_index_df = field_index_df
        self.metadata_df = metadata_df
        self.target_size = target_size
        self.colorspace = colorspace
        self.use_masked = use_masked
        self.normalize_colorspace = normalize_colorspace
        self.mask_dir = mask_dir
        self.cache_images = cache_images
        self.image_cache: Dict[int, torch.Tensor] = {}

        # Validate field_index has required columns
        if 'plotNumber' not in self.field_index_df.columns:
            raise ValueError("field_index_df must have 'plotNumber' column")
        if 'genotype' not in self.field_index_df.columns:
            raise ValueError("field_index_df must have 'genotype' column")

        # Create plotNumber to genotype mapping
        self.plot_to_genotype = dict(zip(
            self.field_index_df['plotNumber'].astype(str),
            self.field_index_df['genotype']
        ))

        # Create a mapping from image basename to metadata row (if metadata provided)
        if self.metadata_df is not None:
            self.metadata_df['_basename'] = self.metadata_df['image_name'].apply(
                lambda x: os.path.basename(x) if pd.notna(x) else ''
            )

            # Verify all images have metadata
            missing_metadata = []
            for img_path in self.image_paths:
                basename = os.path.basename(img_path)
                # Remove file extension and leaf number suffix for matching
                base_name_no_ext = self._get_base_name(basename)
                matching_rows = self.metadata_df[
                    self.metadata_df['image_name'].str.startswith(base_name_no_ext)
                ]
                if len(matching_rows) == 0:
                    missing_metadata.append(basename)

            if missing_metadata and len(missing_metadata) < 10:
                print(f"Warning: {len(missing_metadata)} images missing metadata: {missing_metadata[:5]}")
            elif missing_metadata:
                print(f"Warning: {len(missing_metadata)} images missing metadata")

        # Verify all images have genotype mapping
        missing_genotypes = []
        for img_path in self.image_paths:
            plot_number = self._get_plot_number_from_name(img_path)
            if plot_number not in self.plot_to_genotype:
                missing_genotypes.append(os.path.basename(img_path))

        if missing_genotypes and len(missing_genotypes) < 10:
            print(f"Warning: {len(missing_genotypes)} images missing genotype mapping: {missing_genotypes[:5]}")
        elif missing_genotypes:
            print(f"Warning: {len(missing_genotypes)} images missing genotype mapping")

    def _get_base_name(self, filename: str) -> str:
        """
        Extract base name from filename for metadata matching.

        Example: '1201_LeafPhotoA_2025-09-08 10_44_12.793-05_00_0.png'
                 -> '1201_LeafPhotoA_2025-09-08 10_44_12.793-05_00'

        Args:
            filename: Image filename

        Returns:
            Base name without extension and leaf number
        """
        # Remove .png extension
        name_no_ext = os.path.splitext(filename)[0]

        # Remove trailing _N suffix (leaf number)
        parts = name_no_ext.rsplit('_', 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0]
        return name_no_ext

    def _get_plot_number_from_name(self, filename: str) -> str:
        """
        Extract plot number from filename.

        Plot number is the first part before the first underscore.
        Example: '1201_LeafPhotoA_...' -> '1201'

        Args:
            filename: Image filename or path

        Returns:
            Plot number identifier
        """
        basename = os.path.basename(filename)
        return basename.split('_')[0]

    def _get_genotype_from_name(self, filename: str) -> str:
        """
        Extract genotype ID from filename using field index mapping.

        First extracts plot number from filename, then maps to genotype.
        Example: '1201_LeafPhotoA_...' -> '1201' -> looks up genotype in field_index

        Args:
            filename: Image filename or path

        Returns:
            Genotype identifier (or plot number if not found in mapping)
        """
        plot_number = self._get_plot_number_from_name(filename)
        # Return genotype from mapping, or plot_number if not found
        return self.plot_to_genotype.get(plot_number, plot_number)

    def __len__(self) -> int:
        """Return the number of images in the dataset."""
        return len(self.image_paths)

    def __getitem__(self, idx: int) -> Dict[str, Any]:
        """
        Load and return an image with its metadata.

        Args:
            idx: Index of the sample to load

        Returns:
            Dictionary containing:
                - 'image': torch.Tensor of shape (3, H, W) with values in [0, 1]
                - 'metadata': Dictionary with all metadata fields
                - 'image_path': Path to the image file
                - 'genotype': Genotype identifier
        """
        # Check cache first
        if self.cache_images and idx in self.image_cache:
            image_tensor = self.image_cache[idx]
        else:
            # Load and process image
            img_path = self.image_paths[idx]

            # Check file extension to determine loading method
            if img_path.endswith('.npz'):
                # Load normalized image from NPZ file (already in [0, 1] range and RGB)
                image = np.load(img_path)['image']
                if image is None or image.size == 0:
                    raise ValueError(f"Failed to load image: {img_path}")
            elif img_path.endswith('.npy'):
                # Load normalized image from NPY file (already in [0, 1] range and RGB)
                image = np.load(img_path)
                if image is None or image.size == 0:
                    raise ValueError(f"Failed to load image: {img_path}")
            else:
                # Load from PNG/JPG using OpenCV
                image = cv2.imread(img_path, cv2.IMREAD_COLOR)

                if image is None:
                    raise ValueError(f"Failed to load image: {img_path}")

                # Convert BGR to RGB
                image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

            # Apply mask if needed (normalize_colorspace=False and use_masked=True)
            if not self.normalize_colorspace and self.use_masked and self.mask_dir is not None:
                # Construct path to mask file
                basename = os.path.basename(img_path)

                # Handle case where mask_dir is just "masks_cropped" (relative)
                if self.mask_dir == "masks_cropped":
                    # Determine mask path from image path
                    img_dir = os.path.dirname(img_path)
                    parent_dir = os.path.dirname(img_dir)
                    mask_path = os.path.join(parent_dir, "masks_cropped", basename)
                else:
                    # mask_dir is absolute path
                    mask_path = os.path.join(self.mask_dir, basename)

                # Load mask
                if os.path.exists(mask_path):
                    mask = cv2.imread(mask_path, cv2.IMREAD_COLOR)
                    if mask is not None:
                        # Convert BGR to RGB
                        mask = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)

                        # Resize mask to match image if needed
                        if mask.shape[:2] != image.shape[:2]:
                            mask = cv2.resize(
                                mask,
                                (image.shape[1], image.shape[0]),
                                interpolation=cv2.INTER_NEAREST
                            )

                        # Convert mask to binary (0 or 1)
                        # Mask is [255, 255, 255] for leaf pixels, [0, 0, 0] for background
                        mask_binary = (mask / 255.0).astype(np.float32)

                        # Apply mask by element-wise multiplication
                        image = image.astype(np.float32) * mask_binary
                    else:
                        print(f"Warning: Failed to load mask: {mask_path}")
                else:
                    print(f"Warning: Mask not found: {mask_path}")

            # Resize to target size
            if image.shape[:2] != self.target_size:
                image = cv2.resize(
                    image,
                    (self.target_size[1], self.target_size[0]),  # cv2 uses (width, height)
                    interpolation=cv2.INTER_LINEAR
                )

            # Convert to float32 and normalize to [0, 1] if needed
            if img_path.endswith('.npz') or img_path.endswith('.npy'):
                # NPZ/NPY files are already normalized float32 [0, 1]
                image = image.astype(np.float32)
            elif image.dtype == np.uint8:
                # PNG/JPG files need normalization
                image = image.astype(np.float32) / 255.0
            else:
                image = image.astype(np.float32)

            # Convert to torch tensor and rearrange to (C, H, W)
            image_tensor = torch.from_numpy(image).permute(2, 0, 1)

            # Cache if enabled
            if self.cache_images:
                self.image_cache[idx] = image_tensor

        # Get metadata
        img_path = self.image_paths[idx]
        basename = os.path.basename(img_path)
        base_name_no_ext = self._get_base_name(basename)
        plot_number = self._get_plot_number_from_name(basename)

        # Build metadata dictionary
        metadata_row = {'image_name': base_name_no_ext, 'plotNumber': plot_number}

        # Add field index data (genotype and any other columns)
        if plot_number in self.plot_to_genotype:
            # Find the row in field_index for this plot
            field_row = self.field_index_df[
                self.field_index_df['plotNumber'].astype(str) == plot_number
            ]
            if len(field_row) > 0:
                field_data = field_row.iloc[0].to_dict()
                metadata_row.update(field_data)

        # Add additional metadata if provided
        if self.metadata_df is not None:
            matching_rows = self.metadata_df[
                self.metadata_df['image_name'].str.startswith(base_name_no_ext)
            ]
            if len(matching_rows) > 0:
                extra_metadata = matching_rows.iloc[0].to_dict()
                # Remove internal columns
                extra_metadata.pop('_basename', None)
                # Update with additional metadata (don't overwrite field_index data)
                for key, value in extra_metadata.items():
                    if key not in metadata_row:
                        metadata_row[key] = value

        # Extract genotype
        genotype = self._get_genotype_from_name(basename)

        return {
            'image': image_tensor,
            'metadata': metadata_row,
            'image_path': img_path,
            'genotype': genotype,
        }

    def get_genotypes(self) -> List[str]:
        """
        Get list of all unique genotypes in the dataset.

        Returns:
            List of unique genotype identifiers
        """
        genotypes = [self._get_genotype_from_name(path) for path in self.image_paths]
        return list(set(genotypes))

    def get_samples_by_genotype(self, genotype: str) -> List[int]:
        """
        Get indices of all samples belonging to a specific genotype.

        Args:
            genotype: Genotype identifier

        Returns:
            List of sample indices
        """
        indices = []
        for idx, path in enumerate(self.image_paths):
            if self._get_genotype_from_name(path) == genotype:
                indices.append(idx)
        return indices
