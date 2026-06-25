"""
Utility functions for creating data splits and DataLoaders.

This module provides functions to:
- Load configuration from YAML
- Create genotype-based train/val/test splits
- Create PyTorch DataLoaders for each split
"""

import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

import yaml
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader, Subset

from dataset import LeafImageDataset


def load_config(config_path: str = "src/autoencoder/config.yaml") -> Dict[str, Any]:
    """
    Load configuration from YAML file.

    Args:
        config_path: Path to configuration file

    Returns:
        Dictionary containing configuration parameters
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_image_directory(
    base_dir: str,
    environment: str,
    subdir: str,
    colorspace: str,
    use_masked: bool,
    normalize_colorspace: bool = True
) -> List[str]:
    """
    Construct path(s) to processed image directory/directories.

    Expected base layout: {base_dir}/{environment}/cropped[_{colorspace}_normalized][_masked]/
    e.g., data/generatable/images_processed/ne2025/cropped/

    Supports glob patterns in subdir to match multiple subdirectories.

    Args:
        base_dir: Base data directory (e.g., 'data/generatable/images_processed')
        environment: Environment name (e.g., 'aamu2025', 'fvsu2025', 'ne2025')
        subdir: Subdirectory within environment (leave empty; glob patterns also supported)
        colorspace: Colorspace name (e.g., 'RGB')
        use_masked: Whether to use masked images
        normalize_colorspace: Whether images are color-normalized

    Returns:
        List of paths to image directories (multiple if subdir is a glob pattern)
    """
    from glob import glob

    # Build directory name based on normalize_colorspace
    if normalize_colorspace:
        # Use normalized directories: cropped_{colorspace}_normalized[_masked]
        dir_name = f"cropped_{colorspace}_normalized"
        if use_masked:
            dir_name += "_masked"
    else:
        # Use raw cropped directory (masking will be applied in dataset if use_masked=True)
        dir_name = "cropped"

    # Build full path(s)
    if subdir:
        # Check if subdir contains glob pattern characters
        if any(char in subdir for char in ['*', '?', '[', ']']):
            # Remove trailing slash if present
            subdir_pattern = subdir.rstrip('/')
            # Construct glob pattern
            pattern = os.path.join(base_dir, environment, subdir_pattern, dir_name)
            paths = sorted(glob(pattern))
            if not paths:
                # Try without the dir_name to see if subdirs exist
                base_pattern = os.path.join(base_dir, environment, subdir_pattern)
                matching_subdirs = sorted(glob(base_pattern))
                if matching_subdirs:
                    # Return the full paths even if cropped dirs don't exist yet
                    paths = [os.path.join(sd, dir_name) for sd in matching_subdirs]
        else:
            # Single subdirectory
            paths = [os.path.join(base_dir, environment, subdir, dir_name)]
    else:
        # No subdirectory
        paths = [os.path.join(base_dir, environment, dir_name)]

    return paths


def load_image_paths(image_dirs: List[str], extension: str = '.png', normalize_colorspace: bool = False) -> List[str]:
    """
    Load all image paths from one or more directories.

    Args:
        image_dirs: Directory or list of directories containing images
        extension: Image file extension (default: '.png', or '.npz' for normalized images)
        normalize_colorspace: If True, looks for .npz files; if False, uses extension parameter

    Returns:
        List of absolute paths to image files from all directories
    """
    # Handle both single directory (str) and multiple directories (list)
    if isinstance(image_dirs, str):
        image_dirs = [image_dirs]

    # Override extension if normalized colorspace (NPZ files)
    if normalize_colorspace:
        extension = '.npz'

    all_image_paths = []

    for image_dir in image_dirs:
        image_dir_path = Path(image_dir)

        if not image_dir_path.exists():
            print(f"Warning: Image directory does not exist: {image_dir_path}")
            continue

        image_paths = sorted(image_dir_path.glob(f"*{extension}"))
        image_paths = [str(p.absolute()) for p in image_paths]

        if len(image_paths) == 0:
            print(f"Warning: No images found in {image_dir_path} with extension {extension}")
        else:
            all_image_paths.extend(image_paths)

    if len(all_image_paths) == 0:
        raise ValueError(f"No images found in any of the directories: {image_dirs}")

    return sorted(all_image_paths)


def create_genotype_splits(
    image_paths: List[str],
    field_index_df: pd.DataFrame,
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    random_seed: int = 42,
    verbose: bool = True
) -> Tuple[List[int], List[int], List[int]]:
    """
    Create train/val/test splits based on genotypes.

    Each genotype appears in only ONE split to prevent data leakage.

    Args:
        image_paths: List of image file paths
        field_index_df: DataFrame with plotNumber to genotype mapping
        train_ratio: Fraction of genotypes for training
        val_ratio: Fraction of genotypes for validation
        test_ratio: Fraction of genotypes for testing
        random_seed: Random seed for reproducibility
        verbose: Whether to print split statistics

    Returns:
        Tuple of (train_indices, val_indices, test_indices)
    """
    # Validate ratios
    if not np.isclose(train_ratio + val_ratio + test_ratio, 1.0):
        raise ValueError(
            f"Split ratios must sum to 1.0, got {train_ratio + val_ratio + test_ratio}"
        )

    # Create plotNumber to genotype mapping
    plot_to_genotype = dict(zip(
        field_index_df['plotNumber'].astype(str),
        field_index_df['genotype']
    ))

    # Extract genotypes from filenames using field index mapping
    # Plot number is the first part before the first underscore
    plot_numbers = [os.path.basename(path).split('_')[0] for path in image_paths]
    genotypes = [plot_to_genotype.get(plot, plot) for plot in plot_numbers]

    # Get unique genotypes
    unique_genotypes = sorted(set(genotypes))
    n_genotypes = len(unique_genotypes)

    if verbose:
        print(f"Total images: {len(image_paths)}")
        print(f"Unique genotypes: {n_genotypes}")

    # Shuffle genotypes
    rng = np.random.RandomState(random_seed)
    shuffled_genotypes = unique_genotypes.copy()
    rng.shuffle(shuffled_genotypes)

    # Calculate split sizes
    n_train = int(n_genotypes * train_ratio)
    n_val = int(n_genotypes * val_ratio)
    # Ensure all genotypes are assigned
    n_test = n_genotypes - n_train - n_val

    # Split genotypes
    train_genotypes = set(shuffled_genotypes[:n_train])
    val_genotypes = set(shuffled_genotypes[n_train:n_train + n_val])
    test_genotypes = set(shuffled_genotypes[n_train + n_val:])

    if verbose:
        print(f"\nGenotype split:")
        print(f"  Train: {len(train_genotypes)} genotypes ({len(train_genotypes)/n_genotypes*100:.1f}%)")
        print(f"  Val:   {len(val_genotypes)} genotypes ({len(val_genotypes)/n_genotypes*100:.1f}%)")
        print(f"  Test:  {len(test_genotypes)} genotypes ({len(test_genotypes)/n_genotypes*100:.1f}%)")

    # Assign indices to splits based on genotype
    train_indices = []
    val_indices = []
    test_indices = []

    for idx, genotype in enumerate(genotypes):
        if genotype in train_genotypes:
            train_indices.append(idx)
        elif genotype in val_genotypes:
            val_indices.append(idx)
        elif genotype in test_genotypes:
            test_indices.append(idx)

    if verbose:
        print(f"\nImage split:")
        print(f"  Train: {len(train_indices)} images ({len(train_indices)/len(image_paths)*100:.1f}%)")
        print(f"  Val:   {len(val_indices)} images ({len(val_indices)/len(image_paths)*100:.1f}%)")
        print(f"  Test:  {len(test_indices)} images ({len(test_indices)/len(image_paths)*100:.1f}%)")

    return train_indices, val_indices, test_indices


def create_datasets_from_config(
    config: Optional[Dict[str, Any]] = None,
    config_path: str = "src/autoencoder/config.yaml"
) -> Tuple[LeafImageDataset, LeafImageDataset, LeafImageDataset]:
    """
    Create train, val, and test datasets from configuration.

    Args:
        config: Configuration dictionary (if None, will load from config_path)
        config_path: Path to configuration file

    Returns:
        Tuple of (train_dataset, val_dataset, test_dataset)
    """
    # Load config if not provided
    if config is None:
        config = load_config(config_path)

    # Extract parameters
    base_dir = config['data']['base_dir']
    environment = config['data']['environment']
    subdir = config['data'].get('subdir', '')
    field_index_csv = config['data']['field_index_csv']
    metadata_csv = config['data'].get('metadata_csv', '')

    colorspace = config['image']['colorspace']
    use_masked = config['image']['use_masked']
    normalize_colorspace = config['image'].get('normalize_colorspace', True)
    target_size = tuple(config['image']['target_size'])

    train_ratio = config['split']['train_ratio']
    val_ratio = config['split']['val_ratio']
    test_ratio = config['split']['test_ratio']
    random_seed = config['split']['random_seed']

    cache_images = config['settings']['cache_images']
    verbose = config['settings']['verbose']

    # Get image directory/directories and load paths
    image_dirs = get_image_directory(base_dir, environment, subdir, colorspace, use_masked, normalize_colorspace)

    if verbose:
        if len(image_dirs) == 1:
            print(f"Loading images from: {image_dirs[0]}")
        else:
            print(f"Loading images from {len(image_dirs)} directories:")
            for img_dir in image_dirs:
                print(f"  - {img_dir}")

    # Pass normalize_colorspace to determine file extension (.npz vs .png)
    image_paths = load_image_paths(image_dirs, normalize_colorspace=normalize_colorspace)

    # Load field index (required)
    if verbose:
        print(f"Loading field index from: {field_index_csv}")

    field_index_df = pd.read_csv(field_index_csv)

    # Validate field index has required columns
    if 'plotNumber' not in field_index_df.columns:
        raise ValueError(f"field_index_csv must have 'plotNumber' column. Found: {field_index_df.columns.tolist()}")
    if 'genotype' not in field_index_df.columns:
        raise ValueError(f"field_index_csv must have 'genotype' column. Found: {field_index_df.columns.tolist()}")

    if verbose:
        print(f"  Found {len(field_index_df)} plot-genotype mappings")
        print(f"  Unique genotypes: {field_index_df['genotype'].nunique()}")

    # Load additional metadata (optional)
    metadata_df = None
    if metadata_csv and os.path.exists(metadata_csv):
        if verbose:
            print(f"Loading additional metadata from: {metadata_csv}")
        metadata_df = pd.read_csv(metadata_csv)

    # Create splits
    train_indices, val_indices, test_indices = create_genotype_splits(
        image_paths,
        field_index_df=field_index_df,
        train_ratio=train_ratio,
        val_ratio=val_ratio,
        test_ratio=test_ratio,
        random_seed=random_seed,
        verbose=verbose
    )

    # Create datasets
    train_paths = [image_paths[i] for i in train_indices]
    val_paths = [image_paths[i] for i in val_indices]
    test_paths = [image_paths[i] for i in test_indices]

    # Determine mask directory if needed
    mask_dir = None
    if not normalize_colorspace and use_masked:
        # Build path to masks_cropped directory
        if subdir:
            if any(char in subdir for char in ['*', '?', '[', ']']):
                # For glob patterns, we'll need to construct mask paths per image
                # Set a flag that will be checked in the dataset
                mask_dir = "masks_cropped"  # Relative name, will be resolved per image
            else:
                mask_dir = os.path.join(base_dir, environment, subdir, "masks_cropped")
        else:
            mask_dir = os.path.join(base_dir, environment, "masks_cropped")

    train_dataset = LeafImageDataset(
        image_paths=train_paths,
        field_index_df=field_index_df,
        metadata_df=metadata_df,
        target_size=target_size,
        colorspace=colorspace,
        use_masked=use_masked,
        normalize_colorspace=normalize_colorspace,
        mask_dir=mask_dir,
        cache_images=cache_images
    )

    val_dataset = LeafImageDataset(
        image_paths=val_paths,
        field_index_df=field_index_df,
        metadata_df=metadata_df,
        target_size=target_size,
        colorspace=colorspace,
        use_masked=use_masked,
        normalize_colorspace=normalize_colorspace,
        mask_dir=mask_dir,
        cache_images=cache_images
    )

    test_dataset = LeafImageDataset(
        image_paths=test_paths,
        field_index_df=field_index_df,
        metadata_df=metadata_df,
        target_size=target_size,
        colorspace=colorspace,
        use_masked=use_masked,
        normalize_colorspace=normalize_colorspace,
        mask_dir=mask_dir,
        cache_images=cache_images
    )

    return train_dataset, val_dataset, test_dataset


def create_dataloaders_from_config(
    config: Optional[Dict[str, Any]] = None,
    config_path: str = "src/autoencoder/config.yaml"
) -> Tuple[DataLoader, DataLoader, DataLoader]:
    """
    Create train, val, and test DataLoaders from configuration.

    Args:
        config: Configuration dictionary (if None, will load from config_path)
        config_path: Path to configuration file

    Returns:
        Tuple of (train_loader, val_loader, test_loader)
    """
    # Load config if not provided
    if config is None:
        config = load_config(config_path)

    # Create datasets
    train_dataset, val_dataset, test_dataset = create_datasets_from_config(
        config, config_path
    )

    # Extract DataLoader parameters
    batch_size = config['dataloader']['batch_size']
    num_workers = config['dataloader']['num_workers']
    shuffle_train = config['dataloader']['shuffle_train']
    shuffle_val = config['dataloader']['shuffle_val']
    shuffle_test = config['dataloader']['shuffle_test']
    drop_last = config['dataloader']['drop_last']
    pin_memory = config['dataloader']['pin_memory']

    # Create DataLoaders
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=shuffle_train,
        num_workers=num_workers,
        drop_last=drop_last,
        pin_memory=pin_memory
    )

    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=shuffle_val,
        num_workers=num_workers,
        drop_last=False,  # Never drop last for validation
        pin_memory=pin_memory
    )

    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=shuffle_test,
        num_workers=num_workers,
        drop_last=False,  # Never drop last for test
        pin_memory=pin_memory
    )

    return train_loader, val_loader, test_loader


def verify_no_genotype_leakage(
    train_dataset: LeafImageDataset,
    val_dataset: LeafImageDataset,
    test_dataset: LeafImageDataset
) -> bool:
    """
    Verify that no genotype appears in multiple splits.

    Args:
        train_dataset: Training dataset
        val_dataset: Validation dataset
        test_dataset: Test dataset

    Returns:
        True if no leakage detected, False otherwise
    """
    train_genotypes = set(train_dataset.get_genotypes())
    val_genotypes = set(val_dataset.get_genotypes())
    test_genotypes = set(test_dataset.get_genotypes())

    # Check for overlaps
    train_val_overlap = train_genotypes & val_genotypes
    train_test_overlap = train_genotypes & test_genotypes
    val_test_overlap = val_genotypes & test_genotypes

    has_leakage = False

    if train_val_overlap:
        print(f"WARNING: {len(train_val_overlap)} genotypes in both train and val!")
        print(f"  Overlapping genotypes: {list(train_val_overlap)[:5]}")
        has_leakage = True

    if train_test_overlap:
        print(f"WARNING: {len(train_test_overlap)} genotypes in both train and test!")
        print(f"  Overlapping genotypes: {list(train_test_overlap)[:5]}")
        has_leakage = True

    if val_test_overlap:
        print(f"WARNING: {len(val_test_overlap)} genotypes in both val and test!")
        print(f"  Overlapping genotypes: {list(val_test_overlap)[:5]}")
        has_leakage = True

    if not has_leakage:
        print("✓ No genotype leakage detected!")
        print(f"  Train genotypes: {len(train_genotypes)}")
        print(f"  Val genotypes:   {len(val_genotypes)}")
        print(f"  Test genotypes:  {len(test_genotypes)}")

    return not has_leakage
