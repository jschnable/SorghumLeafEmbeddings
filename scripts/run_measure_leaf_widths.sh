#!/usr/bin/env bash
# Run measure_leaf_width.py over every raw image directory
# (data/externalsourcerequired/{environment}/{subdir}) and write one CSV per
# directory to data/generatable/leaf_widths/{environment}_{subdir}_leaf_widths.csv
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_root="$repo_root/data/externalsourcerequired"
out_dir="$repo_root/data/generatable/leaf_widths"

mkdir -p "$out_dir"

for env_dir in "$raw_root"/*/; do
    environment="$(basename "$env_dir")"

    for subdir in "$env_dir"*/; do
        [ -d "$subdir" ] || continue
        subdir_name="$(basename "$subdir")"
        out_csv="$out_dir/${environment}_${subdir_name}_leaf_widths.csv"

        echo "Measuring leaf widths: $environment/$subdir_name -> $out_csv"
        python "$repo_root/scripts/measure_leaf_width.py" "$subdir" --output "$out_csv"
    done
done
