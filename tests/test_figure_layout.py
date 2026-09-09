"""Distribution checks for standalone figure directories."""
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_bundled_figure_inputs_stay_below_one_mb():
    # Include new files before staging, but exclude ignored generated artifacts.
    paths = subprocess.check_output(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard', 'figures'],
        cwd=ROOT, text=True,
    ).splitlines()
    sizes = {}
    non_input_extensions = {'.png', '.jpg', '.jpeg', '.svg', '.pdf', '.tif', '.tiff',
                            '.webp', '.py', '.r', '.md'}
    for name in set(paths):
        path = ROOT / name
        if not path.is_file() or path.suffix.lower() in non_input_extensions:
            continue
        figure = Path(*Path(name).parts[:3])
        sizes[figure] = sizes.get(figure, 0) + path.stat().st_size
    oversized = {str(figure): size for figure, size in sizes.items() if size >= 1_000_000}
    assert not oversized, f'Figure inputs exceed 1 MB (images excluded): {oversized}'


def test_gallery_renders_from_its_directory_alone(tmp_path):
    source = ROOT / 'figures/supplemental/FigS3_exg_leaf_gallery'
    destination = tmp_path / 'gallery'
    destination.mkdir()
    for name in ['make_grid.py', 'layout.csv']:
        shutil.copyfile(source / name, destination / name)
    shutil.copytree(source / 'panels', destination / 'panels')
    subprocess.run([sys.executable, str(destination / 'make_grid.py')],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    assert (destination / 'grid.png').is_file()
