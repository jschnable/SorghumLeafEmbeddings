"""Render the current three-panel Figure 2 SVG and export its PNG."""
from pathlib import Path
import subprocess

DIR = Path(__file__).resolve().parent
subprocess.run(['Rscript', str(DIR / 'figure2.R')], check=True)
subprocess.run(['inkscape', str(DIR / 'figure2.svg'), '--export-width=1950',
                '--export-background=white',
                '--export-filename=' + str(DIR / 'figure2.png')], check=True)
