"""Refresh Figure 2's statistical panels, retaining its hand-drawn workflow panel."""
from pathlib import Path
import base64
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT/'figures/main/Fig2_embeddings'
SVG = 'http://www.w3.org/2000/svg'
XLINK = 'http://www.w3.org/1999/xlink'
ET.register_namespace('', SVG)
ET.register_namespace('xlink', XLINK)
tree = ET.parse(DIR/'figure2.svg')
layer = next(e for e in tree.getroot() if e.attrib.get('id') == 'layer1')
# The retained groups are exclusively the original workflow boxes/arrows and a/b labels.
keep = {'g12', 'g12-4', 'g14', 'path21-6', 'path21-6-6', 'text877-1', 'text877-5'}
for element in list(layer):
    if element.attrib.get('id') not in keep:
        layer.remove(element)
for name, file, x, y, width, height in [
    ('current_rf_accuracy', 'rf_accuracy.png', 43.3, 0, 83.0, 46.5),
    ('current_correlations', 'sam3_cor.png', 0.15, 46.0, 125.8, 58.45),
]:
    ET.SubElement(layer, '{'+SVG+'}image', {
        'id': name, 'x': str(x), 'y': str(y), 'width': str(width), 'height': str(height),
        'preserveAspectRatio': 'xMidYMid meet',
        '{'+XLINK+'}href': 'data:image/png;base64,'+base64.b64encode((DIR/file).read_bytes()).decode(),
    })
# Keep panel b's original letter in front of its new plot.
label = next(e for e in layer if e.attrib.get('id') == 'text877-1')
layer.remove(label)
layer.append(label)
tree.write(DIR/'figure2.svg', encoding='utf-8', xml_declaration=True)
subprocess.run(['inkscape', str(DIR/'figure2.svg'), '--export-width=1950',
                '--export-background=white', '--export-filename='+str(DIR/'figure2.png')], check=True)
