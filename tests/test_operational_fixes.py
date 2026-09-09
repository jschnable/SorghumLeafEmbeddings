from pathlib import Path
from types import SimpleNamespace
import importlib.util
import json
import sys
import zipfile

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import compute_yellowness_profiles as yellow
import prepare_disease_gwas as disease
import run_embedding_replication as replication
import run_embedding_correlations as correlations
import run_gwas_panicle as gwas


@pytest.mark.parametrize('mode', ['missing', 'rejected', 'exception'])
def test_failed_yellowness_preserves_existing_output(tmp_path, monkeypatch, mode):
    image = tmp_path / 'leaf.jpg'
    if mode != 'missing':
        image.touch()
    meta = tmp_path / 'metadata.csv'
    pd.DataFrame([dict(environment='Nebraska2025', genotype='G1', image_id='leaf',
                       image_path=str(image))]).to_csv(meta, index=False)
    output = tmp_path / 'profiles.csv'
    output.write_text('previous result\n')
    monkeypatch.setattr(yellow, 'META', meta)
    monkeypatch.setattr(yellow, 'EXCLUDE_LIST', None)
    monkeypatch.setattr(sys, 'argv', ['profiles', '--out', str(output)])

    class SerialPool:
        def __init__(self, *args): pass
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def imap_unordered(self, worker, rows, **kwargs): return map(worker, rows)

    def profile(path):
        if mode == 'exception':
            raise ImportError('dependency unavailable')
        return None

    monkeypatch.setattr(yellow, 'Pool', SerialPool)
    monkeypatch.setattr(yellow, 'leaf_profile', profile)
    with pytest.raises((SystemExit, RuntimeError)) as error:
        yellow.main()
    assert 'leaf.jpg' in str(error.value) or 'No usable leaf profiles' in str(error.value)
    assert output.read_text() == 'previous result\n'


def test_unreadable_image_is_an_error(tmp_path):
    corrupt = tmp_path / 'corrupt.jpg'
    corrupt.write_text('not an image')
    with pytest.raises(RuntimeError, match='Could not read image'):
        yellow._worker(('G1', corrupt))


@pytest.mark.parametrize('relative', [False, True])
def test_disease_area_path_reaches_fitting(tmp_path, monkeypatch, relative):
    area = tmp_path / 'area.npz'
    np.savez(area, metadata_json=np.array(json.dumps({
        'columns': ['source_image_path', 'mask_pixels'], 'data': [['leaf.jpg', 100000]]})))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['disease', '--area-npz',
                        'area.npz' if relative else str(area), '--out-dir', str(tmp_path / 'out')])
    class ReachedFitting(Exception): pass
    def stop(args): raise ReachedFitting()
    monkeypatch.setattr(disease.cb, 'load_data', stop)
    with pytest.raises(ReachedFitting):
        disease.main()


def test_zero_hit_outputs_are_readable_by_hotspot_consumers(tmp_path, monkeypatch):
    markers = pd.DataFrame(dict(CHROM=[1], POS=[100], MARKER=['m'], REF=['A'], ALT=['C']))
    hits = tmp_path / 'hits.csv'
    hits.write_text('stale output')
    gwas.write_significant_markers([], markers, 0.05, hits)
    saved = pd.read_csv(hits)
    assert saved.empty
    assert {'trait', 'CHROM', 'POS', 'p_value', 'effect', 'se', 'q_value_within_trait'}.issubset(saved)
    replication.generate_hotspot_tables(hits, hits, tmp_path)
    peaks = tmp_path / 'sam3_peaks_ge10_embeddings.csv'
    assert pd.read_csv(peaks).empty
    assert replication.select_all_hotspot_embedding_pairs(peaks, hits).empty
    monkeypatch.setattr(correlations, 'SIGNIFICANT_MARKERS', {'sam3': hits, 'dino2': hits})
    assert correlations.hotspot_embeddings(pd.DataFrame([dict(
        peak_marker='m', chrom=1, peak_start_bp=0, peak_end_bp=100000)])) == {'m': []}


def test_hidden_ancestors_are_excluded(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location('review_server', ROOT / 'supporting_info/LeafWebScore/server.py')
    server = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(server)
    for name in ['visible/leaf.jpg', '.hidden/leaf.jpg', 'nested/__MACOSX/leaf.jpg',
                 'nested/.backup/leaf.jpg', 'nested/.leaf.jpg']:
        path = tmp_path / 'project' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    monkeypatch.setattr(server, 'IMAGES_DIR', tmp_path)
    assert server.discover_projects()['project']['images'] == ['visible/leaf.jpg']


@pytest.fixture
def blue_run(tmp_path, monkeypatch):
    for name in ['data/provided/field_image_metadata.csv', 'data/provided/image_ids_exclude.csv',
                 'scripts/calculate_blues.py', 'scripts/embedding_io.py', 'scripts/run_embedding_replication.py']:
        p = tmp_path / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('fixture')
    monkeypatch.setattr(replication, 'REPO_ROOT', tmp_path)
    monkeypatch.setattr(replication.subprocess, 'check_output', lambda *a, **kw: 'R and package versions')
    source = tmp_path / 'scores.npz'
    np.savez(source, features=np.array([[1.0]]))
    state = SimpleNamespace(calls=0, fail=False, source=source, output=tmp_path / 'blues')
    def fit(command, **kwargs):
        state.calls += 1
        if state.fail:
            raise RuntimeError('fit failed')
        out = Path(command[command.index('--out-dir') + 1])
        for env in replication.ENVIRONMENTS:
            pd.DataFrame({'genotype': ['G1'], 'embedding_mean_0': [1.0]}).to_csv(out / f'blues_{env}.csv', index=False)
    monkeypatch.setattr(replication.subprocess, 'run', fit)
    state.run = lambda: replication.calculate_blues(source, state.output, ['embedding_mean_0'], 1, True)
    return state


def test_blue_reuse_tracks_inputs_and_validates_outputs(blue_run, tmp_path):
    run = blue_run
    run.run()
    run.run()
    assert run.calls == 1
    # Changing the ZIP container timestamp alone must not invalidate the fit.
    with zipfile.ZipFile(run.source) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
    with zipfile.ZipFile(run.source, 'w') as archive:
        for name, value in members.items():
            archive.writestr(zipfile.ZipInfo(name, date_time=(2025, 1, 1, 0, 0, 0)), value)
    run.run()
    assert run.calls == 1
    np.savez(run.source, features=np.array([[2.0]]))
    run.run()
    assert run.calls == 2
    (tmp_path / 'data/provided/image_ids_exclude.csv').write_text('changed exclusions')
    run.run()
    assert run.calls == 3
    output = run.output / f'blues_{replication.ENVIRONMENTS[0]}.csv'
    output.write_text('genotype,embedding_mean_0\n')
    run.run()
    assert run.calls == 4
    output.write_text('genotype,embedding_mean_0\nG1,42\n')
    run.run()
    assert run.calls == 5
    (run.output / 'reuse_provenance.json').unlink()
    run.run()
    assert run.calls == 6


def test_failed_blue_refit_retains_previous_complete_fit(blue_run):
    run = blue_run
    run.run()
    before = {p.name: p.read_bytes() for p in run.output.iterdir()}
    np.savez(run.source, features=np.array([[3.0]]))
    run.fail = True
    with pytest.raises(RuntimeError, match='fit failed'):
        run.run()
    assert {p.name: p.read_bytes() for p in run.output.iterdir()} == before


@pytest.mark.parametrize('scope', ['within', 'cross'])
def test_empty_correlations_do_not_load_genotypes(tmp_path, monkeypatch, scope):
    hits = tmp_path / 'hits.csv'
    pd.DataFrame(columns=['CHROM', 'POS', 'trait']).to_csv(hits, index=False)
    master = tmp_path / 'hotspots.csv'
    pd.DataFrame(columns=['peak_marker', 'chrom', 'peak_start_bp', 'peak_end_bp']).to_csv(master, index=False)
    output = tmp_path / 'correlations.csv'
    output.write_text('stale result')
    monkeypatch.setattr(correlations, 'SIGNIFICANT_MARKERS', {'sam3': hits, 'dino2': hits})
    monkeypatch.setattr(correlations, 'HOTSPOT_MASTER', master)
    monkeypatch.setattr(correlations, 'OUT_CSV', output)
    getattr(correlations, 'run_' + scope)()
    saved = pd.read_csv(output)
    assert saved.empty and {'response_embedding', 'predictor_embedding', 'partial_r'}.issubset(saved)
