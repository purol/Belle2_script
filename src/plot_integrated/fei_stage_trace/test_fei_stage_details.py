from pathlib import Path
import hashlib

import numpy as np
import pandas as pd
import pytest
import uproot

import trace_fei_stages as t
from fei_stage_details import describe
from test_trace_fei_stages import candidate, fixture, frame, write_root


def detailed_fixture(root):
    args = fixture(root)
    args.details = True
    extra = frame([
        candidate(200, tag_mode=15, probability=.4),
        candidate(200, 1, tag_mode=0, probability=.9, mbc=5.27),
        candidate(201, tag_mode=16, probability=.8, delta_e=.2),
        candidate(202, tag_mode=18, probability=.2, tag_is_signal=0),
        candidate(202, 1, tag_mode=3, probability=.7),
        candidate(203, tag_mode=19, probability=.6),
        candidate(204, tag_mode=22, probability=.3, tag_is_signal=np.nan),
        candidate(204, 1, tag_mode=15, probability=.8, eecl=1.3),
    ])
    refs = pd.read_csv(args.reference / 'selected_events.csv')
    manifest = pd.read_csv(args.reference / 'files.csv')
    additions = []
    for sample in ['DATA', 'SIGNAL']:
        base = (args.data_base / 'SIGNAL_analysis' if sample == 'DATA' else args.mc_base / (sample + '_analysis')) / args.version
        stem = 'input_' + sample
        original, _ = t.read_root(base / t.SNAPSHOTS[0] / (stem + '_' + t.SNAPSHOTS[0] + '.root'))
        data = pd.concat([original, extra], ignore_index=True)
        data['entry'] = np.arange(len(data))
        for stage in t.SNAPSHOTS[:-1]:
            write_root(base / stage / (stem + '_' + stage + '.root'), data)
            data = t.apply_next_cut(data, stage)
        # Preserve old recorded winners and explicitly choose known new winners.
        old = data[data.event.lt(200) & data.entry.isin([1, 3, 4, 5])]
        new = data[((data.event == 200) & (data.candidate == 0)) |
                   ((data.event == 202) & (data.candidate == 1)) |
                   (data.event.isin([203, 204]) & (data.candidate == 0))]
        final = pd.concat([old, new], ignore_index=True)
        write_root(base / 'final_output' / (stem + '_' + args.version + '_final_output.root'), final)
        fid = manifest.loc[manifest['sample'] == sample, 'file_id'].iloc[0]
        manifest.loc[manifest.file_id == fid, 'entries'] = len(final)
        for entry, row in enumerate(final.to_dict('records')):
            if row['event'] >= 200:
                additions.append(dict(file_id=fid, entry=entry, sample=sample,
                                      **{k: row[k] for k in t.EVENT + ['candidate', 'category', 'tag_mode', 'sig_mode']}))
    pd.concat([refs, pd.DataFrame(additions)], ignore_index=True).to_csv(args.reference / 'selected_events.csv', index=False)
    manifest.to_csv(args.reference / 'files.csv', index=False)
    return args


def test_truth_ties_and_other_category_are_distinct():
    g = frame([candidate(1, tag_mode=22, probability=.3, tag_is_signal=np.nan),
               candidate(1, 1, tag_mode=24, probability=.3, tag_is_signal=0),
               candidate(1, 2, tag_mode=23, sig_mode=4, mass=1.2, probability=.8, tag_is_signal=1)])
    row = describe(g, 3, 'mode0', 'SIGNAL')
    assert row['best_target_modes'] == '22;24'
    assert row['truth_class'] == 'unknown_without_correct'
    assert row['n_target'] == 2 and row['n_outside_target'] == 1
    assert not row['other_category'] and not row['target_at_event_maximum']
    assert describe(g, 3, 'mode0', 'DATA')['truth_class'] == 'data'


def test_optional_fields_are_read_or_marked_missing(tmp_path):
    data = frame([candidate(1)])
    path = tmp_path / 'optional.root'
    with uproot.recreate(path) as root:
        for tree, fields in t.FIELDS.items():
            mapping = dict(t.ID_BRANCHES, **fields)
            if tree == 'Btag':
                mapping.update(t.TRUTH)
            arrays = {branch: data[name].to_numpy() for name, branch in mapping.items()}
            if tree == 'Upsilon':
                arrays[t.EXTRA_FIELDS['Upsilon']['eecl_matched']] = np.array([.37])
            root.mktree(tree, {key: value.dtype for key, value in arrays.items()})
            root[tree].extend(arrays)
    actual, _ = t.read_root(path, details=True, capture_schema=True)
    assert actual.eecl_matched.iloc[0] == .37
    assert np.isnan(actual.eecl_unmatched.iloc[0])
    assert actual.attrs['extra_availability']['eecl_matched']
    assert not actual.attrs['extra_availability']['eecl_unmatched']
    assert t.EXTRA_FIELDS['Upsilon']['eecl_matched'] in actual.attrs['branch_schema']['Upsilon']


def test_full_detailed_lineage_and_component_closure(tmp_path):
    args = detailed_fixture(tmp_path)
    input_paths = list(tmp_path.rglob('*.root')) + list(args.reference.glob('*.csv'))
    input_hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in input_paths}
    assert t.run(args) == 0
    assert {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in input_paths} == input_hashes
    assert (args.output / 'COMPLETE').exists()
    rows = pd.read_csv(args.output / 'selection_transitions.csv.gz')
    d = rows[(rows['sample'] == 'DATA') & (rows['view'] == 'mode0')]
    # Removing the higher-scoring competitor at Mbc keeps the C2 target.
    mbc = d[(d.event == 200) & (d.category == 2) & (d.stage_before == 'before_Mbc_cut')].iloc[0]
    assert mbc.before_other_category and not mbc.after_other_category and mbc.retained_target
    delta = d[(d.event == 201) & (d.category == 2) & (d.stage_before == 'before_delE_cut')].iloc[0]
    assert not delta.retained_target and not delta.retained_event
    eecl = d[(d.event == 204) & (d.category == 2) & (d.stage_before == 'before_Eecl_cut')].iloc[0]
    assert not eecl.retained_target and eecl.retained_event
    bcs = d[(d.event == 202) & (d.category == 2) & (d.stage_before == 'after_Dveto_reconstructed')].iloc[0]
    assert not bcs.retained_target and bcs.retained_event
    # An event that later loses the target still has its surviving competitor exported.
    candidates = pd.read_csv(args.output / 'selection_candidate_trace.csv.gz')
    last = candidates[(candidates['sample'] == 'DATA') & (candidates.event == 204) & (candidates.stage == 'final_output')]
    assert last.category.tolist() == [3]
    assert last.matches_final_identity.all()
    events = pd.read_csv(args.output / 'selection_event_trace.csv.gz')
    assert not events.duplicated(['file_id', 'stage'] + t.EVENT + ['category', 'view']).any()
    flow = pd.read_csv(args.output / 'category23_cutflow.csv')
    target = flow[(flow.category == 2) & (flow['view'] == 'mode0')].set_index('stage')
    assert target.loc['before_Mbc_cut', 'data'] == 5
    assert target.loc['before_Dveto_cut', 'data'] == 3
    assert target.loc['final_output', 'data'] == 2
    coverage = pd.read_csv(args.output / 'extra_branch_coverage.csv')
    assert coverage.files_available.sum() == 0
    assert (coverage.files_missing == 1).all()


def test_detailed_parallel_and_serial_agree(tmp_path):
    args = detailed_fixture(tmp_path)
    assert t.run(args) == 0
    serial = args.output
    args.output = tmp_path / 'parallel_details'
    args.workers = 2
    assert t.run(args) == 0
    for name in ['selection_components.csv', 'selection_event_trace.csv.gz',
                 'selection_transitions.csv.gz', 'selection_candidate_trace.csv.gz',
                 'category23_cutflow.csv', 'category23_cut_retention.csv', 'extra_branch_coverage.csv']:
        pd.testing.assert_frame_equal(pd.read_csv(serial / name), pd.read_csv(args.output / name))


def test_details_require_full_scope(tmp_path):
    args = fixture(tmp_path)
    args.scope = 'bcs'
    args.details = True
    with pytest.raises(ValueError, match='requires --scope full'):
        t.run(args)
    assert not args.output.exists()
