"""Boundary, candidate competition and end-to-end ROOT fixture checks."""
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import uproot

import trace_fei_stages as t


def candidate(event, candidate=0, **updates):
    row = dict(experiment=7, run=23, event=event, ncandidates=2, candidate=candidate,
               upsilon_id=0, sig_mode=0, tag_mode=23, mass=.493677, probability=.2,
               eecl=.5, missing_theta=1., sig_p=1.5, mbc=5.28, delta_e=.02,
               tag_is_signal=1., tag_mc_errors=np.nan, tag_mc_pdg=np.nan)
    row.update(updates)
    return row


def frame(rows):
    result = pd.DataFrame(rows, columns=list(candidate(0)))
    result['category'] = np.where(result.upsilon_id == 1, 0,
                                   np.where(result.tag_mode.isin([0, 1, 3, 4]), 1,
                                            np.where(result.tag_mode.isin([15, 16, 18, 19]), 2, 3)))
    result['entry'] = np.arange(len(result))
    return result


def write_root(path, data, misalign=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(path) as root:
        for tree, fields in t.FIELDS.items():
            mapping = dict(t.ID_BRANCHES, **fields)
            if tree == 'Btag':
                mapping.update(t.TRUTH)
            arrays = {branch: data[name].to_numpy(dtype=np.int64 if name in t.IDENTITY else np.float64)
                      for name, branch in mapping.items()}
            if misalign and tree == 'Btag' and len(data):
                arrays['__event__'] = arrays['__event__'] + 1
            root.mktree(tree, {key: value.dtype for key, value in arrays.items()})
            if len(data):
                root[tree].extend(arrays)


def test_strict_cut_boundaries():
    mass = frame([candidate(i, mass=m) for i, m in enumerate([1.84, 1.85, 1.89])])
    assert t.apply_next_cut(mass, 'before_Dveto_cut').mass.tolist() == [1.84, 1.89]
    assert t.apply_next_cut(frame([candidate(1, mbc=5.27)]), 'before_Mbc_cut').empty
    assert t.apply_next_cut(frame([candidate(1, delta_e=-.2)]), 'before_delE_cut').empty
    assert t.apply_next_cut(frame([candidate(1, sig_p=2.96)]), 'before_psig_cut').empty


def test_competing_other_mode_and_actual_random_winner():
    before = frame([
        candidate(100, probability=.1), candidate(100, 1, tag_mode=0, sig_mode=4, mass=1.2, probability=.8),
        candidate(101), candidate(101, 1, tag_mode=0, sig_mode=1, mass=.8),
        candidate(102, probability=.3),
        candidate(103), candidate(103, 1, tag_mode=0, sig_mode=4, mass=1.86, probability=.9),
    ])
    # Actual stored random choice at event 101 is the other category. We do not reseed/reselect.
    final = before.iloc[[1, 3, 4, 5]].copy()
    veto, maxima, issues = t.audit_bcs(before, final)
    assert issues == []
    rows = pd.DataFrame(t.event_trace(before, veto, maxima, final))
    mode0 = rows[(rows.origin_category == 3) & (rows.view == 'mode0')].set_index('event')
    assert mode0.loc[100, 'outcome'] == 'lost_at_max_probability'
    assert mode0.loc[100, 'winner_sig_mode'] == 4
    assert mode0.loc[101, 'outcome'] == 'lost_among_equal_maxima'
    assert mode0.loc[102, 'outcome'] == 'kept_target'
    assert mode0.loc[103, 'outcome'] == 'kept_target'  # high-score competitor removed by D veto
    assert mode0.loc[101, 'all_maximal_candidates'] == 2


def test_loader_grouping_and_nonmaximal_winner_are_not_hidden():
    before = frame([candidate(1, probability=.1), candidate(2, probability=.3),
                    candidate(1, 1, probability=.8)])
    _, maximum, issues = t.audit_bcs(before, before)
    assert len(maximum) == 3  # faithful contiguous grouping, not a global max over physical ID
    assert 'same_Loader_group_in_noncontiguous_blocks' in issues
    changed_n = frame([candidate(1, ncandidates=2), candidate(1, 1, probability=.9, ncandidates=3)])
    _, maximum, issues = t.audit_bcs(changed_n, changed_n)
    assert len(maximum) == 2
    assert 'multiple_ncandidates_values_for_one_physical_event' in issues
    same_group = frame([candidate(1), candidate(1, 1, probability=.8)])
    assert 'final_candidate_not_a_maximum_of_its_Loader_block' in t.audit_bcs(same_group, same_group.iloc[:1])[2]


def test_exact_ties_and_duplicate_candidates():
    before = frame([candidate(1), candidate(1, 1, probability=.2 + 1e-12)])
    assert t.max_probability_candidates(before).candidate.tolist() == [1]
    duplicate = pd.concat([before.iloc[:1], before.iloc[:1]], ignore_index=True)
    assert 'duplicate_candidate_identifiers_before_BCS' in t.audit_bcs(duplicate, duplicate.iloc[:1])[2]


def fixture(root):
    reference = root / 'reference'
    reference.mkdir(parents=True)
    version = 'validation_v009'
    manifest, events = [], []
    original = frame([
        candidate(100, probability=.1), candidate(100, 1, tag_mode=0, probability=.8),
        candidate(101), candidate(101, 1, tag_mode=0, sig_mode=1, mass=.8),
        candidate(102, probability=.3),
        candidate(103), candidate(103, 1, tag_mode=0, sig_mode=4, mass=1.86, probability=.9),
        candidate(104, mbc=5.27), candidate(105, delta_e=.2), candidate(106, eecl=1.3),
        candidate(107, missing_theta=.297), candidate(108, sig_p=.5),
    ])
    for fid, sample in enumerate(t.SAMPLES):
        base = root / ('data/SIGNAL_analysis' if sample == 'DATA' else f'mc/{sample}_analysis') / version
        stem = f'input_{sample}'
        data = original.copy() if sample in ['DATA', 'SIGNAL'] else original.iloc[:0].copy()
        for stage in t.SNAPSHOTS[:-1]:
            write_root(base / stage / (stem + '_' + stage + '.root'), data)
            data = t.apply_next_cut(data, stage)
        final = data[data.entry.isin([1, 3, 4, 5])].copy()
        write_root(base / 'final_output' / (stem + '_' + version + '_final_output.root'), final)
        # Merged files in subdirectories must not enter the cohort.
        (base / 'before_Dveto_cut' / 'Merge').mkdir()
        (base / 'before_Dveto_cut' / 'Merge' / 'merged.root').write_bytes(b'not a real file')
        norm = 1. if sample == 'DATA' else .25
        manifest.append(dict(file_id=fid, sample=sample,
                             path=f'/host/{sample}/MVAoutput_{stem}_{version}_final_output_data_after_MVA.root',
                             entries=len(final), normalization=norm, sync_keys='__event__;'))
        for entry, row in enumerate(final.to_dict('records')):
            events.append(dict(file_id=fid, entry=entry, sample=sample,
                               **{key: row[key] for key in t.EVENT + ['candidate', 'category', 'tag_mode', 'sig_mode']}))
    pd.DataFrame(manifest).to_csv(reference / 'files.csv', index=False)
    pd.DataFrame(events).to_csv(reference / 'selected_events.csv', index=False)
    args = Namespace(reference=reference, mc_base=root / 'mc', data_base=root / 'data',
                     output=root / 'result', scope='full', version=version, preflight_only=False)
    return args


def test_end_to_end_root_snapshots(tmp_path):
    args = fixture(tmp_path)
    assert t.run(args) == 0
    assert (args.output / 'COMPLETE').exists()
    checks = pd.read_csv(args.output / 'file_checks.csv')
    assert checks.selected_candidate_reproduction.all()
    stages = pd.read_csv(args.output / 'stage_summary.csv')
    data = stages[(stages['sample'] == 'DATA') & (stages.category == 3) & (stages.view == 'mode0')].set_index('stage')
    assert data.loc['before_Mbc_cut', 'events'] == 9
    assert data.loc['before_Dveto_cut', 'events'] == 4
    assert data.loc['final_output', 'events'] == 2
    traces = pd.read_csv(args.output / 'bcs_event_trace.csv.gz')
    target = traces[(traces['sample'] == 'DATA') & (traces.origin_category == 3) & (traces.view == 'mode0')]
    assert Counter_for_test(target.outcome) == {'kept_target': 2, 'lost_at_max_probability': 1, 'lost_among_equal_maxima': 1}
    ratios = pd.read_csv(args.output / 'data_mc_stage_ratios.csv')
    last = ratios[(ratios.stage == 'final_output') & (ratios.category == 3) & (ratios.view == 'mode0')].iloc[0]
    assert last.data == 2 and last.mc_normalization_only == .5


def Counter_for_test(values):
    return values.value_counts().to_dict()


def test_empty_samples_bcs_scope_and_preflight(tmp_path):
    args = fixture(tmp_path)
    args.scope = 'bcs'
    args.preflight_only = True
    assert t.run(args) == 0 and not args.output.exists()
    args.preflight_only = False
    assert t.run(args) == 0
    summary = pd.read_csv(args.output / 'stage_summary.csv')
    assert summary.loc[summary['sample'] == 'CHG', 'events'].sum() == 0
    assert 'before_Mbc_cut' not in set(summary.stage)


def test_alignment_and_cohort_errors(tmp_path):
    args = fixture(tmp_path)
    path = tmp_path / 'unaligned.root'
    write_root(path, frame([candidate(1)]), misalign=True)
    with pytest.raises(ValueError, match='alignment'):
        t.read_root(path)
    lost = args.data_base / 'SIGNAL_analysis/validation_v009/final_output/input_DATA_validation_v009_final_output.root'
    lost.unlink()
    with pytest.raises(ValueError, match='cohort differs'):
        t.run(args)
    assert not args.output.exists()


def test_cut_mismatch_requires_review(tmp_path):
    args = fixture(tmp_path)
    path = args.data_base / 'SIGNAL_analysis/validation_v009/before_delE_cut/input_DATA_before_delE_cut.root'
    wrong, _ = t.read_root(path)
    wrong.loc[wrong.index[0], 'mbc'] = 5.29  # identities match; payload does not
    write_root(path, wrong)
    assert t.run(args) == 2
    assert (args.output / 'REVIEW_REQUIRED').exists()
    assert not (args.output / 'COMPLETE').exists()
    issues = pd.read_csv(args.output / 'audit_issues.csv')
    assert 'stored_snapshot_differs_from_supplied_cut_code' in set(issues.issue)


def test_parallel_reads_match_serial_outputs(tmp_path):
    args = fixture(tmp_path)
    args.scope = 'bcs'
    assert t.run(args) == 0
    serial = args.output
    args.output = tmp_path / 'parallel'
    args.workers = 2
    assert t.run(args) == 0
    for name in ['stage_summary.csv', 'data_mc_stage_ratios.csv', 'bcs_migration_summary.csv',
                 'bcs_event_trace.csv.gz', 'bcs_candidate_trace.csv.gz', 'file_checks.csv',
                 'source_files.csv', 'branch_availability.csv']:
        pd.testing.assert_frame_equal(pd.read_csv(serial / name), pd.read_csv(args.output / name))


if __name__ == '__main__':
    raise SystemExit(pytest.main([__file__, '-q']))
