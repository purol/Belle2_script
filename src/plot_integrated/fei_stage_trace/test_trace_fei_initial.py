from argparse import Namespace
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import uproot

import trace_fei_stages as base
import trace_fei_initial as initial
from test_trace_fei_stages import candidate, frame


def rows_for(sample):
    # A mode-30 Data event fails only loose tracks; the corresponding MC event survives.
    rows = [candidate(100, tag_mode=30), candidate(101, tag_mode=31),
            candidate(102, tag_mode=32), candidate(103, tag_mode=19),
            candidate(104, tag_mode=18), candidate(105, tag_mode=16),
            candidate(106, tag_mode=15), candidate(107, tag_mode=30, sig_mode=4, mass=1.4),
            candidate(108, tag_mode=22), candidate(108, 1, tag_mode=0, probability=.9)]
    data = frame(rows)
    for name in initial.INITIAL_NAMES:
        data[name] = .4 if name.endswith('vertex') else 0.
    data.loc[data.event.eq(100), 'roe_tracks_loose'] = 1. if sample == 'DATA' else 0.
    data.loc[data.event.eq(101), 'roe_pi0'] = 1.
    data.loc[data.event.eq(102), 'roe_ks'] = 1.
    data.loc[data.event.eq(103), 'sig_vertex'] = -.5
    data.loc[data.event.eq(104), 'tag_vertex'] = -.5
    data.loc[data.event.eq(105), 'sig_vertex'] = np.nan
    data.loc[data.event.eq(105), 'mass'] = np.nan  # Preserve rather than replace with zero.
    data.loc[data.event.eq(108) & data.candidate.eq(1), 'roe_tracks_loose'] = 1.
    return data


def write_root(path, data, missing_vertex=False, misalign=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(path) as root:
        for tree, fields in base.FIELDS.items():
            mapping = dict(base.ID_BRANCHES, **fields, **initial.INITIAL_FIELDS[tree])
            if tree == 'Btag':
                mapping.update(base.TRUTH)
            if missing_vertex and tree == 'Bsig':
                del mapping['sig_vertex']
            arrays = {branch:data[name].to_numpy(dtype=np.int64 if name in base.IDENTITY else np.float64)
                      for name,branch in mapping.items()}
            if misalign and tree == 'Btag' and len(data):
                arrays['__event__'] = arrays['__event__'] + 1
            root.mktree(tree,{name:values.dtype for name,values in arrays.items()})
            if len(data):
                root[tree].extend(arrays)


def fixture(root):
    reference = root/'reference'
    reference.mkdir()
    args = Namespace(reference=reference,raw_mc_base=root/'raw_mc',raw_data_base=root/'raw_data',
                     mc_base=root/'analysis_mc',data_base=root/'analysis_data',version='validation_v009',
                     output=root/'result',workers=1,preflight_only=False)
    manifest, anchors = [], []
    for fid,sample in enumerate(base.SAMPLES):
        sample_dir = 'SIGNAL' if sample=='DATA' else sample
        stem = 'test_'+sample
        data = rows_for(sample)
        if sample not in ['DATA','SIGNAL']:
            data = data.iloc[:0]
        raw_path = (args.raw_data_base if sample=='DATA' else args.raw_mc_base)/sample_dir/'validation'/(stem+'.root')
        anchor_path = (args.data_base if sample=='DATA' else args.mc_base)/(sample_dir+'_analysis')/args.version/'before_Mbc_cut'/(stem+'_before_Mbc_cut.root')
        write_root(raw_path,data)
        saved = initial.initial_states(data)[0]['before_Mbc_cut']
        write_root(anchor_path,saved)
        norm = 1. if sample=='DATA' else .25
        manifest.append(dict(file_id=fid,sample=sample,path='/old/MVAoutput_'+stem+'_'+args.version+'_final_output_data_after_MVA.root',normalization=norm))
        for category in range(4):
            for view in base.VIEWS:
                selected = saved[saved.category.eq(category)&base.view_mask(saved,view)]
                n = len(selected.drop_duplicates(base.EVENT))
                anchors.append(dict(sample=sample,stage='before_Mbc_cut',category=category,view=view,
                                    events=n,sumw=n*norm,sumw2=n*norm**2,candidates=len(selected)))
    pd.DataFrame(manifest).to_csv(reference/'files.csv',index=False)
    pd.DataFrame(anchors).to_csv(reference/'before_mbc_reference.csv',index=False)
    # Deliberately invalid companion/merge files must not be opened.
    companion = args.raw_data_base/'Jpsi'/'SIGNAL'
    companion.mkdir(parents=True)
    (companion/'test_DATA.root').write_bytes(b'not the Analysis_main input')
    merge = args.raw_data_base/'SIGNAL'/'validation'/'Merge'
    merge.mkdir()
    (merge/'merged.root').write_bytes(b'not a separate source')
    return args


def test_strict_boundaries_and_nan_follow_cpp():
    data = rows_for('DATA')
    data.loc[data.event.eq(100),'roe_tracks_loose'] = .5
    states, failure = initial.initial_states(data)
    assert failure[:7].tolist() == [1,2,3,4,5,4,0]
    assert states['before_Mbc_cut'].event.tolist() == [106,107,108]
    # A competitor fails, while the C3 target in its event is retained.
    assert failure[-2:].tolist() == [0,1]


def test_reader_keeps_nan_and_checks_alignment(tmp_path):
    path=tmp_path/'raw.root';write_root(path,rows_for('DATA'))
    data,availability,schema=initial.read_input(path,True)
    assert np.isnan(data.loc[data.event.eq(105),'sig_vertex'].iloc[0])
    assert np.isnan(data.loc[data.event.eq(105),'mass'].iloc[0])
    assert not availability['n_jpsi']
    assert 'Btag_chiProb' in schema['Btag']
    write_root(path,rows_for('DATA'),missing_vertex=True)
    with pytest.raises(ValueError,match='missing required'):
        initial.read_input(path)
    write_root(path,rows_for('DATA'),misalign=True)
    with pytest.raises(ValueError,match='alignment'):
        initial.read_input(path)


def test_end_to_end_initial_loss_and_unchanged_inputs(tmp_path):
    args=fixture(tmp_path)
    paths=list(tmp_path.rglob('*.root'))+list(args.reference.glob('*.csv'))
    before={path:hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    assert initial.run(args)==0
    assert (args.output/'COMPLETE').exists()
    assert before=={path:hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    checks=pd.read_csv(args.output/'file_checks.csv')
    assert checks.initial_payload_reproduced.all() and checks.issues.sum()==0
    focus=pd.read_csv(args.output/'focus_mode_groups.csv')
    d=focus[focus['sample'].eq('DATA')&focus.mode_group.eq('category3_modes30_32')]
    assert d[d.stage.eq('raw_input')&d['view'].eq('all_modes')].events.iloc[0]==4
    assert d[d.stage.eq('raw_input')&d['view'].eq('control_K_Kpi')].events.iloc[0]==3
    assert d[d.stage.eq('before_Mbc_cut')&d['view'].eq('control_K_Kpi')].events.iloc[0]==0
    assert d[d.stage.eq('before_Mbc_cut')&d['view'].eq('all_modes')].events.iloc[0]==1
    losses=pd.read_csv(args.output/'initial_event_losses.csv.gz')
    lost=losses[losses['sample'].eq('DATA')&losses.event.eq(100)&losses['view'].eq('mode0')].iloc[0]
    assert lost.first_target_loss_stage=='after_loose_tracks' and not lost.retained_before_mbc
    candidates=pd.read_csv(args.output/'initial_candidate_trace.csv.gz')
    comp=candidates[candidates['sample'].eq('DATA')&candidates.event.eq(108)&candidates.candidate.eq(1)].iloc[0]
    assert comp.first_failed_index==1 and not comp.in_saved_before_mbc
    outside=candidates[candidates['sample'].eq('DATA')&candidates.event.eq(107)].iloc[0]
    assert outside.sig_mode==4 and outside.mass>=1 and outside.in_saved_before_mbc
    assert not (args.output/'RUNNING').exists()
    with pytest.raises(ValueError,match='Output exists'):
        initial.run(args)


def test_parallel_results_match_serial(tmp_path):
    args=fixture(tmp_path)
    assert initial.run(args)==0
    serial=args.output
    args.output=tmp_path/'parallel';args.workers=2
    assert initial.run(args)==0
    for name in ['initial_stage_summary.csv','tag_mode_presence.csv','focus_mode_groups.csv',
                 'tag_signal_mass_composition.csv','initial_candidate_trace.csv.gz',
                 'initial_event_trace.csv.gz','initial_event_losses.csv.gz','file_checks.csv',
                 'initial_data_mc_ratios.csv']:
        pd.testing.assert_frame_equal(pd.read_csv(serial/name),pd.read_csv(args.output/name))


def test_wrong_anchor_is_reported_not_completed(tmp_path):
    args=fixture(tmp_path);jobs=initial.preflight(args)
    path=Path(jobs[0]['anchor_path']);saved,_,_=initial.read_input(path)
    saved.loc[saved.index[0],'mbc']=5.29
    write_root(path,saved)
    assert initial.run(args)==2
    assert (args.output/'REVIEW_REQUIRED').exists() and not (args.output/'COMPLETE').exists()
    issues=pd.read_csv(args.output/'audit_issues.csv')
    assert 'initial_cuts_do_not_reproduce_saved_before_Mbc_payload' in set(issues.issue)


def test_cohort_and_anchor_reference_validated_before_scan(tmp_path):
    args=fixture(tmp_path)
    args.preflight_only=True
    assert initial.run(args)==0 and not args.output.exists()
    extra=args.raw_data_base/'SIGNAL'/'validation'/'unexpected.root'
    extra.write_bytes(b'extra source')
    with pytest.raises(ValueError,match='cohort mismatch'):
        initial.run(args)
    assert not args.output.exists()


def test_incomplete_reference_rejected(tmp_path):
    args=fixture(tmp_path)
    path=args.reference/'before_mbc_reference.csv'
    pd.read_csv(path).iloc[1:].to_csv(path,index=False)
    with pytest.raises(ValueError,match='all sample/category/view'):
        initial.run(args)
    assert not args.output.exists()
