"""Audit the five initial cuts using raw ntuples and saved before_Mbc snapshots.

Read-only inputs. No FEI reconstruction or BCS is rerun. All samples/files from
the existing FEI manifest are retained, including files with no final events.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing
from pathlib import Path, PurePosixPath
import sys
import time

import numpy as np
import pandas as pd
import uproot

import trace_fei_stages as base
from fei_stage_details import describe

HERE = Path(__file__).resolve().parent
RAW_MC = '/home/belle2/junewoo/storage_ghi/Ntuple/MizuhashiRD_LS_MC_Jpsi'
RAW_DATA = '/home/belle2/junewoo/storage_ghi/Ntuple/MizuhashiRD_LS_data_Jpsi'
INITIAL_FIELDS = {
    'Upsilon': {'roe_tracks_loose': 'nROE_Tracks__bolooseMask__bc',
                'roe_pi0': 'nROE_ParticlesInList__bopi0__clmyneutralPion__bc',
                'roe_ks': 'nROE_ParticlesInList__boK_S0__clmyKaonshort__bc'},
    'Bsig': {'sig_vertex': 'Bsig_daughter_0_extraInfo_mychiProb'},
    'Btag': {'tag_vertex': 'Btag_chiProb'},
}
CUTS = [('roe_tracks_loose', 'lt', .5, 'after_loose_tracks'),
        ('roe_pi0', 'lt', .5, 'after_roe_pi0'),
        ('roe_ks', 'lt', .5, 'after_roe_ks'),
        ('sig_vertex', 'gt', -.5, 'after_signal_vertex'),
        ('tag_vertex', 'gt', -.5, 'before_Mbc_cut')]
STAGES = ['raw_input'] + [cut[3] for cut in CUTS]
INITIAL_NAMES = [name for fields in INITIAL_FIELDS.values() for name in fields]
COMPARE = base.PAYLOAD + INITIAL_NAMES
FOCUS = {'category3_modes30_32': (3, [30, 31, 32]),
         'category3_mode22': (3, [22]), 'category2_modes18_19': (2, [18, 19])}


def read_input(path, capture_schema=False):
    """Preserve NaNs before vertex cuts; require aligned scalar IDs and mode IDs."""
    columns, available, schema = {}, {}, {}
    with Path(path).open('rb') as handle, uproot.open(handle, object_cache=None, array_cache=None) as root:
        entries = None
        for tree_name, fields in base.FIELDS.items():
            base.require(tree_name in root, f'{path}: missing tree {tree_name}')
            tree = root[tree_name]
            names = set(tree.keys())
            if capture_schema:
                schema[tree_name] = sorted(names)
            if entries is None:
                entries = tree.num_entries
            base.require(tree.num_entries == entries, f'{path}: unequal tree entry counts')
            required = dict(base.ID_BRANCHES, **fields, **INITIAL_FIELDS[tree_name])
            missing = set(required.values()) - names
            base.require(not missing, f'{path}: {tree_name} missing required branches {sorted(missing)}')
            optional = dict(base.EXTRA_FIELDS.get(tree_name, {}))
            if tree_name == 'Btag':
                optional.update(base.TRUTH)
            mapping = dict(required)
            for name, branch in optional.items():
                available[name] = branch in names
                if branch in names:
                    mapping[name] = branch
            arrays = (tree.arrays(list(set(mapping.values())), library='np', how=dict)
                      if entries else {branch: np.array([], dtype=float) for branch in mapping.values()})
            for name, branch in mapping.items():
                values = np.asarray(arrays[branch])
                base.require(values.ndim == 1 and values.dtype.kind in 'biuf',
                             f'{path}: non-scalar/non-numeric branch {branch}')
                if name in base.ID_BRANCHES:
                    base.require(np.isfinite(values).all() and (values == np.round(values)).all(),
                                 f'{path}: invalid identifier {branch}')
                    values = values.astype(np.int64)
                if name in columns:
                    base.require(np.array_equal(values, columns[name]), f'{path}: tree alignment differs for {branch}')
                else:
                    columns[name] = values
            for name, branch in optional.items():
                if name not in columns:
                    columns[name] = np.full(entries, np.nan)
    frame = pd.DataFrame(columns)
    for name in ['upsilon_id', 'sig_mode', 'tag_mode']:
        base.require(np.isfinite(frame[name]).all() and frame[name].eq(frame[name].round()).all(),
                     f'{path}: invalid mode identifier {name}')
        frame[name] = frame[name].astype(np.int64)
    base.require(frame.upsilon_id.isin([0, 1]).all() and frame.tag_mode.ge(0).all(),
                 f'{path}: invalid Upsilon or tag mode')
    base.require((frame.probability.isna() | frame.probability.between(0, 1)).all(),
                 f'{path}: FEI probability outside [0,1]')
    frame['category'] = np.where(frame.upsilon_id.eq(1), 0,
        np.where(frame.tag_mode.isin([0, 1, 3, 4]), 1,
                 np.where(frame.tag_mode.isin([15, 16, 18, 19]), 2, 3)))
    frame['entry'] = np.arange(len(frame))
    return frame, available, schema


def initial_states(raw):
    active = np.ones(len(raw), dtype=bool)
    first_failure = np.zeros(len(raw), dtype=np.int16)
    states = {'raw_input': raw}
    for index, (field, operation, boundary, stage) in enumerate(CUTS, 1):
        # Identical strict comparisons to Loader::Cut; NaN comparisons fail.
        passed = (raw[field].lt(boundary) if operation == 'lt' else raw[field].gt(boundary)).to_numpy()
        first_failure[active & ~passed] = index
        active &= passed
        states[stage] = raw.loc[active].copy()
    return states, first_failure


def signatures(frame):
    values = frame[COMPARE].astype(object)
    values = values.where(pd.notna(values), None)
    return Counter(values.itertuples(index=False, name=None))


def preflight(args):
    manifest = pd.read_csv(args.reference / 'files.csv')
    base.require(set(manifest['sample']) == set(base.SAMPLES), 'Reference must contain all eight samples')
    base.require(not manifest.file_id.duplicated().any() and not manifest.path.duplicated().any(), 'Duplicate reference file')
    base.require(np.isfinite(manifest.normalization).all() and manifest.normalization.gt(0).all(), 'Invalid normalization')
    base.require(manifest.loc[manifest['sample'].eq('DATA'), 'normalization'].eq(1).all(), 'Data weight must be one')
    anchor = pd.read_csv(args.reference / 'before_mbc_reference.csv')
    key_columns = ['sample', 'category', 'view']
    expected_keys = {(sample, category, view) for sample in base.SAMPLES for category in range(4) for view in base.VIEWS}
    base.require(not anchor.duplicated(key_columns).any() and
                 set(anchor[key_columns].itertuples(index=False, name=None)) == expected_keys,
                 'before_mbc_reference.csv must contain all sample/category/view cells exactly once')
    base.require(anchor.stage.eq('before_Mbc_cut').all(), 'Anchor reference has the wrong stage')
    provenance_path = args.reference / 'initial_reference_provenance.json'
    if provenance_path.exists():
        provenance = json.loads(provenance_path.read_text())
        for name, expected_hash in provenance['reference_sha256'].items():
            base.require(hashlib.sha256((args.reference/name).read_bytes()).hexdigest() == expected_hash,
                         'Packaged reference hash differs: ' + name)
    jobs = []
    for sample in base.SAMPLES:
        selected = manifest[manifest['sample'].eq(sample)].copy()
        selected['stem'] = [base.source_stem(PurePosixPath(p).name, 'reference', args.version) for p in selected.path]
        base.require(not selected.stem.duplicated().any(), 'Duplicate source stem for ' + sample)
        expected = set(selected.stem)
        sample_dir = 'SIGNAL' if sample == 'DATA' else sample
        raw_dir = (args.raw_data_base if sample == 'DATA' else args.raw_mc_base) / sample_dir / 'validation'
        analysis_dir = (args.data_base if sample == 'DATA' else args.mc_base) / (sample_dir + '_analysis') / args.version / 'before_Mbc_cut'
        for directory in [raw_dir, analysis_dir]:
            base.require(directory.is_dir(), 'Missing directory: ' + str(directory))
        # Only the specified validation directory, never Jpsi companions/logs/recursive merges.
        raw_map = {p.stem: p for p in sorted(raw_dir.glob('*.root'))}
        anchor_map = {base.source_stem(p.name, 'before_Mbc_cut', args.version): p
                      for p in sorted(analysis_dir.glob('*.root'))}
        for label, paths in [('raw validation', raw_map), ('before_Mbc', anchor_map)]:
            missing, extra = expected - set(paths), set(paths) - expected
            base.require(not missing and not extra,
                f'{sample} {label}: cohort mismatch; missing={len(missing)} {sorted(missing)[:3]}, '
                f'extra={len(extra)} {sorted(extra)[:3]}. No unmatched files are silently dropped.')
        for index, row in enumerate(selected.to_dict('records')):
            jobs.append(dict(row, raw_path=str(raw_map[row['stem']]), anchor_path=str(anchor_map[row['stem']]),
                             capture_schema=index == 0))
    return jobs


def tally(frame):
    return [len(frame.drop_duplicates(base.EVENT)), len(frame)]


def process_job(job):
    raw, raw_available, raw_schema = read_input(job['raw_path'], job['capture_schema'])
    anchor, anchor_available, anchor_schema = read_input(job['anchor_path'], job['capture_schema'])
    frames, first_failure = initial_states(raw)
    reproduced = frames['before_Mbc_cut']
    issues = []
    for name, frame in [('raw_input', raw), ('before_Mbc_saved', anchor)]:
        if frame.duplicated(base.IDENTITY).any():
            issues.append(name + ': duplicate candidate identifiers')
    same = signatures(reproduced) == signatures(anchor)
    if not same:
        issues.append('initial_cuts_do_not_reproduce_saved_before_Mbc_payload')
    core_fields = [name for fields in base.FIELDS.values() for name in fields]
    if not np.isfinite(reproduced[core_fields].to_numpy()).all():
        issues.append('nonfinite_core_observables_survive_initial_cuts')
    counts, mode_counts, focus_counts, signal_counts = {}, {}, {}, {}
    for stage, frame in frames.items():
        for view in base.VIEWS:
            selected = frame.loc[base.view_mask(frame, view)]
            for category in range(4):
                counts[stage, category, view] = tally(selected[selected.category.eq(category)])
            for (category, mode), group in selected.groupby(['category', 'tag_mode'], sort=False):
                mode_counts[stage, int(category), view, int(mode)] = tally(group)
            for label, (category, wanted) in FOCUS.items():
                focus_counts[stage, view, label] = tally(selected[selected.category.eq(category) & selected.tag_mode.isin(wanted)])
        if not frame.empty:
            bins = np.where(~np.isfinite(frame.mass), 'nonfinite', np.where(frame.mass.lt(1), 'mass_lt_1', 'mass_ge_1'))
            grouped = frame.assign(mass_region=bins).groupby(['category', 'tag_mode', 'sig_mode', 'mass_region'], sort=False)
            for (category, mode, sig_mode, region), group in grouped:
                signal_counts[stage, int(category), int(mode), int(sig_mode), region] = tally(group)
    target_mask = raw.category.isin([2, 3]) & base.view_mask(raw, 'control_K_Kpi')
    # Keep mode-30/31/32 events even outside the control selection, so absence
    # in the control sample is not confused with absence in raw production.
    focus_mask = raw.category.eq(3) & raw.tag_mode.isin([30, 31, 32])
    targets = raw[target_mask | focus_mask]
    origins = set(targets[base.EVENT].itertuples(index=False, name=None))
    prefix = dict(file_id=job['file_id'], sample=job['sample'], normalization=job['normalization'], file_audit_ok=not issues)
    anchor_ids = set(anchor[base.IDENTITY].itertuples(index=False, name=None))
    candidate_rows = []
    for row in raw.to_dict('records'):
        if tuple(row[k] for k in base.EVENT) in origins:
            failed = int(first_failure[row['entry']])
            candidate_rows.append(dict(prefix, **row, first_failed_index=failed,
                first_failed_stage=STAGES[failed] if failed else 'passed_all_initial',
                in_saved_before_mbc=tuple(row[k] for k in base.IDENTITY) in anchor_ids))
    event_rows, initial_cohorts, retained_at = [], {}, defaultdict(list)
    for stage in STAGES:
        for key, group in frames[stage].groupby(base.EVENT, sort=False):
            if key not in origins:
                continue
            for category in [2, 3]:
                for view in ['control_K_Kpi', 'mode0']:
                    description = describe(group, category, view, job['sample'])
                    if description is None:
                        continue
                    description['probability_available'] = bool(np.isfinite(description['max_target_probability']))
                    if not description['probability_available']:
                        description['target_at_event_maximum'] = None
                    record = dict(prefix, stage=stage, stage_order=STAGES.index(stage),
                                  **dict(zip(base.EVENT, (int(x) for x in key))), **description)
                    event_rows.append(record)
                    cohort = (key, category, view)
                    retained_at[cohort].append(stage)
                    if stage == 'raw_input':
                        initial_cohorts[cohort] = description
    losses = []
    for (key, category, view), description in initial_cohorts.items():
        present = retained_at[key, category, view]
        missing = next((stage for stage in STAGES if stage not in present), 'passed_all_initial')
        losses.append(dict(prefix, **dict(zip(base.EVENT, (int(x) for x in key))), category=category, view=view,
            raw_target_modes=description['target_modes'], raw_best_target_modes=description['best_target_modes'],
            raw_truth_class=description['truth_class'], raw_other_category=description['other_category'],
            first_target_loss_stage=missing, retained_before_mbc=missing == 'passed_all_initial'))
    return dict(counts=counts, modes=mode_counts, focus=focus_counts, signal=signal_counts,
                events=event_rows, candidates=candidate_rows, losses=losses, issues=issues,
                reproduced=same, raw_rows=len(raw), derived_rows=len(reproduced), anchor_rows=len(anchor),
                data_events=list(set(raw[base.EVENT].itertuples(index=False, name=None))) if job['sample'] == 'DATA' else [],
                availability={'raw_input': raw_available, 'before_Mbc_saved': anchor_available},
                schema={'raw_input': raw_schema, 'before_Mbc_saved': anchor_schema},
                nonfinite={name:int((~np.isfinite(raw[name])).sum()) for name in core_fields + INITIAL_NAMES})


def results(jobs, workers):
    if workers == 1:
        for job in jobs:
            yield job, process_job(job)
        return
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as pool:
        iterator, pending = iter(jobs), deque()
        # One pending result per worker limits memory for larger raw ntuples.
        for _ in range(workers):
            job = next(iterator, None)
            if job is not None:
                pending.append((job, pool.submit(process_job, job)))
        while pending:
            job, future = pending.popleft()
            result = future.result()
            following = next(iterator, None)
            if following is not None:
                pending.append((following, pool.submit(process_job, following)))
            yield job, result


def run(args):
    base.require(args.workers >= 1, 'workers must be positive')
    output = args.output.resolve()
    base.require(not output.exists(), 'Output exists; choose a new directory: ' + str(output))
    jobs = preflight(args)
    print(f'Preflight passed: {len(jobs)} matched raw/anchor pairs; workers={args.workers}.', flush=True)
    if args.preflight_only:
        return 0
    output.mkdir(parents=True)
    (output / 'RUNNING').write_text('Incomplete until COMPLETE.\n')
    names = {'candidates':'initial_candidate_trace.csv.gz', 'events':'initial_event_trace.csv.gz',
             'losses':'initial_event_losses.csv.gz', 'checks':'file_checks.csv', 'issues':'audit_issues.csv',
             'manifest':'source_files.csv', 'schema':'representative_branch_inventory.csv',
             'nonfinite':'raw_nonfinite_counts.csv'}
    streams = {name:base.StreamCSV(output / path) for name,path in names.items()}
    aggregate = {name:defaultdict(lambda:np.zeros(4)) for name in ['counts','modes','focus','signal']}
    coverage, data_seen = defaultdict(lambda:np.zeros(2,dtype=int)), set()
    issue_count, start = 0, time.monotonic()
    try:
        for i, (job, result) in enumerate(results(jobs, args.workers), 1):
            issues = list(result['issues'])
            keys = set(result['data_events'])
            if keys & data_seen:
                issues.append('DATA physical identifiers overlap other input files')
            data_seen.update(keys)
            if issues:
                for name in ['events','candidates','losses']:
                    for row in result[name]:
                        row['file_audit_ok'] = False
            issue_count += len(issues)
            streams['issues'].write([dict(file_id=job['file_id'],sample=job['sample'],issue=x) for x in issues])
            streams['checks'].write([dict(file_id=job['file_id'],sample=job['sample'],raw_rows=result['raw_rows'],
                reconstructed_rows=result['derived_rows'],saved_before_mbc_rows=result['anchor_rows'],
                initial_payload_reproduced=result['reproduced'],issues=len(issues))])
            streams['manifest'].write([dict(file_id=job['file_id'],sample=job['sample'],stem=job['stem'],
                raw_path=job['raw_path'],anchor_path=job['anchor_path'],normalization=job['normalization'])])
            for name in ['events','candidates','losses']:
                streams[name].write(result[name])
            norm = job['normalization']
            for name, table in aggregate.items():
                for key, (events, candidates) in result[name].items():
                    table[(job['sample'],)+key] += [events, events*norm, events*norm**2, candidates]
            for location, fields in result['availability'].items():
                for name, present in fields.items():
                    coverage[job['sample'],location,name] += [int(present),int(not present)]
            for location, trees in result['schema'].items():
                for tree, branches in trees.items():
                    streams['schema'].write([dict(file_id=job['file_id'],sample=job['sample'],location=location,
                                                 tree=tree,branch=branch) for branch in branches])
            streams['nonfinite'].write([dict(file_id=job['file_id'],sample=job['sample'],variable=name,count=n)
                                        for name,n in result['nonfinite'].items() if n])
            if i % 100 == 0 or i == len(jobs):
                elapsed = time.monotonic()-start
                print(f'{i}/{len(jobs)} files; issues={issue_count}; elapsed={elapsed:.0f}s; '
                      f'rough remaining={elapsed/i*(len(jobs)-i)/3600:.2f}h',flush=True)
    finally:
        for stream in streams.values():
            stream.close()
    keys_by_table = {'counts':['sample','stage','category','view'],
                     'modes':['sample','stage','category','view','tag_mode'],
                     'focus':['sample','stage','view','mode_group'],
                     'signal':['sample','stage','category','tag_mode','sig_mode','mass_region']}
    filenames = {'counts':'initial_stage_summary.csv','modes':'tag_mode_presence.csv',
                 'focus':'focus_mode_groups.csv','signal':'tag_signal_mass_composition.csv'}
    tables = {}
    for name, table in aggregate.items():
        cols = keys_by_table[name]
        rows = [dict(zip(cols,k),events=int(v[0]),sumw=v[1],sumw2=v[2],candidates=int(v[3]),stage_order=STAGES.index(k[1]))
                for k,v in table.items()]
        tables[name] = pd.DataFrame(rows,columns=cols+['events','sumw','sumw2','candidates','stage_order'])
        tables[name].to_csv(output / filenames[name],index=False)
    counts = tables['counts']
    anchor_reference = pd.read_csv(args.reference / 'before_mbc_reference.csv')
    reference_ok = True
    for row in anchor_reference.to_dict('records'):
        observed = counts[counts['sample'].eq(row['sample']) & counts.stage.eq('before_Mbc_cut') &
                          counts.category.eq(row['category']) & counts['view'].eq(row['view'])]
        if len(observed) != 1 or not np.allclose(observed[['events','sumw','sumw2','candidates']].iloc[0],
                                                [row[k] for k in ['events','sumw','sumw2','candidates']]):
            reference_ok = False
    if not reference_ok:
        issue_count += 1
        with (output / 'aggregate_audit.txt').open('w',encoding='utf-8') as handle:
            handle.write('Reconstructed before_Mbc totals differ from the previously supplied full trace.\n')
    ratios = []
    for (stage,cat,view), group in counts.groupby(['stage','category','view']):
        data = group.loc[group['sample'].eq('DATA'),'events'].sum()
        mc = group.loc[group['sample'].ne('DATA'),'sumw'].sum()
        ratios.append(dict(stage=stage,stage_order=STAGES.index(stage),category=cat,view=view,
                           data=int(data),mc_normalization_only=mc,ratio=data/mc if mc else np.nan))
    ratios = pd.DataFrame(ratios).sort_values(['category','view','stage_order'])
    ratios.to_csv(output / 'initial_data_mc_ratios.csv',index=False)
    pd.DataFrame([dict(sample=k[0],location=k[1],variable=k[2],files_available=int(v[0]),files_missing=int(v[1]))
                  for k,v in coverage.items()]).to_csv(output / 'branch_coverage.csv',index=False)
    metadata = dict(source_files=len(jobs),workers=args.workers,audit_issues=issue_count,
        before_mbc_matches_previous_full_trace=reference_ok,stages=STAGES,
        normalization='Original reference normalization only; no detector/FEI weights.',
        raw_roots=dict(mc=str(args.raw_mc_base),data=str(args.raw_data_base)),
        reference_sha256={name:hashlib.sha256((args.reference/name).read_bytes()).hexdigest()
                          for name in ['files.csv','before_mbc_reference.csv']},
        code_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                     for name in ['trace_fei_initial.py','trace_fei_stages.py','fei_stage_details.py']},
        versions=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,uproot=uproot.__version__),
        caveats=['Raw input is after upstream FEI/skim reconstruction, not generator-level data.',
                 'No BCS is rerun. in_saved_before_mbc refers to the pre-BCS anchor, not a BCS winner.',
                 'Category/view/mode-presence event counts overlap; do not add them.',
                 'Unknown truth is preserved; no physical cause or significance is automatically assigned.'])
    (output/'run_info.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    lines = ['Initial FEI selection trace',f'Files: {len(jobs)}',f'Audit issues: {issue_count}',
             f'Before-Mbc totals match previous full trace: {reference_ok}',
             'Read initial_data_mc_ratios.csv for categories 2/3, view=control_K_Kpi.',
             'Read focus_mode_groups.csv for modes30_32, both all_modes and control_K_Kpi.',
             'First target loss: initial_event_losses.csv.gz; candidate losses: initial_candidate_trace.csv.gz.',
             'All yields use source normalization only. This does not identify a physical cause automatically.']
    (output/'summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (output/'RUNNING').unlink()
    marker = 'REVIEW_REQUIRED' if issue_count else 'COMPLETE'
    (output/marker).write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(f'{marker}: {output}',flush=True)
    return 2 if issue_count else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-mc-base',type=Path,default=Path(RAW_MC))
    parser.add_argument('--raw-data-base',type=Path,default=Path(RAW_DATA))
    parser.add_argument('--mc-base',type=Path,default=Path(base.MC_BASE))
    parser.add_argument('--data-base',type=Path,default=Path(base.DATA_BASE))
    parser.add_argument('--reference',type=Path,default=HERE/'reference')
    parser.add_argument('--version',default='validation_v009')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--output',type=Path,default=Path('fei_initial_trace_results'))
    parser.add_argument('--preflight-only',action='store_true')
    try:
        return run(parser.parse_args())
    except Exception as error:
        print('ERROR: '+str(error),file=sys.stderr,flush=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
