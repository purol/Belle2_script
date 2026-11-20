"""Read existing J/psi intermediate ROOT files and audit cuts and BCS.

No reconstruction, random reselection, or changes to the physics files.
Default reference/ contains the user's 2026-11-20 FEI export manifest.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import hashlib
import json
import math
import multiprocessing
from pathlib import Path, PurePosixPath
import sys
import time

import numpy as np
import pandas as pd
import uproot

HERE = Path(__file__).resolve().parent
MC_BASE = '/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_MC_Jpsi'
DATA_BASE = '/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_data_Jpsi'
SAMPLES = ['DATA', 'SIGNAL', 'CHG', 'MIX', 'UUBAR', 'DDBAR', 'SSBAR', 'CHARM']
SNAPSHOTS = ['before_Mbc_cut', 'before_delE_cut', 'before_Eecl_cut',
             'before_missing_momentum_theta_cut', 'before_psig_cut', 'before_Dveto_cut', 'final_output']
EVENT = ['experiment', 'run', 'event']
GROUP = EVENT + ['ncandidates']
IDENTITY = GROUP + ['candidate']
ID_BRANCHES = {name: '__' + name + '__' for name in IDENTITY}
FIELDS = {
    'Upsilon': {'upsilon_id': 'extraInfo__bodecayModeID__bc',
                'eecl': 'extraInfo__boEeclv200__bc',
                'missing_theta': 'missingMomentumOfEvent_theta'},
    'Bsig': {'sig_mode': 'Bsig_daughter_0_extraInfo_decayModeID',
             'mass': 'Bsig_M', 'sig_p': 'Bsig_useCMSFrame_p'},
    'Btag': {'tag_mode': 'Btag_extraInfo_decayModeID',
             'probability': 'Btag_extraInfo_SignalProbability',
             'mbc': 'Btag_Mbc', 'delta_e': 'Btag_deltaE'},
}
TRUTH = {'tag_is_signal': 'Btag_isSignal', 'tag_mc_errors': 'Btag_mcErrors',
         'tag_mc_pdg': 'Btag_mcPDG'}
# Optional observables present in the supplied Loader.h SaveFile implementation.
# Absence is recorded, never interpreted as a measured zero.
EXTRA_FIELDS = {'Upsilon': {
    'roe_tracks_clean': 'nROE_Tracks__bocleanMask__bc',
    'roe_tracks_loose': 'nROE_Tracks__bolooseMask__bc',
    'roe_neutral_ecl': 'nROE_NeutralECLClusters__bocleanMask__bc',
    'roe_pi0': 'nROE_ParticlesInList__bopi0__clmyneutralPion__bc',
    'roe_ks': 'nROE_ParticlesInList__boK_S0__clmyKaonshort__bc',
    'remaining_tracks': 'nRemainingTracksInEvent',
    'missing_p': 'missingMomentumOfEvent',
    'missing_energy_cms': 'missingEnergyOfEventCMS',
    'eecl_matched': 'extraInfo__boEeclv200_matched__bc',
    'eecl_unmatched': 'extraInfo__boEeclv200_unmatched__bc',
    'n_jpsi': 'nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc',
}}
PAYLOAD = IDENTITY + [name for fields in FIELDS.values() for name in fields]
VIEWS = ['all_modes', 'control_K_Kpi', 'mode0']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_stem(name, stage, version='validation_v009'):
    """Strip only the exact filename suffixes used by Analysis_main.cc."""
    require(name.endswith('.root'), 'Not a ROOT filename: ' + name)
    stem = name[:-5]
    if stage == 'reference':
        require(stem.startswith('MVAoutput_'), 'Unexpected MVA filename: ' + name)
        stem = stem[len('MVAoutput_'):]
        suffix = '_' + version + '_final_output_data_after_MVA'
    elif stage == 'final_output':
        suffix = '_' + version + '_final_output'
    else:
        require(stage in SNAPSHOTS, 'Unknown stage ' + stage)
        suffix = '_' + stage
    require(stem.endswith(suffix), 'Unexpected suffix for ' + stage + ': ' + name)
    return stem[:-len(suffix)]


def read_root(path, details=False, capture_schema=False):
    columns = {}
    availability = {}
    extra_availability, schema = {}, {}
    tag = None
    # These inputs are local/mounted files. A file handle avoids fsspec's
    # unnecessary asyncio/socket setup, including on restricted Windows hosts.
    with Path(path).open('rb') as handle, uproot.open(handle, object_cache=None, array_cache=None) as root:
        n = None
        for tree_name, fields in FIELDS.items():
            require(tree_name in root, f'{path}: missing {tree_name}')
            tree = root[tree_name]
            branch_names = set(tree.keys())
            if capture_schema:
                schema[tree_name] = sorted(branch_names)
            if tree_name == 'Btag':
                tag = tree
            if n is None:
                n = tree.num_entries
            require(tree.num_entries == n, f'{path}: different tree entry counts')
            mapping = dict(ID_BRANCHES, **fields)
            missing = set(mapping.values()) - branch_names
            require(not missing, f'{path}: {tree_name} missing branches {sorted(missing)}')
            arrays = (tree.arrays(list(mapping.values()), library='np', how=dict)
                      if n else {branch: np.array([], dtype=float) for branch in mapping.values()})
            for name, branch in mapping.items():
                values = np.asarray(arrays[branch])
                require(values.ndim == 1 and values.dtype.kind in 'biuf',
                        f'{path}: not a numeric scalar branch {branch}')
                require(np.isfinite(values).all(), f'{path}: nonfinite mandatory {branch}')
                if name in ID_BRANCHES:
                    require((values == np.round(values)).all(), f'{path}: noninteger {branch}')
                    values = values.astype(np.int64)
                if name in columns:
                    require(np.array_equal(columns[name], values),
                            f'{path}: {tree_name} event/candidate alignment differs for {branch}')
                else:
                    columns[name] = values
            if details:
                extras = EXTRA_FIELDS.get(tree_name, {})
                present = [branch for branch in extras.values() if branch in branch_names]
                extra_arrays = (tree.arrays(present, library='np', how=dict)
                                if present and n else {})
                for name, branch in extras.items():
                    extra_availability[name] = branch in branch_names
                    values = np.asarray(extra_arrays[branch]) if branch in extra_arrays else np.full(n, np.nan)
                    require(values.ndim == 1 and values.dtype.kind in 'biuf',
                            f'{path}: non-scalar optional branch {branch}')
                    columns[name] = values
        # Reuse the already deserialized TTree. With object_cache=None, looking
        # it up again would reread all its metadata for every snapshot.
        for name, branch in TRUTH.items():
            availability[name] = branch in tag
            columns[name] = tag[branch].array(library='np') if branch in tag else np.full(n, np.nan)
    frame = pd.DataFrame(columns)
    for field in ['upsilon_id', 'sig_mode', 'tag_mode']:
        require((frame[field] == np.round(frame[field])).all(), f'{path}: noninteger {field}')
        frame[field] = frame[field].astype(np.int64)
    require(frame.upsilon_id.isin([0, 1]).all(), f'{path}: invalid Upsilon ID')
    require((frame.tag_mode >= 0).all(), f'{path}: invalid tag mode')
    require(frame.probability.between(0, 1).all(), f'{path}: probability outside [0,1]')
    frame['category'] = np.where(frame.upsilon_id == 1, 0,
                                np.where(frame.tag_mode.isin([0, 1, 3, 4]), 1,
                                         np.where(frame.tag_mode.isin([15, 16, 18, 19]), 2, 3)))
    frame['entry'] = np.arange(len(frame))
    if details:
        frame.attrs['extra_availability'] = extra_availability
        frame.attrs['branch_schema'] = schema
    return frame, availability


def load_job(job):
    """Worker only reads files; all output writing stays in the main process."""
    return {stage: read_root(path, details=job.get('details', False),
                             capture_schema=job.get('capture_schema', False))
            for stage, path in job['paths'].items()}


def loaded_jobs(jobs, workers=1):
    """Bound memory and queued reads; preserve exactly the serial file order."""
    require(workers >= 1, '--workers must be a positive integer')
    if workers == 1:
        for job in jobs:
            yield job, load_job(job)
        return
    # Spawn avoids inheriting Uproot handles or partially initialized threads.
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as pool:
        iterator = iter(jobs)
        pending = deque()
        for _ in range(2 * workers):
            job = next(iterator, None)
            if job is None:
                break
            pending.append((job, pool.submit(load_job, job)))
        while pending:
            job, future = pending.popleft()
            result = future.result()  # Exceptions propagate; never skip a failed file.
            next_job = next(iterator, None)
            if next_job is not None:
                pending.append((next_job, pool.submit(load_job, next_job)))
            yield job, result


def apply_next_cut(frame, stage):
    masks = {
        'before_Mbc_cut': frame.mbc > 5.27,
        'before_delE_cut': (frame.delta_e > -.2) & (frame.delta_e < .2),
        'before_Eecl_cut': frame.eecl < 1.3,
        'before_missing_momentum_theta_cut': (frame.missing_theta > .297) & (frame.missing_theta < 2.618),
        'before_psig_cut': (frame.sig_p > .5) & (frame.sig_p < 2.96),
        # Strict interior is removed, so exactly 1.84 and 1.89 are retained.
        'before_Dveto_cut': ~((frame.mass > 1.84) & (frame.mass < 1.89)),
    }
    return frame.loc[masks[stage]].copy()


def max_probability_candidates(frame):
    """Follow Loader's CONTIGUOUS (exp,run,event,ncandidates) blocks exactly.

    Do not prefilter tag categories or signal modes. Do not introduce a tolerance
    for ties: Loader compares double values with exact equality.
    """
    if frame.empty:
        return frame.copy()
    block = frame[GROUP].ne(frame[GROUP].shift()).any(axis=1).cumsum()
    maximum = frame.groupby(block, sort=False).probability.transform('max')
    return frame.loc[frame.probability == maximum].copy()


def view_mask(frame, view):
    if view == 'all_modes':
        return np.ones(len(frame), dtype=bool)
    if view == 'mode0':
        return (frame.sig_mode == 0) & (frame.mass < 1.)
    return frame.sig_mode.isin([0, 1, 2]) & (frame.mass < 1.)


def signatures(frame, columns=PAYLOAD):
    if frame.empty:
        return Counter()
    return Counter(frame[columns].itertuples(index=False, name=None))


def audit_bcs(before, final):
    after_veto = apply_next_cut(before, 'before_Dveto_cut')
    highest = max_probability_candidates(after_veto)
    issues = []
    final_sigs, veto_sigs, max_sigs = signatures(final), signatures(after_veto), signatures(highest)
    if final_sigs - veto_sigs:
        issues.append('final_candidate_not_in_after_Dveto')
    if final_sigs - max_sigs:
        issues.append('final_candidate_not_a_maximum_of_its_Loader_block')
    all_events = set(after_veto[EVENT].itertuples(index=False, name=None))
    final_events = set(final[EVENT].itertuples(index=False, name=None))
    if all_events != final_events:
        issues.append('physical_event_set_changed_across_BCS')
    if final.duplicated(EVENT).any():
        issues.append('multiple_final_rows_for_one_physical_event')
    if after_veto.duplicated(IDENTITY).any():
        issues.append('duplicate_candidate_identifiers_before_BCS')
    if not after_veto.empty:
        variants = after_veto.groupby(EVENT).ncandidates.nunique()
        if (variants > 1).any():
            issues.append('multiple_ncandidates_values_for_one_physical_event')
        block = after_veto[GROUP].ne(after_veto[GROUP].shift()).any(axis=1).cumsum()
        starts = after_veto.loc[~block.duplicated(), GROUP]
        if starts.duplicated(GROUP).any():
            issues.append('same_Loader_group_in_noncontiguous_blocks')
    return after_veto, highest, issues


def event_trace(before, after_veto, highest, final):
    """Overlapping origin cohorts: event has a target candidate before D veto."""
    if before.empty:
        return []
    grouped = [{key: g for key, g in f.groupby(EVENT, sort=False)}
               for f in [after_veto, highest, final]]
    empty = before.iloc[:0]
    rows = []
    for key, original in before.groupby(EVENT, sort=False):
        veto, maxima, winners = [d.get(key, empty) for d in grouped]
        for cat in range(4):
            for view in ['control_K_Kpi', 'mode0']:
                def target(g):
                    return g[(g.category == cat) & view_mask(g, view)]
                initial = target(original)
                if initial.empty:
                    continue
                v, h, w = target(veto), target(maxima), target(winners)
                if len(winners) != 1 or signatures(winners) - signatures(maxima):
                    outcome = 'audit_issue_no_unique_valid_winner' if not veto.empty else 'no_event_after_Dveto'
                elif not w.empty:
                    outcome = 'kept_target'
                elif v.empty:
                    outcome = 'target_removed_by_Dveto'
                elif h.empty:
                    outcome = 'lost_at_max_probability'
                else:
                    outcome = 'lost_among_equal_maxima'
                winner = winners.iloc[0] if len(winners) == 1 else None
                rows.append(dict(zip(EVENT, (int(v) for v in key)), origin_category=cat, view=view,
                                 origin_candidates=len(initial), after_veto_target_candidates=len(v),
                                 maximal_target_candidates=len(h), all_maximal_candidates=len(maxima),
                                 original_event_candidates=len(original), final_event_rows=len(winners),
                                 outcome=outcome,
                                 winner_category=int(winner.category) if winner is not None else -1,
                                 winner_sig_mode=int(winner.sig_mode) if winner is not None else -1,
                                 winner_tag_mode=int(winner.tag_mode) if winner is not None else -1,
                                 winner_probability=float(winner.probability) if winner is not None else np.nan,
                                 max_target_probability=float(initial.probability.max()),
                                 max_event_probability=float(veto.probability.max()) if len(veto) else np.nan))
    return rows


def preflight(reference, mc_base, data_base, version, scope):
    manifest = pd.read_csv(reference / 'files.csv')
    events = pd.read_csv(reference / 'selected_events.csv')
    require(set(manifest['sample']) == set(SAMPLES), 'Reference must contain DATA and all seven MC samples.')
    require(not manifest.file_id.duplicated().any(), 'Duplicate file IDs in reference.')
    require(not manifest.path.duplicated().any(), 'Duplicate paths in reference.')
    require(np.isfinite(manifest.normalization).all() and (manifest.normalization > 0).all(),
            'Reference has invalid normalization.')
    require((manifest.loc[manifest['sample'] == 'DATA', 'normalization'] == 1).all(), 'Data weight is not one.')
    stages = SNAPSHOTS if scope == 'full' else ['before_Dveto_cut', 'final_output']
    jobs = []
    for sample in SAMPLES:
        refs = manifest[manifest['sample'] == sample].copy()
        refs['stem'] = [source_stem(PurePosixPath(path).name, 'reference', version) for path in refs.path]
        require(not refs.stem.duplicated().any(), 'Duplicate source stems for ' + sample)
        expected = set(refs.stem)
        base = (data_base / 'SIGNAL_analysis' if sample == 'DATA' else mc_base / (sample + '_analysis')) / version
        maps = {}
        for stage in stages:
            directory = base / stage
            require(directory.is_dir(), 'Missing directory: ' + str(directory))
            pairs = [(source_stem(p.name, stage, version), p) for p in sorted(directory.glob('*.root'))]
            require(len(pairs) == len({k for k, p in pairs}), 'Duplicate file stems: ' + str(directory))
            maps[stage] = dict(pairs)  # Deliberately NOT recursive; no Merge subdirectories.
            missing, extra = expected - set(maps[stage]), set(maps[stage]) - expected
            require(not missing and not extra,
                    f'{directory}: cohort differs from reference; missing={len(missing)} {sorted(missing)[:3]}, '
                    f'extra={len(extra)} {sorted(extra)[:3]}. Do not silently drop unmatched files.')
        for row in refs.to_dict('records'):
            jobs.append(dict(row, paths={stage: str(maps[stage][row['stem']]) for stage in stages}))
    require(set(events.file_id).issubset(set(manifest.file_id)), 'Unknown reference event file ID.')
    require(not events.duplicated(['file_id', 'entry']).any(), 'Duplicate reference rows.')
    return jobs, events, stages


class StreamCSV:
    def __init__(self, path):
        self.handle = gzip.open(path, 'wt', encoding='utf-8', newline='') if str(path).endswith('.gz') else path.open('w', encoding='utf-8', newline='')
        self.writer = None

    def write(self, rows):
        for row in rows:
            if self.writer is None:
                self.writer = csv.DictWriter(self.handle, fieldnames=list(row))
                self.writer.writeheader()
            self.writer.writerow(row)

    def close(self):
        self.handle.close()


def run(args):
    output = args.output.resolve()
    workers = getattr(args, 'workers', 1)
    details_enabled = getattr(args, 'details', False)
    require(not details_enabled or args.scope == 'full', '--details requires --scope full')
    require(workers >= 1, '--workers must be a positive integer')
    require(not output.exists(), 'Output exists; choose a new directory: ' + str(output))
    jobs, ref_events, stages = preflight(args.reference, args.mc_base, args.data_base, args.version, args.scope)
    print(f'Preflight passed: {len(jobs)} source files, {len(stages)} stored snapshots each; workers={workers}.', flush=True)
    if args.preflight_only:
        return 0
    schema_seen = set()
    for job in jobs:
        job['details'] = details_enabled
        job['capture_schema'] = details_enabled and job['sample'] not in schema_seen
        schema_seen.add(job['sample'])
    output.mkdir(parents=True)
    (output / 'RUNNING').write_text('Incomplete until COMPLETE; REVIEW_REQUIRED means audit issues.\n')
    streams = {name: StreamCSV(output / filename) for name, filename in {
        'trace': 'bcs_event_trace.csv.gz', 'candidates': 'bcs_candidate_trace.csv.gz',
        'checks': 'file_checks.csv', 'issues': 'audit_issues.csv', 'manifest': 'source_files.csv',
        'coverage': 'branch_availability.csv'}.items()}
    totals = defaultdict(lambda: np.zeros(4))  # events, event sumw, event sumw2, candidates
    outcomes = defaultdict(lambda: np.zeros(3))  # events, sumw, sumw2
    data_seen = defaultdict(set)
    refs = {fid: g for fid, g in ref_events.groupby('file_id')}
    issue_count = 0
    start = time.monotonic()
    stage_order = (SNAPSHOTS[:6] if args.scope == 'full' else ['before_Dveto_cut']) + [
        'after_Dveto_reconstructed', 'after_max_probability_reconstructed', 'final_output']
    detail_export = None
    if details_enabled:
        from fei_stage_details import DetailExport
        detail_export = DetailExport(output, stage_order)
    for sample in SAMPLES:
        for stage in stage_order:
            for cat in range(4):
                for view in VIEWS:
                    totals[sample, stage, cat, view]  # Ensure explicit zero cells, including empty samples.

    def add_issue(job, stage, message):
        nonlocal issue_count
        issue_count += 1
        streams['issues'].write([dict(file_id=job['file_id'], sample=job['sample'], stem=job['stem'], stage=stage, issue=message)])

    def count_stage(job, stage, frame):
        if frame.empty:
            return
        if job['sample'] == 'DATA':
            keys = set(frame[EVENT].itertuples(index=False, name=None))
            if keys & data_seen[stage]:
                add_issue(job, stage, 'DATA physical identifiers overlap other source files; totals need deduplication')
            data_seen[stage].update(keys)
        for cat in range(4):
            for view in VIEWS:
                subset = frame[(frame.category == cat) & view_mask(frame, view)]
                n = len(subset.drop_duplicates(EVENT))
                norm = job['normalization']
                totals[job['sample'], stage, cat, view] += [n, n * norm, n * norm**2, len(subset)]

    try:
        for job_index, (job, loaded) in enumerate(loaded_jobs(jobs, workers)):
            file_issue_start = issue_count
            frames = {}
            local_issues = []
            for stage in stages:
                frame, available = loaded[stage]
                frames[stage] = frame
                count_stage(job, stage, frame)
                streams['coverage'].write([dict(file_id=job['file_id'], stage=stage, variable=k, available=v)
                                           for k, v in available.items()])
                streams['manifest'].write([dict(file_id=job['file_id'], sample=job['sample'], stem=job['stem'],
                                               stage=stage, path=job['paths'][stage], rows=len(frame), normalization=job['normalization'])])
            if args.scope == 'full':
                for left, right in zip(SNAPSHOTS[:5], SNAPSHOTS[1:6]):
                    if signatures(apply_next_cut(frames[left], left)) != signatures(frames[right]):
                        local_issues.append((left + ' -> ' + right, 'stored_snapshot_differs_from_supplied_cut_code'))
            before, final = frames['before_Dveto_cut'], frames['final_output']
            after_veto, highest, bcs_issues = audit_bcs(before, final)
            local_issues.extend(('BCS', text) for text in bcs_issues)
            count_stage(job, 'after_Dveto_reconstructed', after_veto)
            count_stage(job, 'after_max_probability_reconstructed', highest)
            # The pre-MVA BCS output must reproduce the same selected candidates
            # as the user's previously exported post-MVA (score application) file.
            selected = final[view_mask(final, 'control_K_Kpi')]
            ref = refs.get(job['file_id'], ref_events.iloc[:0])
            keys = EVENT + ['candidate', 'category', 'tag_mode', 'sig_mode']
            final_reproduced = signatures(selected, keys) == signatures(ref, keys)
            if not final_reproduced:
                local_issues.append(('final_output -> FEI export', 'selected_candidate_identity_or_count_mismatch'))
            for stage, text in local_issues:
                add_issue(job, stage, text)
            if detail_export is not None:
                detailed_frames = dict(frames, after_Dveto_reconstructed=after_veto,
                                       after_max_probability_reconstructed=highest)
                for stage, message in detail_export.process(job, detailed_frames,
                                                           issue_count == file_issue_start):
                    add_issue(job, stage, message)
            streams['checks'].write([dict(file_id=job['file_id'], sample=job['sample'],
                                          final_candidates=len(final), reference_selected=len(ref),
                                          selected_candidate_reproduction=final_reproduced,
                                          local_audit_issue_count=issue_count - file_issue_start)])
            traces = event_trace(before, after_veto, highest, final)
            for row in traces:
                record = dict(file_id=job['file_id'], sample=job['sample'], source_stem=job['stem'],
                              normalization=job['normalization'], file_audit_ok=issue_count == file_issue_start, **row)
                streams['trace'].write([record])
                key = (job['sample'], row['origin_category'], row['view'], row['outcome'], row['winner_category'], row['winner_sig_mode'])
                norm = job['normalization']
                outcomes[key] += [1, norm, norm**2]
            # Save all competitors for physical events with a charged control candidate.
            targets = before[(before.category > 0) & view_mask(before, 'control_K_Kpi')]
            event_keys = set(targets[EVENT].itertuples(index=False, name=None))
            veto_entries, max_entries = set(after_veto.entry), set(highest.entry)
            final_ids = set(final[IDENTITY].itertuples(index=False, name=None))
            candidate_records = []
            for row in before.to_dict('records'):
                if tuple(row[k] for k in EVENT) in event_keys:
                    candidate_records.append(dict(file_id=job['file_id'], sample=job['sample'],
                                                  normalization=job['normalization'], **row,
                                                  passes_Dveto=row['entry'] in veto_entries,
                                                  is_maximal=row['entry'] in max_entries,
                                                  matches_final_identity=tuple(row[k] for k in IDENTITY) in final_ids))
            streams['candidates'].write(candidate_records)
            if (job_index + 1) % 100 == 0 or job_index + 1 == len(jobs):
                elapsed = time.monotonic() - start
                eta = elapsed / (job_index + 1) * (len(jobs) - job_index - 1)
                print(f'{job_index + 1}/{len(jobs)} files; issues={issue_count}; '
                      f'elapsed={elapsed:.0f}s; rough_remaining={eta/3600:.2f}h', flush=True)
    finally:
        for stream in streams.values():
            stream.close()
        if detail_export is not None:
            detail_export.close()
    rows = [dict(sample=k[0], stage=k[1], stage_order=stage_order.index(k[1]), category=k[2], view=k[3],
                 events=int(v[0]), sumw=v[1], sumw2=v[2], candidates=int(v[3])) for k, v in totals.items()]
    summary = pd.DataFrame(rows).sort_values(['stage_order', 'category', 'view', 'sample'])
    summary.to_csv(output / 'stage_summary.csv', index=False)
    if detail_export is not None:
        detail_export.finish(summary)
    ratios = []
    for (stage, order, cat, view), group in summary.groupby(['stage', 'stage_order', 'category', 'view']):
        d = int(group.loc[group['sample'] == 'DATA', 'events'].sum())
        mc = group[group['sample'] != 'DATA']
        total, square = mc.sumw.sum(), mc.sumw2.sum()
        ratios.append(dict(stage=stage, stage_order=order, category=cat, view=view, data=d,
                           mc_normalization_only=total, mc_sumw2=square,
                           ratio=d / total if total > 0 else np.nan,
                           ratio_error_delta=math.sqrt(d/total**2 + d*d*square/total**4) if total > 0 and d > 0 else np.nan))
    pd.DataFrame(ratios).sort_values(['stage_order', 'category', 'view']).to_csv(output / 'data_mc_stage_ratios.csv', index=False)
    pd.DataFrame([dict(sample=k[0], origin_category=k[1], view=k[2], outcome=k[3],
                       winner_category=k[4], winner_sig_mode=k[5], events=int(v[0]), sumw=v[1], sumw2=v[2])
                  for k, v in outcomes.items()]).to_csv(output / 'bcs_migration_summary.csv', index=False)
    metadata = dict(scope=args.scope, workers=workers, details=details_enabled,
                    source_files=len(jobs), audit_issues=issue_count,
                    normalization='Only per-source normalization from previous FEI export; no detector or FEI weights.',
                    reference=str(args.reference.resolve()),
                    reference_sha256={name: hashlib.sha256((args.reference / name).read_bytes()).hexdigest()
                                      for name in ['files.csv', 'selected_events.csv']},
                    source_sha256=json.loads((HERE / 'source_hashes.json').read_text()) if (HERE / 'source_hashes.json').exists() else {},
                    versions=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__, uproot=uproot.__version__),
                    diagnostic_code_sha256={name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                                            for name in ['trace_fei_stages.py', 'fei_stage_details.py']
                                            if (HERE / name).exists()},
                    stage_order=stage_order,
                    caveats=[
                        'Before BCS, an event can appear in more than one category/view. Do not add category counts.',
                        'MC physical event namespace includes source file; identical IDs in distinct productions are not merged.',
                        'DATA identifier overlaps across files are flagged, not silently counted as validated independent events.',
                        'Reconstructed maxima use exact contiguous Loader grouping and include all competing modes.',
                        'The actual saved winner is used; random seeds/hash behavior are not regenerated.',
                        'final_output includes IsBCSValid handling. Ambiguous/missing/nonmaximal winners require review.',
                        'Earliest saved snapshot follows the initial ROE/vertex and reconstruction/skim selections.',
                        'Stage ratios share events. Errors are statistical only; no p-values or physical cause are assigned.',
                        'If audit issues exist, all totals and migration attributions are provisional.'])
    (output / 'run_info.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    lines = ['FEI stage trace', f'Scope: {args.scope}', f'Files: {len(jobs)}', f'Audit issues: {issue_count}',
             'Read data_mc_stage_ratios.csv (view=mode0, category=3), then bcs_migration_summary.csv.',
             'All stage MC yields use normalization only; they are not the original corrected FEI factors.',
             'BCS origin cohorts overlap. Compare event survival/migration within each cohort.',
             'No random BCS was rerun. Actual final_output candidates were used.',
             'If issues > 0: read audit_issues.csv before interpreting physics.',
             'This scan cannot identify losses before before_Mbc_cut.']
    if detail_export is not None:
        lines += ['Detailed category 2/3 cutflow: category23_cutflow.csv and category23_cut_retention.csv.',
                  'Composition: selection_components.csv and selection_event_trace.csv.gz.',
                  'All competitors and available ROE observables: selection_candidate_trace.csv.gz.',
                  'No physical cause is automatically assigned; upstream production settings are not in these snapshots.']
    (output / 'summary.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    (output / 'RUNNING').unlink()
    marker = 'REVIEW_REQUIRED' if issue_count else 'COMPLETE'
    (output / marker).write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f'{marker}: {output}', flush=True)
    return 2 if issue_count else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mc-base', type=Path, default=Path(MC_BASE))
    parser.add_argument('--data-base', type=Path, default=Path(DATA_BASE))
    parser.add_argument('--version', default='validation_v009')
    parser.add_argument('--reference', type=Path, default=HERE / 'reference')
    parser.add_argument('--output', type=Path, default=Path('fei_stage_trace_results'))
    parser.add_argument('--scope', choices=['bcs', 'full'], default='full')
    parser.add_argument('--details', action='store_true',
                        help='Export category 2/3 event/candidate histories and conditional cut retention (full scope).')
    parser.add_argument('--workers', type=int, default=1,
                        help='Concurrent read processes (default 1; try 4 for large local/mounted datasets).')
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    try:
        return run(args)
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr, flush=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
