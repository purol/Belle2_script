"""Detailed category-2/3 histories from snapshots already loaded by the tracer.

No additional ROOT opens, input mutations, new BCS choices, or fitted weights.
"""
from collections import defaultdict
import numpy as np
import pandas as pd

from trace_fei_stages import EVENT, IDENTITY, StreamCSV, require, view_mask

VIEWS = ['control_K_Kpi', 'mode0']
COMPONENT_KEYS = ['sample', 'stage', 'stage_order', 'category', 'view',
                  'other_category', 'best_target_modes', 'truth_class', 'high_probability']


def modes(values):
    return ';'.join(str(int(x)) for x in sorted(set(values)))


def describe(group, category, view, sample):
    target_mask = group.category.eq(category) & view_mask(group, view)
    target = group[target_mask]
    if target.empty:
        return None
    other = group[group.category != category]
    outside = group[~target_mask]
    correct = int(target.tag_is_signal.eq(1).sum())
    incorrect = int(target.tag_is_signal.eq(0).sum())
    unknown = len(target) - correct - incorrect
    truth = ('data' if sample == 'DATA' else 'any_correct' if correct else
             'unknown_without_correct' if unknown else 'all_incorrect')
    best_probability = float(target.probability.max())
    best = target[target.probability == best_probability]
    return dict(category=category, view=view, n_all=len(group), n_target=len(target),
                n_other_category=len(other), other_category=bool(len(other)),
                n_outside_target=len(outside), target_modes=modes(target.tag_mode),
                best_target_modes=modes(best.tag_mode), target_maxima=len(best),
                n_correct_target=correct, n_incorrect_target=incorrect, n_unknown_target=unknown,
                truth_class=truth, max_target_probability=best_probability,
                max_other_category_probability=float(other.probability.max()) if len(other) else np.nan,
                max_outside_target_probability=float(outside.probability.max()) if len(outside) else np.nan,
                max_event_probability=float(group.probability.max()),
                target_at_event_maximum=bool(best_probability == group.probability.max()),
                high_probability=best_probability > .01,
                mode22_present=bool(target.tag_mode.eq(22).any()))


def divide(a, b):
    return float(a / b) if b > 0 else np.nan


class DetailExport:
    def __init__(self, output, stage_order):
        self.output = output
        self.stage_order = stage_order
        self.streams = {name: StreamCSV(output / filename) for name, filename in {
            'events': 'selection_event_trace.csv.gz',
            'candidates': 'selection_candidate_trace.csv.gz',
            'transitions': 'selection_transitions.csv.gz',
            'schema': 'representative_branch_inventory.csv'}.items()}
        self.components = defaultdict(lambda: np.zeros(4))
        self.availability = defaultdict(lambda: np.zeros(2, dtype=int))

    def process(self, job, frames, audit_ok):
        issues = []
        for stage in self.stage_order:
            frame = frames[stage]
            if frame.duplicated(IDENTITY).any():
                issues.append((stage, 'duplicate_candidate_identity_in_detailed_snapshot'))
            if stage in job['paths']:
                for field, available in frame.attrs.get('extra_availability', {}).items():
                    self.availability[job['sample'], stage, field] += [int(available), int(not available)]
                for tree, branches in frame.attrs.get('branch_schema', {}).items():
                    self.streams['schema'].write([
                        dict(file_id=job['file_id'], sample=job['sample'], stage=stage,
                             path=job['paths'][stage], tree=tree, branch=branch) for branch in branches])
        initial = frames[self.stage_order[0]]
        targets = initial[initial.category.isin([2, 3]) & view_mask(initial, 'control_K_Kpi')]
        origins = set(targets[EVENT].itertuples(index=False, name=None))
        if not origins:
            return issues
        prefix = dict(file_id=job['file_id'], sample=job['sample'],
                      normalization=job['normalization'], file_audit_ok=bool(audit_ok and not issues))
        final_ids = set(frames['final_output'][IDENTITY].itertuples(index=False, name=None))
        states = {}
        grouped = {}
        empty = initial.iloc[:0]
        for order, stage in enumerate(self.stage_order):
            groups = {key: g for key, g in frames[stage].groupby(EVENT, sort=False) if key in origins}
            grouped[stage] = groups
            for key, group in groups.items():
                identity = dict(zip(EVENT, (int(x) for x in key)))
                candidate_rows = [dict(prefix, stage=stage, stage_order=order, **row,
                                       matches_final_identity=tuple(row[k] for k in IDENTITY) in final_ids)
                                  for row in group.to_dict('records')]
                self.streams['candidates'].write(candidate_rows)
                for category in [2, 3]:
                    for view in VIEWS:
                        record = describe(group, category, view, job['sample'])
                        if record is None:
                            continue
                        row = dict(prefix, stage=stage, stage_order=order, **identity, **record)
                        self.streams['events'].write([row])
                        states[stage, key, category, view] = record
                        norm = job['normalization']
                        component = tuple(row[k] for k in COMPONENT_KEYS)
                        self.components[component] += [1, norm, norm**2, record['n_target']]
        for left, right in zip(self.stage_order[:-1], self.stage_order[1:]):
            for key, before in grouped[left].items():
                after = grouped[right].get(key, empty)
                for category in [2, 3]:
                    for view in VIEWS:
                        origin = states.get((left, key, category, view))
                        if origin is None:
                            continue
                        following = states.get((right, key, category, view))
                        row = dict(prefix, stage_before=left, stage_after=right,
                                   **dict(zip(EVENT, (int(x) for x in key))), category=category, view=view,
                                   before_n_target=origin['n_target'],
                                   after_n_target=following['n_target'] if following else 0,
                                   retained_target=following is not None, retained_event=bool(len(after)),
                                   before_other_category=origin['other_category'],
                                   before_best_target_modes=origin['best_target_modes'],
                                   before_truth_class=origin['truth_class'],
                                   before_max_target_probability=origin['max_target_probability'],
                                   after_categories=modes(after.category),
                                   after_other_category=bool(after.category.ne(category).any()))
                        self.streams['transitions'].write([row])
        return issues

    def close(self):
        for stream in self.streams.values():
            stream.close()

    def finish(self, original_summary):
        rows = [dict(zip(COMPONENT_KEYS, key), events=int(value[0]), sumw=value[1],
                     sumw2=value[2], candidates=int(value[3]))
                for key, value in self.components.items()]
        components = pd.DataFrame(rows, columns=COMPONENT_KEYS + ['events', 'sumw', 'sumw2', 'candidates'])
        components.to_csv(self.output / 'selection_components.csv', index=False)
        subset = original_summary[original_summary.category.isin([2, 3]) & original_summary['view'].isin(VIEWS)]
        for row in subset.to_dict('records'):
            group = components[(components['sample'] == row['sample']) & (components.stage == row['stage']) &
                               (components.category == row['category']) & (components['view'] == row['view'])]
            for field in ['events', 'sumw', 'sumw2', 'candidates']:
                require(np.isclose(group[field].sum(), row[field]),
                        'Detailed component totals disagree with full cutflow: ' + str((row['sample'], row['stage'], row['category'], row['view'], field)))
        availability = [dict(sample=k[0], stage=k[1], variable=k[2],
                             files_available=int(v[0]), files_missing=int(v[1]))
                        for k, v in self.availability.items()]
        pd.DataFrame(availability).to_csv(self.output / 'extra_branch_coverage.csv', index=False)
        cutflow = []
        for (stage, order, cat, view), group in subset.groupby(['stage', 'stage_order', 'category', 'view']):
            data = group[group['sample'] == 'DATA']
            mc = group[group['sample'] != 'DATA']
            d, m = data.events.sum(), mc.sumw.sum()
            cutflow.append(dict(stage=stage, stage_order=order, category=cat, view=view,
                                data=int(d), mc_events=int(mc.events.sum()), mc_sumw=m,
                                mc_sumw2=mc.sumw2.sum(), data_mc=divide(d, m)))
        cutflow = pd.DataFrame(cutflow).sort_values(['category', 'view', 'stage_order'])
        cutflow.to_csv(self.output / 'category23_cutflow.csv', index=False)
        # Closed gzip files can be streamed; do not reload all candidate histories.
        transition_totals = defaultdict(lambda: np.zeros(4))
        path = self.output / 'selection_transitions.csv.gz'
        try:
            for chunk in pd.read_csv(path, chunksize=100000):
                chunk['kept_weight'] = chunk.normalization * chunk.retained_target
                for key, group in chunk.groupby(['sample', 'stage_before', 'stage_after', 'category', 'view']):
                    transition_totals[key] += [len(group), group.retained_target.sum(),
                                               group.normalization.sum(), group.kept_weight.sum()]
        except pd.errors.EmptyDataError:
            pass
        # Paired-event transition totals must close to BOTH adjacent snapshots.
        for row in subset[subset.stage != self.stage_order[-1]].to_dict('records'):
            left = row['stage']
            right = self.stage_order[self.stage_order.index(left) + 1]
            numbers = transition_totals[row['sample'], left, right, row['category'], row['view']]
            following = subset[(subset['sample'] == row['sample']) & (subset.stage == right) &
                               (subset.category == row['category']) & (subset['view'] == row['view'])].iloc[0]
            require(np.allclose(numbers, [row['events'], following.events, row['sumw'], following.sumw]),
                    'Event lineage does not close between snapshots: ' + str((left, right, row['sample'], row['category'], row['view'])))
        effects = []
        for (category, view), group in cutflow.groupby(['category', 'view']):
            ordered = group.sort_values('stage_order').to_dict('records')
            for left, right in zip(ordered[:-1], ordered[1:]):
                dr = divide(right['data'], left['data'])
                mr = divide(right['mc_sumw'], left['mc_sumw'])
                effects.append(dict(category=category, view=view, stage_before=left['stage'],
                                    stage_after=right['stage'], data_before=left['data'], data_after=right['data'],
                                    mc_before=left['mc_sumw'], mc_after=right['mc_sumw'],
                                    data_retention=dr, mc_retention=mr,
                                    relative_retention=divide(dr, mr),
                                    data_mc_before=left['data_mc'], data_mc_after=right['data_mc']))
        effects = pd.DataFrame(effects)
        effects.to_csv(self.output / 'category23_cut_retention.csv', index=False)
        lines = ['Category 2/3 descriptive diagnostics. All yields use source normalization only.',
                 'No physical cause or statistical significance is automatically assigned.',
                 'Largest observed ratio drops below are exploratory and use overlapping cohorts.',
                 'Pre-BCS target truth means any target candidate; it is not necessarily the winner.',
                 'Unknown truth and unavailable optional branches are not converted to zero.', '']
        for (category, view), group in effects.groupby(['category', 'view']):
            finite = group[np.isfinite(group.relative_retention)]
            if finite.empty:
                continue
            low = finite.loc[finite.relative_retention.idxmin()]
            label = 'largest decrease' if low.relative_retention < 1 else 'no observed decrease'
            lines.append(f'Category {category}, {view}: {label}; {low.stage_before} -> {low.stage_after}; '
                         f'Data/MC {low.data_mc_before:.5g} -> {low.data_mc_after:.5g}; '
                         f'relative retention {low.relative_retention:.5g}.')
        lines += ['', 'If the first saved stage already has a deficit, investigate upstream production/reconstruction.',
                  'Branch inventory is representative, not proof that production configurations match.',
                  'Category-2/3 components and paired transitions close to the original full cutflow.']
        (self.output / 'category23_diagnostic_summary.txt').write_text('\n'.join(lines)+'\n', encoding='utf-8')
