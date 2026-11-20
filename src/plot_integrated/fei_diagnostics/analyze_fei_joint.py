"""Joint diagnostics from an existing FEI export; no ROOT rerun is needed.

Usage: python analyze_fei_joint.py EXPORT_DIRECTORY OUTPUT_DIRECTORY
The output is exploratory and does not fit a new tag efficiency correction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import analyze_fei_diagnostics as diag


def joint_table(events):
    frame = events.copy()
    frame['probability_region'] = np.where(
        ~np.isfinite(frame.tag_probability), 'unknown',
        np.where(frame.tag_probability > .01, 'gt_0.01', 'le_0.01'))
    rows = []
    for (cat, mode, region), group in frame.groupby(
            ['category', 'sig_mode', 'probability_region']):
        mc = group[group['sample'] != 'DATA']
        n = int((group['sample'] == 'DATA').sum())
        stats = diag.ratio_stats(n, mc.weight)
        labels = diag.truth_labels(mc.tag_is_signal)
        truth = {f'mc_{kind}': float(mc.loc[labels == kind, 'weight'].sum())
                 for kind in ['correct', 'incorrect', 'unknown']}
        rows.append(dict(category=int(cat), sig_mode=int(mode),
                         probability_region=region,
                         mc_minus_data=stats['mc'] - n, **truth, **stats))
    return pd.DataFrame(rows)


def application_table(events):
    """Apply the thesis's stated prescription to this same control sample.

    This is a correlated, algebraic illustration, NOT an independent closure
    test or a fit. Unknown truth stays in the non-correct component.
    """
    rows = []
    for cat, group in events.groupby('category'):
        mc = group[group['sample'] != 'DATA']
        n = int((group['sample'] == 'DATA').sum())
        total = float(mc.weight.sum())
        if total <= 0:
            raise ValueError('Application diagnostic requires positive MC yield.')
        ratio = n / total
        correct = mc.tag_is_signal.eq(1)
        correct_yield = float(mc.loc[correct, 'weight'].sum())
        official_rest = float((mc.loc[~correct, 'weight'] *
                               mc.loc[~correct, 'official_fei']).sum())
        hybrid = ratio * correct_yield + official_rest
        rows.append(dict(category=int(cat), data=n, mc_original=total,
                         inclusive_ratio=ratio, correct_mc_original=correct_yield,
                         official_mc_noncorrect_including_unknown=official_rest,
                         mc_after_ratio_on_all=ratio * total,
                         mc_after_ratio_on_correct_official_on_rest=hybrid,
                         data_over_hybrid=n / hybrid if hybrid > 0 else np.nan))
    return pd.DataFrame(rows)


def plot_joint(table, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    selected = table[table.category == 3].set_index(
        ['sig_mode', 'probability_region'])
    keys = [(mode, region) for mode in range(3)
            for region in ['le_0.01', 'gt_0.01']]
    if not all(key in selected.index for key in keys):
        return  # Do not fabricate unavailable cells.
    selected = selected.loc[keys]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(figsize=(10.2, 5.5), layout='constrained')
    bottom = np.zeros(len(keys))
    for name, label, color in [('correct', 'MC: correct tag', '#4c91bb'),
                               ('incorrect', 'MC: incorrect tag', '#e49b50'),
                               ('unknown', 'MC: unknown truth', '#b5b5b5')]:
        value = selected['mc_' + name].to_numpy()
        ax.bar(x, value, bottom=bottom, width=.64, label=label, color=color)
        bottom += value
    ax.errorbar(x, selected.mc, yerr=selected.mc_stat_error,
                fmt='none', color='#555555', capsize=3, label='MC statistical error')
    ax.errorbar(x, selected.data,
                yerr=np.array([selected.data - selected.data_low,
                               selected.data_high - selected.data]),
                fmt='o', color='black', capsize=3, label='Data (68% Garwood)')
    ax.set_xticks(x, [f'Mode {mode}\npFEI {"<= " if region == "le_0.01" else "> "}0.01'
                     for mode, region in keys])
    ax.set_ylabel('Selected yield (weighted MC)')
    ax.set_title('Charged category 3: signal-side mode and tag probability')
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9)
    ax.text(.99, .61, 'Exploratory split; no systematic errors\nMC has no FEI calibration applied',
            transform=ax.transAxes, ha='right', va='top', fontsize=9)
    fig.savefig(output / 'category3_joint.png', dpi=160)
    fig.savefig(output / 'category3_joint.pdf')
    plt.close(fig)


def analyze(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output:
        raise ValueError('Use a separate output directory to preserve the export.')
    names = ['selected_events', 'files', 'reference', 'input_cutflow', 'branch_availability']
    inputs = {name: pd.read_csv(source / (name + '.csv')) for name in names}
    e, f = inputs['selected_events'], inputs['files']
    checks = diag.validate(e, inputs['reference'], f, inputs['input_cutflow'])
    joint = joint_table(e)
    original = diag.category_table(e).query('scenario == "original"')
    # Disjoint partitions must reproduce all three moments, not just ratios.
    sums = joint.groupby('category')[['data', 'mc', 'mc_sumw2']].sum()
    for row in original.itertuples():
        for quantity in ['data', 'mc', 'mc_sumw2']:
            if not np.isclose(sums.loc[row.category, quantity], getattr(row, quantity),
                              rtol=1e-10, atol=1e-10):
                raise ValueError('Joint partition does not conserve ' + quantity)
    if not np.allclose(joint[['mc_correct', 'mc_incorrect', 'mc_unknown']].sum(axis=1),
                       joint.mc, rtol=1e-10, atol=1e-10):
        raise ValueError('Truth partition does not conserve MC yield.')
    application = application_table(e)
    identifiers = []
    for sample, g in e.groupby('sample'):
        keys = ['experiment', 'run', 'event']
        complete = bool(np.isfinite(g[keys].to_numpy()).all())
        identifiers.append(dict(sample=sample, rows=len(g), keys=';'.join(keys),
                                keys_complete=complete,
                                repeated_rows=int(g.duplicated(keys, keep=False).sum()) if complete else None,
                                production_complete=bool(g.production.notna().all())))
    selected_availability = inputs['branch_availability'].loc[
        inputs['branch_availability'].file_id.isin(e.file_id)]
    availability = selected_availability.groupby('variable').available.agg(
        available_files='sum', selected_files='size').reset_index()
    hotspot = e[(e.category == 3) & (e.sig_mode == 0) & (e.tag_probability > .01)]
    hot_modes = diag.mode_table(hotspot, 'tag_mode')
    hot_events = hotspot.merge(f[['file_id', 'path']], on='file_id', validate='many_to_one')
    output.mkdir(parents=True, exist_ok=True)
    tables = {'reproduction_check': checks, 'joint_signal_mode_probability': joint,
              'application_diagnostic': application,
              'identifier_check': pd.DataFrame(identifiers),
              'selected_branch_availability': availability,
              'category3_mode0_high_probability_tag_modes': hot_modes,
              'category3_mode0_high_probability_events': hot_events,
              'category3_data_events': e[(e.category == 3) & (e['sample'] == 'DATA')]
                                      .merge(f[['file_id', 'path']], on='file_id', validate='many_to_one')}
    for name, table in tables.items():
        table.to_csv(output / (name + '.csv'), index=False)
    metadata = dict(source=str(source), selected_rows=len(e), input_files=len(f),
                    reproduced_checks=int(checks.passed.sum()),
                    partition_checks='passed',
                    complete_marker_present=(source / 'COMPLETE').exists(),
                    run_info_present=(source / 'run_info.txt').exists(),
                    inputs_sha256={name: hashlib.sha256((source / (name + '.csv')).read_bytes()).hexdigest()
                                   for name in names},
                    caveats=[
                        'The raw ROOT files were not rerun or independently inspected here.',
                        'Missing COMPLETE/run_info limits run provenance; copied numerical inputs were validated.',
                        'Ratios use total selected MC, with no control-background subtraction.',
                        'Exploratory mode/probability splits do not establish a physical cause or significance.',
                        'MC truth categories have no observable data equivalents.',
                        'Application table reuses the calibration data; no independent closure claim.',
                        'Repeated identifier checks omit production when unavailable.',
                        'Earlier skim selections and their efficiency are not recoverable from this export.'])
    (output / 'validation.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    plot_joint(joint, output)
    return tables


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('export_directory', type=Path)
    parser.add_argument('output_directory', type=Path)
    args = parser.parse_args()
    result = analyze(args.export_directory, args.output_directory)
    print('Validated original yields and joint partitions; wrote', args.output_directory)
    print(result['joint_signal_mode_probability'].query('category == 3')[
        ['sig_mode', 'probability_region', 'data', 'mc', 'mc_correct']].to_string(index=False))
