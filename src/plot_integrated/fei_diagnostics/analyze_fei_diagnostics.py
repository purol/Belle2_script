"""Explain small FEI yield ratios using the instrumented original calculation.

Run: python3 analyze_fei_diagnostics.py fei_diagnostics_output
Errors are counting/MC-statistical diagnostics, not a new calibration fit.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import math

import numpy as np
import pandas as pd
from scipy.stats import chi2

SAMPLES = ['DATA', 'SIGNAL', 'CHG', 'MIX', 'UUBAR', 'DDBAR', 'SSBAR', 'CHARM']
FACTORS = ['normalization', 'input_fei', 'pi0', 'kaon_id', 'pion_id', 'fake_pid', 'ks0',
           'knn', 'xnn', 'photon_multiplicity', 'kpklkl', 'ksklkl', 'kstarklkl', 'xklkl', 'b_to_dkl']
LABELS = {0: 'Neutral tags (combined)', 1: 'Charged category 1', 2: 'Charged category 2', 3: 'Charged category 3'}
NAN = float('nan')


def moments(weights):
    w = np.asarray(weights, dtype=float)
    return float(w.sum()), float(np.square(w).sum())


def ratio_stats(data_count, weights):
    """Garwood data limits plus MC sumw2, combined in quadrature for plotting.

    This is not a likelihood significance or an exact interval for a ratio.
    The symmetric delta-method error matches the original TH1 division for n>0.
    """
    n = int(data_count)
    total, sumw2 = moments(weights)
    alpha = 1 - 0.682689492137
    lower = 0.0 if n == 0 else 0.5 * chi2.ppf(alpha/2, 2*n)
    upper = 0.5 * chi2.ppf(1-alpha/2, 2*(n+1))
    result = dict(data=n, mc=total, mc_sumw2=sumw2, mc_stat_error=math.sqrt(sumw2),
                  mc_raw_count=len(weights), mc_neff=total**2/sumw2 if sumw2 else NAN,
                  data_low=lower, data_high=upper, ratio=NAN, ratio_error_delta=NAN,
                  ratio_error_low=NAN, ratio_error_high=NAN, ratio_status='nonpositive_mc')
    if total > 0:
        mc_term=n*math.sqrt(sumw2)/total**2
        result.update(ratio=n/total, ratio_status='ok',
                      ratio_error_delta=math.sqrt(n/total**2+mc_term**2) if n > 0 else NAN,
                      ratio_error_low=math.hypot((n-lower)/total, mc_term),
                      ratio_error_high=math.hypot((upper-n)/total, mc_term))
    return result


def truth_labels(values):
    values=np.asarray(values, dtype=float)
    return np.where(values == 1, 'correct', np.where(values == 0, 'incorrect', 'unknown'))


def weight_without(frame, factor):
    # Multiplication, not division: handles zero correction factors correctly.
    return frame[['f_'+x for x in FACTORS if x != factor]].prod(axis=1).to_numpy()


def validate(events, reference, files, cuts):
    expected={'file_id','entry','sample','category','upsilon_id','tag_mode','sig_mode','selection_mass',
              'weight','official_fei','tag_is_signal','sig_is_signal'} | {'f_'+x for x in FACTORS}
    missing=expected-set(events.columns)
    if missing: raise ValueError(f'Missing columns: {sorted(missing)}')
    if events.empty: raise ValueError('No selected rows; check input paths and cuts.')
    if set(files['sample']) != set(SAMPLES): raise ValueError('Expected input files for DATA and all 7 MC components.')
    if files['file_id'].duplicated().any(): raise ValueError('Duplicate file IDs in manifest.')
    if files['path'].duplicated().any(): raise ValueError('The same input file is listed more than once.')
    if not set(events['sample']).issubset(SAMPLES): raise ValueError('Unknown sample label.')
    if events.duplicated(['file_id','entry']).any(): raise ValueError('A selected input row was exported twice.')
    info=files.set_index('file_id')
    for fid,group in events.groupby('file_id'):
        if fid not in info.index: raise ValueError('Unknown input file ID.')
        if not (group['sample'] == info.loc[fid,'sample']).all(): raise ValueError('Sample/file mismatch.')
        entry=group['entry'].to_numpy()
        if not (np.isfinite(entry).all() and (entry == np.round(entry)).all() and
                (entry >= 0).all() and (entry < info.loc[fid,'entries']).all()):
            raise ValueError('Invalid input entry index.')
    numeric=['weight','official_fei','upsilon_id','tag_mode','sig_mode','selection_mass']+['f_'+x for x in FACTORS]
    if not np.isfinite(events[numeric].to_numpy()).all(): raise ValueError('Non-finite mandatory value.')
    if not events['upsilon_id'].isin([0,1]).all(): raise ValueError('Invalid tag charge ID.')
    for name in ['tag_mode','sig_mode']:
        if not (events[name] == np.round(events[name])).all(): raise ValueError(f'Noninteger {name}.')
    if not (events['tag_mode'] >= 0).all(): raise ValueError('Negative tag mode ID.')
    expected_cat=np.where(events['upsilon_id'] == 1, 0,
                         np.where(events['tag_mode'].isin([0,1,3,4]),1,
                                  np.where(events['tag_mode'].isin([15,16,18,19]),2,3)))
    if not np.array_equal(expected_cat,events['category']): raise ValueError('Category mapping mismatch.')
    if not ((events['sig_mode'] > -0.5)&(events['sig_mode'] < 2.5)&(events['selection_mass'] < 1.0)).all():
        raise ValueError('Export contains a row failing the original local selection.')
    product=events[['f_'+x for x in FACTORS]].prod(axis=1)
    if not np.allclose(product,events['weight'],rtol=1e-10,atol=1e-12): raise ValueError('Correction product mismatch.')
    data=events['sample'] == 'DATA'
    if not np.all(events.loc[data,'weight'] == 1): raise ValueError('Data must have unit weight.')
    bb=events['sample'].isin(['SIGNAL','CHG','MIX'])
    if not np.all(events.loc[bb,'f_input_fei'] == 1): raise ValueError('Unexpected pre-existing BB FEI correction.')
    if not (events['official_fei'] > 0).all(): raise ValueError('Invalid official factor.')
    if sorted(reference['category'].tolist()) != [0,1,2,3]: raise ValueError('Expected 4 original reference rows.')
    checks=[]
    for row in reference.to_dict('records'):
        selected=events[events['category'] == row['category']]
        mc=selected[selected['sample'] != 'DATA']
        observed={'data':int((selected['sample'] == 'DATA').sum()),
                  'mc':float(mc['weight'].sum()),'mc_sumw2':float(np.square(mc['weight']).sum())}
        for key,value in observed.items():
            ok=math.isclose(value,float(row[key]),rel_tol=1e-9,abs_tol=1e-10)
            checks.append(dict(category=row['category'],quantity=key,original=row[key],exported=value,passed=ok))
            if not ok: raise ValueError(f'Original calculation not reproduced: category {row["category"]} {key}.')
    final=cuts[cuts['stage'] == 3].groupby(['sample','category']).agg(n=('n','sum'),sumw=('sumw','sum'))
    selected=events.groupby(['sample','category']).agg(n=('weight','size'),sumw=('f_normalization','sum'))
    for key in set(final.index) | set(selected.index):
        a=final.loc[key] if key in final.index else pd.Series({'n':0,'sumw':0})
        b=selected.loc[key] if key in selected.index else pd.Series({'n':0,'sumw':0})
        if a['n'] != b['n'] or not math.isclose(a['sumw'],b['sumw'],rel_tol=1e-9,abs_tol=1e-10):
            raise ValueError(f'Cutflow and selected rows disagree: {key}')
    return pd.DataFrame(checks)


def category_table(events):
    rows=[]
    for cat in range(4):
        selected=events[events['category'] == cat]
        data=int((selected['sample'] == 'DATA').sum())
        mc=selected[selected['sample'] != 'DATA']
        for scenario,weights in [('original',mc['weight']),('normalization_only',mc['f_normalization']),
                                 ('official_fei_applied',mc['weight']*mc['official_fei'])]:
            rows.append(dict(category=cat,scenario=scenario,**ratio_stats(data,weights)))
    return pd.DataFrame(rows)


def decomposition(events, columns):
    mc=events[events['sample'] != 'DATA'].copy()
    mc['tag_truth']=truth_labels(mc['tag_is_signal'])
    mc['sig_truth']=truth_labels(mc['sig_is_signal'])
    rows=[]
    for key,group in mc.groupby(columns,dropna=False):
        key=(key,) if len(columns)==1 else key
        total,square=moments(group['weight'])
        denominator=mc.loc[mc['category'] == group['category'].iloc[0],'weight'].sum()
        rows.append(dict(zip(columns,key),n=len(group),sumw=total,sumw2=square,
                         mc_stat_error=math.sqrt(square),fraction_of_category=total/denominator if denominator>0 else NAN))
    return pd.DataFrame(rows)


def mode_table(events, variable):
    rows=[]
    for (cat,mode),group in events.groupby(['category',variable]):
        data=int((group['sample'] == 'DATA').sum())
        mc=group[group['sample'] != 'DATA']
        stats=ratio_stats(data,mc['weight'])
        rows.append(dict(category=cat,**{variable:int(mode)},mc_minus_data=stats['mc']-data,**stats))
    return pd.DataFrame(rows)


def weight_tables(events):
    sensitivity=[]; distributions=[]
    for cat,group in events.groupby('category'):
        data=int((group['sample'] == 'DATA').sum())
        mc=group[group['sample'] != 'DATA']
        if mc.empty: continue
        nominal=float(mc['weight'].sum())
        for factor in FACTORS[2:]:
            weights=weight_without(mc,factor)
            stats=ratio_stats(data,weights)
            sensitivity.append(dict(category=cat,removed_factor=factor,nominal_mc=nominal,
                                    change_in_mc=stats['mc']-nominal,**stats))
        for sample,sample_group in mc.groupby('sample'):
            for factor in FACTORS:
                values=sample_group['f_'+factor]
                distributions.append(dict(category=cat,sample=sample,factor=factor,n=len(values),
                                          minimum=float(values.min()),median=float(values.median()),maximum=float(values.max()),
                                          zeros=int((values==0).sum()),negative=int((values<0).sum())))
    return pd.DataFrame(sensitivity),pd.DataFrame(distributions)


def cutflow_table(cuts):
    rows=[]
    for (cat,stage),group in cuts.groupby(['category','stage']):
        data=int(group.loc[group['sample']=='DATA','n'].sum())
        mc=group[group['sample']!='DATA']
        total=float(mc['sumw'].sum()); square=float(mc['sumw2'].sum())
        rows.append(dict(category=cat,stage=stage,data=data,mc_normalization_only=total,
                         mc_raw_count=int(mc['n'].sum()),mc_sumw2=square,
                         ratio=data/total if total>0 else NAN,
                         stage_description=['input_ntuple','signal_mode','signal_mode_and_mass','original_final'][int(stage)]))
    return pd.DataFrame(rows)


def scan_table(events):
    rows=[]
    settings={'mva_bb':[0.0,0.2,0.4,0.6,0.8,0.86,0.9,0.95],
              'tag_probability':[0.001,0.003,0.01,0.03,0.1,0.3]}
    for cat,group in events.groupby('category'):
        for variable,thresholds in settings.items():
            finite=np.isfinite(group[variable])
            if not finite.all():
                rows.append(dict(category=cat,variable=variable,threshold=NAN,status='unavailable_partial_coverage',
                                 missing_rows=int((~finite).sum())))
                continue
            for threshold in thresholds:
                selected=group[group[variable]>threshold]
                n=int((selected['sample']=='DATA').sum())
                stats=ratio_stats(n,selected.loc[selected['sample']!='DATA','weight'])
                rows.append(dict(category=cat,variable=variable,threshold=threshold,status='tightening_within_input',
                                 missing_rows=0,**stats))
    return pd.DataFrame(rows)


def identity_table(events):
    rows=[]
    keys=['experiment','run','event','production']
    for sample,group in events.groupby('sample'):
        if not np.isfinite(group[keys].to_numpy()).all():
            rows.append(dict(sample=sample,status='identifiers_incomplete',selected_rows=len(group),repeated_event_rows=NAN))
        else:
            unique=len(group[keys].drop_duplicates())
            rows.append(dict(sample=sample,status='identifiers_available',selected_rows=len(group),
                             repeated_event_rows=len(group)-unique))
    return pd.DataFrame(rows)


def plots(events,summary,destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.grid':False})
    nominal=summary[summary['scenario']=='original'].sort_values('category')
    official=summary[summary['scenario']=='official_fei_applied'].sort_values('category')
    fig,ax=plt.subplots(figsize=(9,4.5))
    for table,offset,label in [(nominal,-0.07,'Original yield ratio'),(official,0.07,'After official FEI factors')]:
        valid=table['ratio_status']=='ok'
        v=table[valid]
        ax.errorbar(v['category']+offset,v['ratio'],yerr=v[['ratio_error_low','ratio_error_high']].to_numpy().T,
                    fmt='o',capsize=3,label=label)
    ax.axhline(1,color='gray',linestyle='--')
    ax.set_xticks(range(4),['Neutral','Charged 1','Charged 2','Charged 3'])
    ax.set_ylabel('Data / MC'); ax.legend(); fig.tight_layout()
    fig.savefig(destination/'category_ratios.png',dpi=160); plt.close(fig)
    ranges={'mva_bb':np.linspace(0,1,16),'tag_probability':np.linspace(-4,0,17),
            'tag_mbc':np.linspace(5.20,5.30,21),'tag_delta_e':np.linspace(-0.3,0.3,21),
            'selection_mass':np.linspace(0,1,16)}
    coverage=[]
    for cat,group in events.groupby('category'):
        for variable,bins in ranges.items():
            if ((group['sample']!='DATA')&(group['weight']<0)).any():
                coverage.append(dict(category=cat,variable=variable,status='skipped_signed_weights'))
                continue
            x=group[variable].to_numpy(float)
            if not np.isfinite(x).all() or (variable=='tag_probability' and not (x>0).all()):
                coverage.append(dict(category=cat,variable=variable,status='skipped_incomplete_or_invalid'))
                continue
            if variable=='tag_probability': x=np.log10(x)
            sample=group['sample'].to_numpy(); w=group['weight'].to_numpy()
            data=np.histogram(x[sample=='DATA'],bins=bins)[0]
            mc=np.histogram(x[sample!='DATA'],bins=bins,weights=w[sample!='DATA'])[0]
            mc2=np.histogram(x[sample!='DATA'],bins=bins,weights=w[sample!='DATA']**2)[0]
            centers=(bins[1:]+bins[:-1])/2; widths=np.diff(bins)
            fig,(ax,ratio)=plt.subplots(2,1,figsize=(7,6),sharex=True,gridspec_kw={'height_ratios':[3,1]})
            base=np.zeros_like(mc)
            for component in SAMPLES[1:]:
                mask=sample==component
                h=np.histogram(x[mask],bins=bins,weights=w[mask])[0]
                ax.bar(centers,h,width=widths,bottom=base,alpha=0.65,label=component)
                base+=h
            alpha=1-0.682689492137
            lo=np.where(data>0,0.5*chi2.ppf(alpha/2,2*data),0)
            hi=0.5*chi2.ppf(1-alpha/2,2*(data+1))
            ax.errorbar(centers,data,yerr=[data-lo,hi-data],fmt='k.',label='Data')
            ax.stairs(mc+np.sqrt(mc2),bins,color='gray',linewidth=0.6)
            ax.stairs(mc-np.sqrt(mc2),bins,color='gray',linewidth=0.6)
            ok=mc>0
            ratio.errorbar(centers[ok],data[ok]/mc[ok],yerr=[(data-lo)[ok]/mc[ok],(hi-data)[ok]/mc[ok]],fmt='k.')
            relative=np.sqrt(mc2[ok])/mc[ok]
            ratio.bar(centers[ok],2*relative,width=widths[ok],bottom=1-relative,color='gray',alpha=.3)
            ratio.axhline(1,color='gray',linestyle='--'); ratio.set_ylabel('D/MC')
            ratio.set_xlabel('log10(tag probability)' if variable=='tag_probability' else variable)
            ax.set_title(LABELS[int(cat)]); ax.set_ylabel('Weighted entries'); ax.legend(ncol=4,fontsize=8)
            fig.tight_layout(); fig.savefig(destination/f'category_{cat}_{variable}.png',dpi=140); plt.close(fig)
            coverage.append(dict(category=cat,variable=variable,status='plotted',
                                 outside_range=int(((x<bins[0])|(x>bins[-1])).sum())))
    pd.DataFrame(coverage).to_csv(destination/'plot_coverage.csv',index=False)


def analyze(directory: Path, make_plots=True):
    if not (directory/'COMPLETE').is_file(): raise ValueError('No COMPLETE marker: the ROOT export did not finish.')
    events=pd.read_csv(directory/'selected_events.csv')
    reference=pd.read_csv(directory/'reference.csv')
    files=pd.read_csv(directory/'files.csv')
    cuts=pd.read_csv(directory/'input_cutflow.csv')
    coverage=pd.read_csv(directory/'branch_availability.csv')
    checks=validate(events,reference,files,cuts)
    out=directory/'analysis'; out.mkdir(exist_ok=True)
    checks.to_csv(out/'reproduction_check.csv',index=False)
    summary=category_table(events)
    modes=mode_table(events,'tag_mode')
    sensitivity,factor_ranges=weight_tables(events)
    identity=identity_table(events)
    tables={
        'category_summary':summary,
        'tag_mode_ratios':modes,
        'signal_mode_ratios':mode_table(events,'sig_mode'),
        'mc_composition':decomposition(events,['category','sample']),
        'tag_truth':decomposition(events,['category','tag_truth']),
        'truth_by_sample_and_mode':decomposition(events,['category','sample','tag_mode','tag_truth','sig_truth']),
        'weight_sensitivity':sensitivity,
        'weight_ranges':factor_ranges,
        'cutflow_ratios':cutflow_table(cuts),
        'threshold_scans':scan_table(events),
        'event_identifiers':identity,
    }
    for name,table in tables.items(): table.to_csv(out/(name+'.csv'),index=False)
    if make_plots: plots(events,summary,out)
    lines=['FEI calibration diagnostics','', 'Original selected yields and MC sumw2: REPRODUCED.',
           'These outputs localize discrepancies; they do not establish a physical cause.',
           'Ratios are selected data / total MC, not background-subtracted tag efficiencies.',
           'Errors assume independent selected entries. No systematic uncertainty is included.',
           'Asymmetric plot bars combine Garwood data intervals and MC sumw2 approximately.',
           'The two FEI scenarios share events and are correlated; do not compare them as independent measurements.',
           'Correction-removal scenarios are sensitivity checks, not assigned uncertainties.','']
    for cat in range(4):
        row=summary[(summary['category']==cat)&(summary['scenario']=='original')].iloc[0]
        off=summary[(summary['category']==cat)&(summary['scenario']=='official_fei_applied')].iloc[0]
        lines += [LABELS[cat],f"  Data={row['data']}, MC={row['mc']:.6g} +/- {row['mc_stat_error']:.4g} (MC stat)",
                  f"  Original ratio={row['ratio']:.6g}; delta-method error={row['ratio_error_delta']:.4g}",
                  f"  Ratio after applying official FEI factors={off['ratio']:.6g}"]
        truth=tables['tag_truth']
        for rec in truth[truth['category']==cat].to_dict('records'):
            lines.append(f"  MC tag truth {rec['tag_truth']}: {rec['sumw']:.6g} ({100*rec['fraction_of_category']:.1f}% of weighted MC)")
        largest=modes[modes['category']==cat].sort_values('mc_minus_data',ascending=False).head(3)
        for rec in largest.to_dict('records'):
            lines.append(f"  Tag mode {rec['tag_mode']}: data={rec['data']}, MC={rec['mc']:.6g}, MC-data={rec['mc_minus_data']:.6g}")
        local=sensitivity[sensitivity['category']==cat] if not sensitivity.empty else sensitivity
        if not local.empty:
            rec=local.loc[local['change_in_mc'].abs().idxmax()]
            lines.append(f"  Largest one-factor removal effect: {rec['removed_factor']}, delta MC={rec['change_in_mc']:.6g}")
        lines.append('')
    unavailable=coverage[coverage['available']==0]
    lines+=['Checks and limitations',f'  Missing optional branch/file combinations: {len(unavailable)} (see branch_availability.csv).',
            '  Missing/nonbinary truth values are unknown, never assigned to incorrect tags.',
            '  Cutflow stages 0/1/2/3: input ntuple / local signal mode / local mass / final tag option.',
            '  Cutflow weights contain normalization only. Earlier skim/MVA selections cannot be undone.',
            '  Threshold scans only tighten cuts within the selected input; partial-coverage variables are skipped.',
            '  MC truth partitions have no observable data equivalents; no data/correct-tag-only ratio is computed.',
            '  The Jpsi tree may have different candidate ordering and is not joined by row number.']
    negative=int(((events['sample']!='DATA')&(events['weight']<0)).sum())
    if negative: lines.append(f'  WARNING: {negative} negative MC weights; stacked plots/fractions need signed-weight interpretation.')
    for rec in identity.to_dict('records'):
        if rec['status']=='identifiers_incomplete':
            lines.append(f"  {rec['sample']}: event multiplicity could not be checked (missing identifiers).")
        elif rec['repeated_event_rows']:
            lines.append(f"  {rec['sample']}: {int(rec['repeated_event_rows'])} repeated event-ID rows. Review multiple candidates/dataset provenance; row-based errors may not represent independent events.")
    if files['sync_keys'].isna().any() or (files['sync_keys'].fillna('')=='').any():
        lines.append('  Some input files lack shared tree identifiers; only their entry counts could be checked for alignment.')
    lines+=['','Read category_summary.csv first; then tag_mode_ratios.csv, truth_by_sample_and_mode.csv,',
            'weight_sensitivity.csv and threshold_scans.csv to identify which comparison changes the deficit.']
    (out/'summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('\n'.join(lines[:8])); print(f'Full diagnostic report: {out / "summary.txt"}')
    return tables


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args()
    analyze(args.directory,make_plots=not args.no_plots)
