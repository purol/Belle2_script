"""Numerical/integrity regression tests; no ROOT installation required."""
from pathlib import Path
import hashlib
import re
import tempfile

import numpy as np
import pandas as pd
import pytest

import analyze_fei_diagnostics as diag
from build_diagnostic_macro import generate


def fixture(directory):
    directory.mkdir(parents=True,exist_ok=True)
    rows=[]; files=[]; coverage=[]; cuts=[]
    optional=['tag_is_signal','sig_is_signal','tag_mc_errors','tag_mc_pdg','tag_probability','tag_mbc',
              'tag_delta_e','mva_bb','xs_mass','tag_wrong_daughters','tag_missing_gamma_ok',
              'experiment','run','event','production','candidate']
    for fid,sample in enumerate(diag.SAMPLES):
        norm=1 if sample=='DATA' else 2
        files.append(dict(file_id=fid,sample=sample,path=f'/synthetic/{sample}.root',entries=4,
                          normalization=norm,sync_keys='__event__;'))
        for cat,mode in [(0,0),(1,3),(2,18),(3,31)]:
            row=dict(file_id=fid,entry=cat,sample=sample,category=cat,upsilon_id=int(cat==0),tag_mode=mode,
                     sig_mode=cat%3,selection_mass=0.5+0.1*(cat%3),weight=1.,official_fei=0.8 if sample in ['SIGNAL','CHG','MIX'] else 1.)
            row.update({name:float('nan') for name in optional})
            row.update(tag_is_signal=1 if sample=='SIGNAL' else (0 if sample in ['CHG','MIX'] else np.nan),
                       sig_is_signal=1 if sample=='SIGNAL' else 0,
                       tag_probability=.02+cat*.02,tag_mbc=5.279,tag_delta_e=.01,mva_bb=.65+cat*.1,
                       experiment=7,run=23,event=cat,production=fid,candidate=0)
            row.update({'f_'+name:1. for name in diag.FACTORS})
            row['f_normalization']=norm
            row['f_pi0']=1 if sample=='DATA' else .5
            rows.append(row)
            for stage in range(4):
                cuts.append(dict(sample=sample,category=cat,tag_mode=mode,stage=stage,n=1,sumw=norm,sumw2=norm**2))
        for variable in optional:
            coverage.append(dict(file_id=fid,sample=sample,variable=variable,tree='T',branch=variable,
                                 available=int(variable not in ['xs_mass','tag_wrong_daughters'])))
    events=pd.DataFrame(rows)
    reference=pd.DataFrame([dict(category=c,data=1.,mc=7.,mc_sumw2=7.) for c in range(4)])
    tables={'selected_events':events,'files':pd.DataFrame(files),'reference':reference,
            'input_cutflow':pd.DataFrame(cuts),'branch_availability':pd.DataFrame(coverage)}
    for name,table in tables.items(): table.to_csv(directory/(name+'.csv'),index=False)
    (directory/'COMPLETE').write_text('synthetic fixture only\n')
    (directory/'run_info.txt').write_text('schema_version=1\nSYNTHETIC TEST DATA\n')
    return events,reference,pd.DataFrame(files),pd.DataFrame(cuts)


def test_ratio_and_zero_count():
    result=diag.ratio_stats(7,[10,10,10])
    assert result['ratio']==pytest.approx(7/30)
    assert result['mc_sumw2']==300
    assert result['ratio_error_delta']==pytest.approx(np.sqrt(7/900+49*300/30**4))
    zero=diag.ratio_stats(0,[2,3])
    assert zero['ratio']==0
    assert zero['ratio_error_high']>0
    assert np.isnan(zero['ratio_error_delta'])
    assert diag.ratio_stats(2,[1,-1])['ratio_status']=='nonpositive_mc'


def test_full_report_reproduces_and_decomposes(tmp_path):
    events,reference,files,cuts=fixture(tmp_path)
    checks=diag.validate(events,reference,files,cuts)
    assert checks['passed'].all()
    results=diag.analyze(tmp_path,make_plots=False)
    summary=results['category_summary']
    assert np.allclose(summary.loc[summary.scenario=='original','ratio'],1/7)
    assert np.allclose(summary.loc[summary.scenario=='official_fei_applied','ratio'],1/6.4)
    truth=results['tag_truth']
    assert (truth.loc[truth.tag_truth=='unknown','sumw']==4).all()
    assert (truth.loc[truth.tag_truth=='incorrect','sumw']==2).all()
    removed=results['weight_sensitivity']
    assert (removed.loc[removed.removed_factor=='pi0','mc']==14).all()
    assert (tmp_path/'analysis/summary.txt').is_file()


def test_reproduction_rejects_weight_change(tmp_path):
    events,reference,files,cuts=fixture(tmp_path)
    events.loc[events['sample']=='CHG','weight']=2
    events.loc[events['sample']=='CHG','f_pi0']=1
    with pytest.raises(ValueError,match='not reproduced'): diag.validate(events,reference,files,cuts)


def test_reproduction_rejects_duplicated_row(tmp_path):
    events,reference,files,cuts=fixture(tmp_path)
    events=pd.concat([events,events.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match='exported twice'): diag.validate(events,reference,files,cuts)


def test_missing_truth_not_incorrect_and_scans_skip_partial_coverage(tmp_path):
    events,_,_,_=fixture(tmp_path)
    assert list(diag.truth_labels([1,0,np.nan,-1,2]))==['correct','incorrect','unknown','unknown','unknown']
    events.loc[(events['category']==2)&(events['sample']=='CHG'),'mva_bb']=np.nan
    scans=diag.scan_table(events)
    row=scans[(scans['category']==2)&(scans['variable']=='mva_bb')]
    assert len(row)==1
    assert row['status'].iloc[0]=='unavailable_partial_coverage'


def test_zero_factor_can_be_removed_without_division(tmp_path):
    events,_,_,_=fixture(tmp_path)
    mc=events[events['sample']=='CHG'].copy()
    mc['f_pi0']=0; mc['weight']=0
    assert np.all(diag.weight_without(mc,'pi0')==2)
    assert np.all(diag.weight_without(mc,'ks0')==0)


def test_incomplete_export_and_bad_cutflow_fail(tmp_path):
    events,reference,files,cuts=fixture(tmp_path)
    cuts.loc[(cuts['stage']==3)&(cuts['sample']=='DATA'),'n']=2
    with pytest.raises(ValueError,match='Cutflow'): diag.validate(events,reference,files,cuts)
    (tmp_path/'COMPLETE').unlink()
    with pytest.raises(ValueError,match='did not finish'): diag.analyze(tmp_path,make_plots=False)


def test_generated_macro_preserves_weight_expression_and_source(tmp_path):
    here=Path(__file__).resolve().parent
    source=(here/'source/THStack_Jpsi_FEI_efficiency.c').read_bytes()
    assert hashlib.sha256(source).hexdigest()=='b0a7b899aa0bb63289157b840500f01b8fe9affa9b3276857d639606cf40768b'
    generated=generate(source)
    assert generated==(here/'THStack_Jpsi_FEI_diagnostics.C').read_text(encoding='utf-8')
    expression=lambda s:re.findall(r'double weights = ([^;]+);',s)
    assert expression(generated)==expression(source.decode('utf-8'))
    assert generated.count('feidiag::output().input(')==2
    assert generated.count('feidiag::output().selected(')==2
    assert 'fopen("NEWFEIcal.txt"' not in generated
    selection='if ((Bsig_ID > -0.5) && (Bsig_ID < 2.5) && (Bsig_M < 1.0)) {}'
    assert generated.count(selection)==source.decode('utf-8').count(selection)==2
    assert 'FEI_DIAG_MC_BASE' in generated and 'FEI_DIAG_DATA_BASE' in generated


if __name__=='__main__':
    import sys
    if len(sys.argv)==3 and sys.argv[1]=='--make-fixture':
        fixture(Path(sys.argv[2]))
    else:
        raise SystemExit(pytest.main([__file__,'-q']))
