"""Check observed baseline, missingness, selection, and survival-status tables."""
from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT/'results/analysis'
COHORTS = ['MIMIC-IV', 'MIMIC-III', 'NWICU']


def read(name):
    return pd.read_csv(D/name, keep_default_na=False)


def validate_aggregates():
    table1 = read('00_CR_ALL_table1_observed_baseline.csv')
    missing = read('00_CR_ALL_Table1_variable_missingness.csv')
    selection = read('00_CR_ALL_included_excluded_comparison.csv')
    survival = read('00_CR_M4_survivor_nonsurvivor_baseline.csv')
    counts = read('00_CR_ALL_included_excluded_counts.csv')
    fields = read('00_CR_ALL_baseline_source_fields.csv')
    assert len(table1)==42 and table1.variable.nunique()==42
    assert len(survival)==42 and survival.variable.tolist()==table1.variable.tolist()
    assert survival.cohort.eq('MIMIC-IV').all()
    assert survival.overall.tolist()==table1.MIMICIV.tolist()
    assert survival.total_n.eq(3056).all()
    assert survival.group_1_n.eq(2735).all() and survival.group_2_n.eq(321).all()
    assert (survival.group_1_n+survival.group_2_n).eq(survival.total_n).all()
    assert (survival.group_1_nonmissing_n+survival.group_2_nonmissing_n).eq(survival.overall_nonmissing_n).all()
    assert survival.overall_nonmissing_n.between(0,3056).all()
    assert survival.group_1_nonmissing_n.between(0,2735).all()
    assert survival.group_2_nonmissing_n.between(0,321).all()
    assert survival.group_1_label.eq('Hospital survivor').all()
    assert survival.group_2_label.eq('Hospital non-survivor').all()
    assert survival.loc[survival.type.eq('continuous'),'comparison_test'].eq('Mann-Whitney U').all()
    assert survival.loc[survival.type.eq('binary'),'comparison_test'].isin(['Pearson chi-square','Fisher exact']).all()
    assert len(missing)==126 and len(selection)==144 and len(fields)==144
    assert not missing.duplicated(['cohort','variable']).any()
    assert not selection.duplicated(['cohort','variable']).any()
    assert counts.cohort.tolist()==COHORTS
    assert counts.included_n.tolist()==[3056,1133,1556]
    assert counts.excluded_n.tolist()==[48508,30989,15795]
    assert counts.eligible_n.eq(counts.included_n+counts.excluded_n).all()
    for co, column in zip(COHORTS,['MIMICIV','MIMICIII','NWICU']):
        m=missing[missing.cohort.eq(co)]
        s=selection[selection.cohort.eq(co)]
        assert m.variable.tolist()==table1.variable.tolist()
        assert s.variable.tolist()[:42]==table1.variable.tolist()
        assert s.group_1.iloc[:42].tolist()==table1[column].tolist()
        if co=='MIMIC-IV':
            assert survival.overall_nonmissing_n.tolist()==m.nonmissing_n.astype(int).tolist()
        n=int(counts.loc[counts.cohort.eq(co),'included_n'].iloc[0])
        assert m.total_n.eq(n).all() and s.group_1_n.eq(n).all()
        unavailable=m.variable.eq('longterm_antidiabetics') & (co=='MIMIC-III')
        observed=m[~unavailable]
        mn=observed.missing_n.astype(int); non=observed.nonmissing_n.astype(int)
        assert (mn+non).eq(n).all()
        assert np.allclose(observed.missing_percent.astype(float),100*mn/n)
        for r in observed.itertuples():
            assert r.display==f'{int(r.missing_n)} ({float(r.missing_percent):.1f}%)'
        assert observed.loc[observed.variable.isin(['gender','age','sofa','charlson','Diabetes']), 'missing_n'].astype(int).eq(0).all()
        assert observed.loc[observed.variable.isin(['TG','Glucose','HbA1c','HDL','SHR','TyG','AIP']),'missing_n'].astype(int).eq(0).all()
    unavailable=missing[missing.cohort.eq('MIMIC-III')&missing.variable.eq('longterm_antidiabetics')].iloc[0]
    assert unavailable.display=='NA' and unavailable.missing_n==''
    nw=missing[missing.cohort.eq('NWICU')&missing.variable.isin(['gcs','INR'])]
    assert len(nw)==2 and nw.missing_n.astype(int).eq(1556).all()
    for data in [table1,selection,survival]:
        numeric=pd.to_numeric(data.p_value,errors='coerce')
        assert numeric.dropna().between(0,1).all()
        for p,display in zip(numeric,data.p_display):
            assert display==('NA' if pd.isna(p) else '<0.001' if p<.001 else f'{p:.3f}')
    assert table1.comparison_scope[table1.variable.isin(['gcs','INR'])].eq('MIMIC-IV;MIMIC-III').all()
    return {'status':'PASS','baseline_variables':42,'missingness_cohort_variable_rows':126,
            'survival_comparisons':42,'survivors':2735,'non_survivors':321,
            'selection_comparisons':144,'cohort_order':COHORTS,'included_n':[3056,1133,1556],
            'excluded_n':[48508,30989,15795]}


if __name__=='__main__':
    print(json.dumps(validate_aggregates(),indent=2))
