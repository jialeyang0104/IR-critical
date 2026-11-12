from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
import csv,hashlib,re,json,subprocess
import pandas as pd
R=Path(__file__).resolve().parents[2];D=R/'results/analysis'
files=list(D.glob('01_CR_*_regression_results.csv'));assert len(files)==12
risk=pd.read_csv(D/'00_CR_risk_set_audit.csv')
for p in files:
 a=pd.read_csv(p);assert len(a)==108 and set(a.exposure)=={'SHR','TyG','AIP'}
 e=risk[(risk.cohort==a.cohort.iloc[0])&(risk.time_origin==a.time_origin.iloc[0])].set_index('outcome').included_n
 assert (a.n==a.outcome.map(e)).all() and (a.m==0).all()
 assert ((a.conf_low>0)&(a.conf_low<=a.estimate)&(a.estimate<=a.conf_high)).all()
 assert ((a.outcome=='mortality_28d')==(a.effect_measure=='HR')).all()
public=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=R).decode().split('\0')
for name in public:
 if not name:continue
 p=R/name
 if not p.is_file():continue
 assert name.isascii(),name
 assert not name.startswith(('data/restricted/','reports/internal/','src/reporting/')),name
 assert p.suffix.lower() not in {'.rds','.rdata','.mids','.env','.log','.docx','.sqlite'},name
 if p.suffix.lower()=='.csv':
  with p.open(encoding='utf-8-sig',newline='')as f:h=next(csv.reader(f))
  assert not set(x.lower()for x in h)&{'subject_id','patient_id','hadm_id','stay_id','icustay_id','source_row','original_source_row','charttime','intime','outtime','icu_intime','icu_outtime','admittime','dischtime','deathtime','dod'},name
 if p.suffix.lower()in{'.r','.py','.sql','.js','.md','.txt','.json'} and p.name!='01_validate_public.py':
  text=p.read_text();assert '/Users/'not in text,name
  assert not re.search(r'(?:postgres(?:ql)?://[^\s:]+:[^\s@]+@|(?:password|api_key|access_token)\s*[=:]\s*["\'][^"\']+["\'])',text,re.I),name
for line in (R/'reports/reproducibility/CHECKSUMS.sha256').read_text().splitlines():
 digest,name=line.split('  ',1);assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
f=R/'deliverables/figures';m=pd.read_csv(f/'figure_export_manifest.csv');assert len(m)==16
for n in m.figure:
 for ext in ['.png','.pdf','.svg']:assert(f/(n+ext)).is_file()
primary=pd.concat([pd.read_csv(D/f'01_CR_{co}_AO_M5_regression_results.csv') for co in ['M4','M3','NW']])
primary=primary[~(primary.cohort.eq('NWICU')&primary.outcome.eq('mortality_28d'))]
keys=['cohort','exposure','outcome','model','term']
columns=keys+['exposure_scale','n','events','effect_measure','estimate','conf_low','conf_high','p_value']
for name,scale,count in [('Figure2_All_cohorts','per_unit',72),('FigureS01_Tertile_all_cohorts','tertile',144)]:
 actual=pd.read_csv(f/'source_data'/f'{name}_source.csv')
 expected=primary[primary.exposure_scale.eq(scale)]
 assert len(actual)==len(expected)==count and not actual.duplicated(keys).any(),name
 pd.testing.assert_frame_equal(actual[columns].sort_values(keys).reset_index(drop=True),
  expected[columns].sort_values(keys).reset_index(drop=True),check_dtype=False,check_exact=False,rtol=1e-12,atol=1e-14)
 assert all(color in (f/f'{name}.svg').read_text().lower() for color in ['#4e79a7','#f28e2b','#2e8b57'])
assert not list(f.glob('FigureS01_Primary_MIMICIV.*'))
a=pd.read_csv(D/'05_CR_M4_AO_subgroup_interactions.csv');assert len(a)==27 and a.convergence.all() and a.n.eq(3056).all()
s=pd.read_csv(D/'05_CR_M4_AO_M5_exploratory_subgroups.csv');assert len(s)==54 and s.groupby('exposure').index_sd.nunique().eq(1).all()
t=pd.read_csv(D/'10_CR_table_result_figure_mapping.csv');assert t.table.tolist()==['S'+str(i)for i in range(1,17)]
assert t.loc[t.table.eq('S6'),'rows'].item()==288
spec=spec_from_file_location('baseline_validation',R/'src/validation/03_validate_baseline_tables.py')
baseline_validation=module_from_spec(spec);spec.loader.exec_module(baseline_validation)
baseline_validation.validate_aggregates()
for n in ['11_CR_paired_perSD_effect_contrasts.csv','12_CR_patient_bootstrap_CV_increments.csv','12_CR_paired_CV_increment_contrasts.csv']:
 a=pd.read_csv(D/n);assert(a.bootstrap_success==200).all()
print('PASS: patient-level matrices, risk-set denominators, complete covariates, 16 figures, 200 bootstraps, file-field audit and checksums.')
