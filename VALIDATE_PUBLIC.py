from pathlib import Path
import csv,hashlib,re,json
import pandas as pd
R=Path(__file__).resolve().parent;D=R/'08_CR_竞争风险修订/03_汇总结果'
files=list(D.glob('01_CR_*_regression_results.csv'));assert len(files)==12
risk=pd.read_csv(D/'00_CR_risk_set_audit.csv')
for p in files:
 a=pd.read_csv(p);assert len(a)==108 and set(a.exposure)=={'SHR','TyG','AIP'}
 e=risk[(risk.cohort==a.cohort.iloc[0])&(risk.time_origin==a.time_origin.iloc[0])].set_index('outcome').included_n
 assert (a.n==a.outcome.map(e)).all() and (a.m==0).all()
 assert ((a.conf_low>0)&(a.conf_low<=a.estimate)&(a.estimate<=a.conf_high)).all()
 assert ((a.outcome=='mortality_28d')==(a.effect_measure=='HR')).all()
for p in R.rglob('*'):
 if not p.is_file():continue
 if '.git' in p.relative_to(R).parts:continue
 assert '.restricted_input'not in p.parts and p.suffix.lower()not in{'.rds','.rdata','.mids','.env','.log','.docx','.sqlite'}
 if p.suffix.lower()=='.csv':
  with p.open(encoding='utf-8-sig',newline='')as f:h=next(csv.reader(f))
  assert not set(x.lower()for x in h)&{'subject_id','patient_id','hadm_id','stay_id','icustay_id','source_row','original_source_row','charttime','intime','outtime','icu_intime','icu_outtime','admittime','dischtime','deathtime','dod'}
 if p.suffix.lower()in{'.r','.py','.sql','.js','.md','.txt','.json'} and p.name!='VALIDATE_PUBLIC.py':
  s=p.read_text();assert '/Users/'not in s
for line in (R/'CHECKSUMS.sha256').read_text().splitlines():
 digest,name=line.split('  ',1);assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
f=R/'08_CR_竞争风险修订/01_最终投稿文件/figures';m=pd.read_csv(f/'figure_export_manifest.csv');assert len(m)==16
for n in m.figure:
 for ext in ['.png','.pdf','.svg']:assert(f/(n+ext)).is_file()
a=pd.read_csv(D/'05_CR_M4_AO_subgroup_interactions.csv');assert len(a)==27 and a.convergence.all() and a.n.eq(3056).all()
s=pd.read_csv(D/'05_CR_M4_AO_M5_exploratory_subgroups.csv');assert len(s)==54 and s.groupby('exposure').index_sd.nunique().eq(1).all()
t=pd.read_csv(D/'10_CR_table_result_figure_mapping.csv');assert t.table.tolist()==['S'+str(i)for i in range(1,16)]
for n in ['11_CR_paired_perSD_effect_contrasts.csv','12_CR_patient_bootstrap_CV_increments.csv','12_CR_paired_CV_increment_contrasts.csv']:
 a=pd.read_csv(D/n);assert(a.bootstrap_success==200).all()
print('PASS: patient-level matrices, risk-set denominators, complete covariates, 16 figures, 200 bootstraps, file-field audit and checksums.')
