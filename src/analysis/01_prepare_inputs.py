"""Rebuild restricted patient-level inputs from preserved extracts and source ICU tables.

Ranking precedes every eligibility filter. No database objects are changed.
Only aggregate QA is written outside data/restricted.
"""
from pathlib import Path
import hashlib
import io
import json
import os
import subprocess
import pandas as pd

ROOT = Path(os.environ['IR_REVISION_ROOT'])
RES = ROOT / 'results/analysis'
RES.mkdir(parents=True, exist_ok=True)
(ROOT / 'reports/internal').mkdir(parents=True, exist_ok=True)
PSQL = os.environ.get('IR_PSQL', 'psql')

def query(db, port, sql):
    p = subprocess.run([PSQL, '-X', '-v', 'ON_ERROR_STOP=1', '-h', 'localhost',
                        '-p', str(port), '-d', db, '--csv', '-c', sql],
                       check=True, capture_output=True, text=True)
    return pd.read_csv(io.StringIO(p.stdout))

def category(code, version):
    code = str(code).replace('.', '').upper()
    if code in ('NAN', '', 'NONE'):
        return 'Unknown'
    if version == 9:
        if not code[:3].isdigit():
            return 'Other'
        n = int(code[:3])
        for lo, hi, label in [(1,139,'Infectious'),(140,239,'Neoplasm'),(240,279,'Endocrine/metabolic'),
                              (390,459,'Circulatory'),(460,519,'Respiratory'),(520,579,'Digestive'),
                              (580,629,'Genitourinary'),(800,999,'Injury/poisoning')]:
            if lo <= n <= hi:
                return label
        return 'Other'
    return {'A':'Infectious','B':'Infectious','C':'Neoplasm','E':'Endocrine/metabolic',
            'I':'Circulatory','J':'Respiratory','K':'Digestive','N':'Genitourinary',
            'S':'Injury/poisoning','T':'Injury/poisoning'}.get(code[0], 'Other')

flows, audit, manifest = [], [], []
specs = [('MIMIC-IV','mimiciv',5432,'mimiciv_icu','mimiciv_hosp','stay_id','mimiciv_analysis.csv'),
         ('MIMIC-III','mimiciii',5433,'mimiciii','mimiciii','icustay_id','mimiciii_external_validation.csv'),
         ('NWICU','nwicu',5432,'icu','hosp','stay_id','nwicu_external_validation.csv')]
for cohort, db, port, icu, hosp, sid, filename in specs:
    diagnosis = ('dx.icd9_code AS primary_code, 9 AS primary_version' if db == 'mimiciii'
                 else 'dx.icd_code AS primary_code, dx.icd_version AS primary_version')
    codecol = 'icd9_code' if db == 'mimiciii' else 'icd_code'
    source = 'i.dbsource' if db == 'mimiciii' else "NULL::text AS dbsource"
    if db == 'mimiciii':
        ages = 'dd.admission_age'
        agejoin = f'LEFT JOIN mimiciii.icustay_detail dd ON dd.icustay_id=i.{sid}'
    elif db == 'mimiciv':
        ages = 'dd.admission_age'
        agejoin = 'LEFT JOIN mimiciv_derived.icustay_detail dd ON dd.stay_id=i.stay_id'
    else:
        ages = 'dd.age'
        agejoin = 'LEFT JOIN derived.cohort_all_icu_lab_final_select dd ON dd.stay_id=i.stay_id'
    meta = query(db, port, f'''
        WITH ranked AS (SELECT *, ROW_NUMBER() OVER
          (PARTITION BY subject_id ORDER BY intime, {sid}) AS patient_icu_rank
          FROM {icu}.icustays), primary_dx AS (
          SELECT DISTINCT ON (hadm_id) * FROM {hosp}.diagnoses_icd
          ORDER BY hadm_id,seq_num NULLS LAST,{codecol})
        SELECT i.subject_id,i.hadm_id,i.{sid} AS stay_key,i.patient_icu_rank,
          i.intime,i.outtime,a.dischtime,a.admission_type,i.first_careunit,
          {source}, {ages} AS source_age, {diagnosis}
        FROM ranked i JOIN {hosp}.admissions a ON a.hadm_id=i.hadm_id
        {agejoin}
        LEFT JOIN primary_dx dx ON dx.hadm_id=i.hadm_id
        ''')
    assert not meta.stay_key.duplicated().any()
    src = ROOT / 'data/restricted' / filename
    raw = pd.read_csv(src, low_memory=False)
    raw['original_source_row'] = range(1, len(raw)+1)
    # Match one source stay to its rank in the complete database, never rank selected rows.
    newcols = ['patient_icu_rank','admission_type','first_careunit','dbsource','primary_code','primary_version']
    raw = raw.drop(columns=[c for c in newcols if c in raw], errors='ignore')
    joined = raw.merge(meta[['stay_key'] + newcols], left_on=sid, right_on='stay_key',
                       how='left', validate='one_to_one', sort=False)
    assert joined.patient_icu_rank.notna().all()
    selected = joined.loc[joined.patient_icu_rank.eq(1)].copy()
    assert not selected.subject_id.duplicated().any()
    selected['diagnosis_category'] = [category(c,v) for c,v in zip(selected.primary_code, selected.primary_version)]
    sepsis9 = "icd9_code IN ('99591','99592','78552')"
    head9 = "substring(icd9_code,1,3) BETWEEN '800' AND '804' OR substring(icd9_code,1,3) BETWEEN '850' AND '854'"
    if db == 'mimiciii':
        sepsis, head = sepsis9, head9
    else:
        sepsis = "(icd_version=9 AND icd_code IN ('99591','99592','78552')) OR (icd_version=10 AND (replace(icd_code,'.','') ~ '^(A40|A41|R652|R572)'))"
        head = "(icd_version=9 AND (substring(icd_code,1,3) BETWEEN '800' AND '804' OR substring(icd_code,1,3) BETWEEN '850' AND '854')) OR (icd_version=10 AND icd_code ~ '^S0[0-9]')"
    flags = query(db, port, f'''SELECT hadm_id,MAX(({sepsis})::int) AS coded_sepsis,
       MAX(({head})::int) AS coded_head_trauma FROM {hosp}.diagnoses_icd GROUP BY hadm_id''')
    selected = selected.merge(flags, on='hadm_id', how='left', validate='many_to_one')
    # Preserved original age and an explicit indicator distinguish deidentification from missingness.
    agecol = 'admission_age' if 'admission_age' in selected else 'age'
    original_age = selected['age_original'] if 'age_original' in selected else selected[agecol]
    selected['age_deidentified'] = ((original_age >= 300) & (db == 'mimiciii')).astype(int)
    selected['age_original'] = original_age
    selected.loc[selected.age_deidentified.eq(1), agecol] = 90
    if db == 'nwicu':
        covid = query(db, port, """SELECT hadm_id,
          MAX((replace(upper(icd_code),'.','')='U071' AND icd_version=10)::int) AS covid_u071,
          MAX((replace(upper(icd_code),'.','')='B9729' AND icd_version=10)::int) AS coronavirus_b9729
          FROM hosp.diagnoses_icd GROUP BY hadm_id""")
        # Exact first 24 hours; only documented administration, systemic routes.
        steroids = query(db, port, """SELECT DISTINCT i.stay_id,1 AS steroid_24h
          FROM icu.icustays i JOIN hosp.emar e ON e.hadm_id=i.hadm_id
            AND e.charttime>=i.intime AND e.charttime<i.intime+interval '24 hours'
          JOIN hosp.prescriptions p ON p.hadm_id=e.hadm_id AND p.pharmacy_id=e.pharmacy_id
          WHERE lower(e.event_txt) IN ('confirmed','applied')
            AND lower(e.medication) ~ '(hydrocortisone|dexamethasone|methylprednisolone|prednisone|prednisolone|cortisone|betamethasone)'
            AND upper(p.route) IN ('IV','IV PUSH','IVPB','PO','ORAL','NG','GT','IM','INTRAVENOUS','INTRAMUSCULAR','ENTERAL')""")
        selected = selected.merge(covid, on='hadm_id', how='left', validate='many_to_one')
        selected['diagnosis_available'] = selected.covid_u071.notna().astype(int)
        selected[['covid_u071','coronavirus_b9729']] = selected[['covid_u071','coronavirus_b9729']].fillna(0)
        selected = selected.merge(steroids, on='stay_id', how='left', validate='one_to_one')
        selected['steroid_24h'] = selected.steroid_24h.fillna(0).astype(int)
    selected = selected.sort_values('original_source_row').drop(columns=['stay_key','primary_code','primary_version'])
    target = ROOT / 'data/restricted' / f'{db}_patient_level.csv'
    tmp = target.with_suffix('.tmp')
    selected.to_csv(tmp, index=False)
    os.chmod(tmp, 0o600)
    tmp.replace(target)
    first = meta.loc[meta.patient_icu_rank.eq(1)].copy()
    adult = first.source_age.ge(18)
    duration = ((pd.to_datetime(first.outtime)-pd.to_datetime(first.intime)).dt.total_seconds().ge(86400) &
                (pd.to_datetime(first.dischtime)-pd.to_datetime(first.intime)).dt.total_seconds().ge(86400))
    flows.append(dict(cohort=cohort,all_icu_stays=len(meta),first_patient_stays=len(first),
                      source_covariates_available=int(first.source_age.notna().sum()),
                      excluded_source_covariates=int(first.source_age.isna().sum()),
                      adult_first_stays=int(adult.sum()),eligible_pre_laboratory=len(selected),
                      excluded_nonfirst=len(meta)-len(first),excluded_age=int((first.source_age.notna() & ~adult).sum()),
                      excluded_duration=int((adult & ~duration).sum())))
    audit.append(dict(cohort=cohort,source_rows=len(raw),source_unique_patients=raw.subject_id.nunique(),
                      patient_rows=len(selected),unique_patients=selected.subject_id.nunique(),
                      duplicate_patients=int(selected.subject_id.duplicated().sum()),
                      deidentified_age_rows=int(selected.age_deidentified.sum())))
    manifest.append(dict(cohort=cohort,input=filename,input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
                         output=target.name,output_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                         rule='Full ICU table: subject_id / intime / stay ID; rank=1 before eligibility',
                         main_elderly_age=90,sensitivity_elderly_age=91.4))
pd.DataFrame(flows).to_csv(RES/'00_CR_patient_flow.csv', index=False)
pd.DataFrame(audit).to_csv(RES/'00_CR_patient_selection_audit.csv', index=False)
(ROOT/'reports/internal/00_CR_input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(pd.DataFrame(audit).to_string(index=False))
