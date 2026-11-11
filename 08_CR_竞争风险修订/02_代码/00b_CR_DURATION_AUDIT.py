"""Aggregate audit of the two duration gates after full-table patient ranking."""
from pathlib import Path
import io, os, subprocess
import pandas as pd
R=Path(os.environ['IR_REVISION_ROOT'])
psql=os.environ.get('IR_PSQL','psql')
rows=[]
for cohort,db,port,icu,hosp,sid,age,join in [
 ('MIMIC-IV','mimiciv',5432,'mimiciv_icu','mimiciv_hosp','stay_id','d.admission_age','LEFT JOIN mimiciv_derived.icustay_detail d ON d.stay_id=i.stay_id'),
 ('MIMIC-III','mimiciii',5433,'mimiciii','mimiciii','icustay_id','d.admission_age','LEFT JOIN mimiciii.icustay_detail d ON d.icustay_id=i.icustay_id'),
 ('NWICU','nwicu',5432,'icu','hosp','stay_id','d.age','LEFT JOIN derived.cohort_all_icu_lab_final_select d ON d.stay_id=i.stay_id')]:
 sql=f'''WITH ranked AS (SELECT *,row_number() OVER(PARTITION BY subject_id ORDER BY intime,{sid}) AS rn FROM {icu}.icustays),
 eligible AS (SELECT i.intime,i.outtime,a.dischtime,a.hospital_expire_flag,
   (COALESCE(a.hospital_expire_flag,0)=1 OR (a.deathtime>=a.admittime AND a.deathtime<=a.dischtime)) AS death
   FROM ranked i JOIN {hosp}.admissions a ON a.hadm_id=i.hadm_id {join}
   WHERE i.rn=1 AND {age}>=18), flags AS (SELECT *,outtime>=intime+interval '24 hours' AS icu24,
   dischtime>=intime+interval '24 hours' AS hospital24 FROM eligible)
 SELECT count(*) AS adult_n,count(*) FILTER(WHERE icu24) AS after_icu24,
 count(*) FILTER(WHERE icu24 AND hospital24) AS after_both24,
 count(*) FILTER(WHERE icu24 AND hospital24 IS NOT TRUE) AS excluded_hospital24,
 count(*) FILTER(WHERE icu24 AND hospital24 IS NOT TRUE AND death) AS excluded_hospital24_deaths,
 count(*) FILTER(WHERE icu24 AND hospital24 IS NOT TRUE AND dischtime<outtime) AS excluded_discharge_before_icuout,
 count(*) FILTER(WHERE icu24 AND dischtime IS NULL) AS excluded_missing_discharge FROM flags'''
 p=subprocess.run([psql,'-X','-v','ON_ERROR_STOP=1','-h','localhost','-p',str(port),'-d',db,'--csv','-c',sql],check=True,capture_output=True,text=True)
 row=pd.read_csv(io.StringIO(p.stdout)).iloc[0].to_dict();row={'cohort':cohort,**row};rows.append(row)
a=pd.DataFrame(rows);a['excluded_hospital24_death_percent']=100*a.excluded_hospital24_deaths/a.excluded_hospital24
out=R/'08_CR_竞争风险修订/03_汇总结果/00_CR_duration_gate_audit.csv';a.to_csv(out,index=False)
print(a.to_string(index=False))
