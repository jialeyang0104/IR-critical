-- Read-only supplement for observed baseline tables. No regression input is changed.
-- Vital signs: first 24 hours. Routine laboratories: -6 to +24 hours.
WITH ranked AS (
    SELECT i.*, row_number() OVER (PARTITION BY subject_id ORDER BY intime, icustay_id) AS rn
    FROM mimiciii.icustays i
), eligible AS MATERIALIZED (
    SELECT r.*, a.admittime, a.dischtime
    FROM ranked r
    JOIN mimiciii.icustay_detail d USING (icustay_id)
    JOIN mimiciii.admissions a ON a.hadm_id = r.hadm_id
    WHERE r.rn = 1 AND d.admission_age >= 18
      AND r.outtime >= r.intime + interval '24 hours'
      AND a.dischtime >= r.intime + interval '24 hours'
), labs AS (
    SELECT e.icustay_id,
        avg(l.valuenum) FILTER (WHERE l.itemid IN (51300,51301)) AS wbc_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 51265) AS platelets_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid IN (50811,51222)) AS hemoglobin_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 51006) AS bun_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 50912) AS creatinine_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 50983) AS sodium_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 50971) AS potassium_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 51237) AS inr_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 51274) AS pt_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 50861) AS alt_mean,
        avg(l.valuenum) FILTER (WHERE l.itemid = 50878) AS ast_mean
    FROM eligible e
    LEFT JOIN mimiciii.labevents l ON l.hadm_id = e.hadm_id
      AND l.charttime BETWEEN e.intime - interval '6 hours' AND e.intime + interval '24 hours'
      AND l.itemid IN (51300,51301,51265,50811,51222,51006,50912,50983,50971,51237,51274,50861,50878)
      AND l.valuenum IS NOT NULL AND l.valuenum <> 9999999
    GROUP BY e.icustay_id
), diagnoses AS (
    SELECT e.icustay_id,
        max((substr(d.icd9_code,1,3) IN ('410','412'))::int) AS myocardial_infarct,
        max((substr(d.icd9_code,1,3)='428'
            OR substr(d.icd9_code,1,5) IN ('39891','40201','40211','40291','40401','40403','40411','40413','40491','40493')
            OR substr(d.icd9_code,1,4) BETWEEN '4254' AND '4259')::int) AS congestive_heart_failure,
        max((substr(d.icd9_code,1,3) BETWEEN '490' AND '505'
            OR substr(d.icd9_code,1,4) IN ('4168','4169','5064','5081','5088'))::int) AS chronic_pulmonary_disease,
        max((d.icd9_code='V5867')::int) AS longterm_insulin
    FROM eligible e LEFT JOIN mimiciii.diagnoses_icd d ON d.hadm_id=e.hadm_id
    GROUP BY e.icustay_id
), steroids AS (
    SELECT e.icustay_id, max((lower(p.drug) ~
        '(dexamethasone|hydrocortisone|methylprednisolone|prednisone|prednisolone|cortisone|betamethasone|triamcinolone|fludrocortisone)')::int)
        AS glucocorticoids_1d
    FROM eligible e LEFT JOIN mimiciii.prescriptions p ON p.hadm_id=e.hadm_id
      AND p.startdate <= e.intime + interval '24 hours'
      AND coalesce(p.enddate,p.startdate) >= e.intime - interval '6 hours'
    GROUP BY e.icustay_id
), vaso AS (
    SELECT e.icustay_id, max((v.icustay_id IS NOT NULL)::int) AS vasopressor_day1
    FROM eligible e LEFT JOIN (
        SELECT icustay_id,starttime,endtime FROM mimiciii.vasopressor_durations
        UNION ALL SELECT icustay_id,starttime,endtime FROM mimiciii.vasopressin_durations
    ) v ON v.icustay_id=e.icustay_id
      AND v.starttime <= e.intime + interval '24 hours'
      AND coalesce(v.endtime,v.starttime) >= e.intime - interval '6 hours'
    GROUP BY e.icustay_id
)
SELECT e.subject_id,e.hadm_id,e.icustay_id,
    v.heartrate_mean AS heart_rate_mean,v.resprate_mean AS resp_rate_mean,
    v.meanbp_mean AS mbp_mean,v.tempc_mean AS temperature_mean,
    l.wbc_mean,l.platelets_mean,l.hemoglobin_mean,l.bun_mean,l.creatinine_mean,
    l.sodium_mean,l.potassium_mean,l.inr_mean,l.pt_mean,l.alt_mean,l.ast_mean,
    dx.myocardial_infarct,dx.congestive_heart_failure,dx.chronic_pulmonary_disease,
    dx.longterm_insulin,
    -- ICD-9 has no equivalent dedicated code for long-term oral antidiabetic use.
    NULL::int AS longterm_antidiabetics,
    coalesce(s.glucocorticoids_1d,0) AS glucocorticoids_1d,
    va.vasopressor_day1
FROM eligible e
LEFT JOIN mimiciii.vitals_first_day v USING (icustay_id)
LEFT JOIN labs l USING (icustay_id)
LEFT JOIN diagnoses dx USING (icustay_id)
LEFT JOIN steroids s USING (icustay_id)
LEFT JOIN vaso va USING (icustay_id)
ORDER BY e.icustay_id;
