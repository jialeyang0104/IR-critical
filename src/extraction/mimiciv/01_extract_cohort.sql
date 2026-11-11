/*
MIMIC-IV v3.1 local PostgreSQL restricted extract.

Run without creating database objects:
  psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d mimiciv \
    --csv -f "01_extract_cohort.sql" \
    > "<revision v3>/data/restricted/mimiciv_analysis_extract.csv"

The four exposure-defining laboratory measurements are exported as raw
component summaries. SHR, TyG and AIP are deliberately not calculated here.
*/

WITH patient_order AS (
  SELECT stay_id, ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY intime, stay_id) AS rn
  FROM mimiciv_icu.icustays
), eligible AS (
  SELECT
    d.subject_id,
    d.hadm_id,
    d.stay_id,
    d.gender,
    d.admission_age,
    d.race,
    d.admittime,
    d.dischtime,
    a.deathtime,
    p.dod,
    COALESCE(a.deathtime, CASE WHEN d.hospital_expire_flag = 1 THEN d.dischtime END) AS death_event_time,
    d.icu_intime,
    d.icu_outtime,
    d.hospital_expire_flag,
    d.los_hospital AS los_hospital_recorded,
    d.los_icu AS los_icu_recorded,
    EXTRACT(EPOCH FROM (d.dischtime - d.icu_intime)) / 86400.0 AS hospital_time_days,
    EXTRACT(EPOCH FROM (d.icu_outtime - d.icu_intime)) / 86400.0 AS icu_time_days
  FROM mimiciv_derived.icustay_detail AS d
  INNER JOIN mimiciv_hosp.admissions AS a
    ON a.hadm_id = d.hadm_id
  INNER JOIN mimiciv_hosp.patients AS p
    ON p.subject_id = d.subject_id
  INNER JOIN patient_order r ON r.stay_id=d.stay_id
  WHERE r.rn=1
    AND d.admission_age >= 18
    /* Exact, inclusive, dual 24-hour eligibility from ICU admission. */
    AND d.icu_outtime >= d.icu_intime + INTERVAL '24 hours'
    AND d.dischtime >= d.icu_intime + INTERVAL '24 hours'
),
outcomes AS (
  SELECT
    e.*,
    CASE WHEN e.hospital_expire_flag = 1 THEN 1 ELSE 0 END AS hospital_outcome,
    CASE
      WHEN e.death_event_time IS NOT NULL
       AND e.death_event_time >= e.icu_intime
       AND e.death_event_time <= e.icu_outtime THEN 1
      ELSE 0
    END AS icu_outcome,
    /* 28-day all-cause mortality: exact admission death time first; otherwise
       use the date-level patient DOD. DOD never changes hospital/ICU outcomes. */
    CASE
      WHEN e.deathtime IS NOT NULL THEN
        CASE WHEN e.deathtime >= e.icu_intime
                   AND e.deathtime <= e.icu_intime + INTERVAL '28 days' THEN 1 ELSE 0 END
      WHEN e.dod IS NOT NULL THEN
        CASE WHEN e.dod >= e.icu_intime::date
                   AND e.dod <= e.icu_intime::date + 28 THEN 1 ELSE 0 END
      ELSE 0
    END AS mortality_28d,
    CASE
      WHEN e.deathtime IS NOT NULL
       AND e.deathtime >= e.icu_intime
       AND e.deathtime <= e.icu_intime + INTERVAL '28 days'
        THEN GREATEST(EXTRACT(EPOCH FROM (e.deathtime - e.icu_intime)) / 86400.0, 1.0 / 86400.0)
      WHEN e.deathtime IS NULL
       AND e.dod IS NOT NULL
       AND e.dod >= e.icu_intime::date
       AND e.dod <= e.icu_intime::date + 28
        THEN GREATEST((e.dod - e.icu_intime::date)::numeric, 1.0 / 86400.0)
      ELSE 28.0
    END AS survival_time_28d,
    CASE
      WHEN e.deathtime IS NOT NULL
       AND e.deathtime >= e.icu_intime
       AND e.deathtime <= e.icu_intime + INTERVAL '28 days' THEN 'admissions.deathtime_exact'
      WHEN e.deathtime IS NULL
       AND e.dod IS NOT NULL
       AND e.dod >= e.icu_intime::date
       AND e.dod <= e.icu_intime::date + 28 THEN 'patients.dod_date'
      ELSE 'censored_at_28_days'
    END AS mortality_28d_time_source,
    LEAST(
      7.0,
      CASE
        WHEN e.death_event_time IS NOT NULL AND e.death_event_time >= e.icu_intime
          THEN EXTRACT(EPOCH FROM (e.death_event_time - e.icu_intime)) / 86400.0
        ELSE 7.0
      END
    ) AS survival_time_7d,
    LEAST(
      90.0,
      CASE
        WHEN e.death_event_time IS NOT NULL AND e.death_event_time >= e.icu_intime
          THEN EXTRACT(EPOCH FROM (e.death_event_time - e.icu_intime)) / 86400.0
        ELSE 90.0
      END
    ) AS survival_time_90d,
    CASE
      WHEN e.hospital_expire_flag = 1
        THEN EXTRACT(EPOCH FROM (e.death_event_time - e.icu_intime)) / 86400.0
      ELSE EXTRACT(EPOCH FROM (e.dischtime - e.icu_intime)) / 86400.0
    END AS survival_time_crude
  FROM eligible AS e
),
early_lab_events AS (
  SELECT
    o.stay_id,
    o.icu_intime,
    le.charttime,
    le.itemid,
    le.valuenum,
    le.valueuom
  FROM outcomes AS o
  /* LATERAL preserves use of the local (hadm_id, charttime, itemid) index. */
  CROSS JOIN LATERAL (
    SELECT x.charttime, x.itemid, x.valuenum, x.valueuom
    FROM mimiciv_hosp.labevents AS x
    WHERE x.hadm_id = o.hadm_id
      AND x.charttime >= o.icu_intime - INTERVAL '6 hours'
      AND x.charttime <= o.icu_intime + INTERVAL '24 hours'
      AND x.valuenum IS NOT NULL
      AND x.itemid IN (
     50852, 51631,                         /* HbA1c */
     50904,                                /* HDL-C */
     51000,                                /* triglycerides */
     50809, 50931, 52027, 52569,           /* blood glucose */
     51221, 51638, 51639, 52028,           /* hematocrit */
     50811, 51222, 51640,                  /* hemoglobin */
     51265, 53189,                         /* platelet count */
     51301, 51755, 51756,                  /* WBC */
     50862, 53085, 50868, 52500, 50882,
     50893, 50902, 52535, 50912, 52546,
     50930, 50971, 52610, 50983, 52623,
     51006, 52647, 51237, 51675, 51274,
     52921, 51275, 52923, 50861, 50878,
     50813, 52442, 53154
      )
    OFFSET 0
  ) AS le
),
routine_labs AS (
  SELECT
    stay_id,
    AVG(valuenum) FILTER (WHERE itemid IN (51221,51638,51639,52028)) AS hematocrit_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (50811,51222,51640)) AS hemoglobin_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (51265,53189)) AS platelets_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (51301,51755,51756)) AS wbc_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (50862,53085)) AS albumin_mean,
    AVG(valuenum) FILTER (WHERE itemid = 50930) AS globulin_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (50868,52500)) AS aniongap_mean,
    AVG(valuenum) FILTER (WHERE itemid = 50882) AS bicarbonate_mean,
    AVG(valuenum) FILTER (WHERE itemid = 51006 OR itemid = 52647) AS bun_mean,
    AVG(valuenum) FILTER (WHERE itemid = 50893) AS calcium_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (50902,52535)) AS chloride_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (50912,52546)) AS creatinine_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (50983,52623)) AS sodium_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (50971,52610)) AS potassium_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (51237,51675)) AS inr_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (51274,52921)) AS pt_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (51275,52923)) AS ptt_mean,
    AVG(valuenum) FILTER (WHERE itemid = 50861) AS alt_mean,
    AVG(valuenum) FILTER (WHERE itemid = 50878) AS ast_mean,
    AVG(valuenum) FILTER (WHERE itemid IN (50813,52442,53154)) AS lactate_mean
  FROM early_lab_events
  GROUP BY stay_id
),
exposure_long AS (
  SELECT
    stay_id,
    icu_intime,
    charttime,
    valuenum,
    CASE
      WHEN itemid IN (50809,50931,52027,52569) THEN 'glucose'
      WHEN itemid IN (50852,51631) THEN 'hba1c'
      WHEN itemid = 51000 THEN 'tg'
      WHEN itemid = 50904 THEN 'hdl'
    END AS analyte
  FROM early_lab_events
  WHERE itemid IN (50809,50931,52027,52569,50852,51631,51000,50904)
),
exposure_ranked AS (
  SELECT
    x.*,
    ROW_NUMBER() OVER (PARTITION BY stay_id, analyte ORDER BY charttime, valuenum) AS first_rn,
    ROW_NUMBER() OVER (
      PARTITION BY stay_id, analyte
      ORDER BY ABS(EXTRACT(EPOCH FROM (charttime - icu_intime))), charttime, valuenum
    ) AS closest_rn
  FROM exposure_long AS x
  WHERE analyte IS NOT NULL
),
exposure_labs AS (
  SELECT
    stay_id,
    AVG(valuenum) FILTER (WHERE analyte = 'glucose') AS glucose,
    AVG(valuenum) FILTER (WHERE analyte = 'hba1c') AS hba1c,
    AVG(valuenum) FILTER (WHERE analyte = 'tg') AS tg,
    AVG(valuenum) FILTER (WHERE analyte = 'hdl') AS hdl,
    COUNT(*) FILTER (WHERE analyte = 'glucose') AS glucose_n,
    COUNT(*) FILTER (WHERE analyte = 'hba1c') AS hba1c_n,
    COUNT(*) FILTER (WHERE analyte = 'tg') AS tg_n,
    COUNT(*) FILTER (WHERE analyte = 'hdl') AS hdl_n,
    MAX(valuenum) FILTER (WHERE analyte = 'glucose' AND first_rn = 1) AS glucose_first_value,
    MAX(charttime) FILTER (WHERE analyte = 'glucose' AND first_rn = 1) AS glucose_first_charttime,
    MAX(valuenum) FILTER (WHERE analyte = 'hba1c' AND first_rn = 1) AS hba1c_first_value,
    MAX(charttime) FILTER (WHERE analyte = 'hba1c' AND first_rn = 1) AS hba1c_first_charttime,
    MAX(valuenum) FILTER (WHERE analyte = 'tg' AND first_rn = 1) AS tg_first_value,
    MAX(charttime) FILTER (WHERE analyte = 'tg' AND first_rn = 1) AS tg_first_charttime,
    MAX(valuenum) FILTER (WHERE analyte = 'hdl' AND first_rn = 1) AS hdl_first_value,
    MAX(charttime) FILTER (WHERE analyte = 'hdl' AND first_rn = 1) AS hdl_first_charttime,
    MAX(valuenum) FILTER (WHERE analyte = 'glucose' AND closest_rn = 1) AS glucose_closest_value,
    MAX(charttime) FILTER (WHERE analyte = 'glucose' AND closest_rn = 1) AS glucose_closest_charttime,
    MAX(valuenum) FILTER (WHERE analyte = 'hba1c' AND closest_rn = 1) AS hba1c_closest_value,
    MAX(charttime) FILTER (WHERE analyte = 'hba1c' AND closest_rn = 1) AS hba1c_closest_charttime,
    MAX(valuenum) FILTER (WHERE analyte = 'tg' AND closest_rn = 1) AS tg_closest_value,
    MAX(charttime) FILTER (WHERE analyte = 'tg' AND closest_rn = 1) AS tg_closest_charttime,
    MAX(valuenum) FILTER (WHERE analyte = 'hdl' AND closest_rn = 1) AS hdl_closest_value,
    MAX(charttime) FILTER (WHERE analyte = 'hdl' AND closest_rn = 1) AS hdl_closest_charttime
  FROM exposure_ranked
  GROUP BY stay_id
),
diagnosis_flags AS (
  SELECT
    o.hadm_id,
    MAX(CASE WHEN (d.icd_version = 9 AND d.icd_code LIKE '401%')
                   OR (d.icd_version = 10 AND LEFT(d.icd_code,3) IN ('I10','I11','I15','I16')) THEN 1 ELSE 0 END) AS hypertension,
    MAX(CASE WHEN (d.icd_version = 9 AND d.icd_code LIKE '250%')
                   OR (d.icd_version = 10 AND LEFT(d.icd_code,3) IN ('E10','E11')) THEN 1 ELSE 0 END) AS diabetes,
    MAX(CASE WHEN (d.icd_version = 9 AND d.icd_code IN ('5722','5728'))
                   OR (d.icd_version = 10 AND LEFT(d.icd_code,3) IN ('K70','K71','K72')) THEN 1 ELSE 0 END) AS liver_failure,
    MAX(CASE WHEN (d.icd_version = 9 AND (LEFT(d.icd_code,3) IN ('585','586') OR LEFT(d.icd_code,5) IN ('40301','40311','40391','40402','40403','40412','40413','40492','40493')))
                   OR (d.icd_version = 10 AND LEFT(d.icd_code,3) IN ('N18','N19')) THEN 1 ELSE 0 END) AS renal_failure,
    MAX(CASE WHEN (d.icd_version = 9 AND d.icd_code = 'V5867')
                   OR (d.icd_version = 10 AND d.icd_code = 'Z794') THEN 1 ELSE 0 END) AS longterm_insulin,
    MAX(CASE WHEN d.icd_version = 10 AND d.icd_code = 'Z7984' THEN 1 ELSE 0 END) AS longterm_antidiabetics
  FROM outcomes AS o
  LEFT JOIN mimiciv_hosp.diagnoses_icd AS d
    ON d.hadm_id = o.hadm_id
  GROUP BY o.hadm_id
),
ventilation_24h AS (
  SELECT
    o.stay_id,
    MAX(CASE WHEN v.ventilation_status IS NOT NULL AND v.ventilation_status <> 'None' THEN 1 ELSE 0 END) AS ventilation
  FROM outcomes AS o
  LEFT JOIN mimiciv_derived.ventilation AS v
    ON v.stay_id = o.stay_id
   AND v.starttime <= o.icu_intime + INTERVAL '24 hours'
   AND COALESCE(v.endtime, v.starttime) >= o.icu_intime - INTERVAL '6 hours'
  GROUP BY o.stay_id
),
vasopressor_24h AS (
  SELECT
    o.stay_id,
    MAX(CASE WHEN va.dopamine IS NOT NULL OR va.epinephrine IS NOT NULL
                   OR va.norepinephrine IS NOT NULL OR va.phenylephrine IS NOT NULL
                   OR va.vasopressin IS NOT NULL OR va.dobutamine IS NOT NULL
                   OR va.milrinone IS NOT NULL THEN 1 ELSE 0 END) AS vasopressor
  FROM outcomes AS o
  LEFT JOIN mimiciv_derived.vasoactive_agent AS va
    ON va.stay_id = o.stay_id
   AND va.starttime <= o.icu_intime + INTERVAL '24 hours'
   AND COALESCE(va.endtime, va.starttime) >= o.icu_intime - INTERVAL '6 hours'
  GROUP BY o.stay_id
),
eligible_prescriptions AS MATERIALIZED (
  /* Restrict the 20M-row medication table to eligible admissions before regex. */
  SELECT p.hadm_id, p.starttime, p.stoptime, p.drug, p.route
  FROM mimiciv_hosp.prescriptions AS p
  INNER JOIN (SELECT DISTINCT hadm_id FROM outcomes) AS h
    ON h.hadm_id = p.hadm_id
),
relevant_prescriptions AS MATERIALIZED (
  SELECT hadm_id, starttime, stoptime, drug, route
  FROM eligible_prescriptions
  WHERE LOWER(drug) ~ '(dexamethasone|hydrocortisone|methylprednisolone|prednisone|prednisolone|cortisone|betamethasone|triamcinolone|fludrocortisone|interferon|interleukin|filgrastim|pegfilgrastim|sargramostim|immunoglobulin|immune globulin|azathioprine|tacrolimus|everolimus|sirolimus|cyclosporine|mycophenolate|abatacept|tofacitinib|anti.thymocyte|antithymocyte|atezolizumab|blinatumomab|ipilimumab|nivolumab|odronextamab|pembrolizumab|heparin)'
),
prescription_flags AS (
  SELECT
    o.stay_id,
    MAX(CASE WHEN LOWER(p.drug) ~ '(dexamethasone|hydrocortisone|methylprednisolone|prednisone|prednisolone|cortisone|betamethasone|triamcinolone|fludrocortisone)' THEN 1 ELSE 0 END) AS glucocorticoids_1d,
    MAX(CASE WHEN LOWER(p.drug) ~ '(interferon|interleukin|filgrastim|pegfilgrastim|sargramostim|immunoglobulin|immune globulin)' THEN 1 ELSE 0 END) AS immune_stimulate,
    MAX(CASE WHEN LOWER(p.drug) ~ '(azathioprine|tacrolimus|everolimus|sirolimus|cyclosporine|mycophenolate|abatacept|tofacitinib|anti.thymocyte|antithymocyte|atezolizumab|blinatumomab|ipilimumab|nivolumab|odronextamab|pembrolizumab)' THEN 1 ELSE 0 END) AS immune_suppress,
    MAX(CASE WHEN LOWER(p.drug) LIKE '%heparin%'
                  AND LOWER(COALESCE(p.route,'')) ~ '(subcut|subq|sq|sc)' THEN 1 ELSE 0 END) AS heparin_1d
  FROM outcomes AS o
  LEFT JOIN relevant_prescriptions AS p
    ON p.hadm_id = o.hadm_id
   AND p.starttime <= o.icu_intime + INTERVAL '24 hours'
   AND COALESCE(p.stoptime, p.starttime) >= o.icu_intime - INTERVAL '6 hours'
  GROUP BY o.stay_id
)
SELECT
  o.subject_id,
  o.hadm_id,
  o.stay_id,
  o.gender,
  o.admission_age,
  o.race,
  o.admittime,
  o.dischtime,
  o.deathtime,
  o.dod,
  o.death_event_time,
  o.icu_intime,
  o.icu_outtime,
  o.los_hospital_recorded,
  o.los_icu_recorded,
  o.hospital_time_days AS los_hospital,
  o.icu_time_days AS los_icu,
  o.hospital_time_days,
  o.icu_time_days,
  o.hospital_outcome,
  o.icu_outcome,
  o.mortality_28d,
  o.mortality_28d_time_source,
  o.survival_time_7d,
  o.survival_time_28d,
  o.survival_time_28d AS time_28d_days,
  o.survival_time_90d,
  o.survival_time_crude,
  fs.sofa AS sofa_score,
  g.gcs_min,
  vs.heart_rate_mean,
  vs.sbp_mean,
  vs.dbp_mean,
  vs.mbp_mean,
  vs.resp_rate_mean,
  vs.temperature_mean,
  vs.spo2_mean,
  rl.hematocrit_mean,
  rl.hemoglobin_mean,
  rl.platelets_mean,
  rl.wbc_mean,
  rl.albumin_mean,
  rl.globulin_mean,
  rl.aniongap_mean,
  rl.bicarbonate_mean,
  rl.bun_mean,
  rl.calcium_mean,
  rl.chloride_mean,
  rl.creatinine_mean,
  rl.sodium_mean,
  rl.potassium_mean,
  rl.inr_mean,
  rl.pt_mean,
  rl.ptt_mean,
  rl.alt_mean,
  rl.ast_mean,
  rl.lactate_mean,
  h.height,
  w.weight,
  u.urineoutput,
  rr.dialysis_active,
  rr.dialysis_type,
  ap.apsiii,
  lo.lods,
  oa.oasis,
  sa.sapsii,
  si.sirs,
  ch.myocardial_infarct,
  ch.congestive_heart_failure,
  ch.peripheral_vascular_disease,
  ch.cerebrovascular_disease,
  ch.dementia,
  ch.chronic_pulmonary_disease,
  ch.rheumatic_disease,
  ch.peptic_ulcer_disease,
  ch.mild_liver_disease,
  ch.diabetes_without_cc,
  ch.diabetes_with_cc,
  ch.paraplegia,
  ch.renal_disease,
  ch.malignant_cancer,
  ch.severe_liver_disease,
  ch.metastatic_solid_tumor,
  ch.aids,
  ch.charlson_comorbidity_index,
  COALESCE(df.hypertension,0) AS hypertension,
  COALESCE(df.diabetes, GREATEST(ch.diabetes_without_cc,ch.diabetes_with_cc), 0) AS diabetes,
  COALESCE(df.liver_failure,0) AS liver_failure,
  COALESCE(df.renal_failure,0) AS renal_failure,
  COALESCE(df.longterm_insulin,0) AS longterm_insulin,
  COALESCE(df.longterm_antidiabetics,0) AS longterm_antidiabetics,
  COALESCE(vn.ventilation,0) AS ventilation,
  COALESCE(va.vasopressor,0) AS vasopressor,
  COALESCE(pf.glucocorticoids_1d,0) AS glucocorticoids_1d,
  COALESCE(pf.immune_stimulate,0) AS immune_stimulate,
  COALESCE(pf.immune_suppress,0) AS immune_suppress,
  COALESCE(pf.heparin_1d,0) AS heparin_1d,
  /* Raw exposure-component summaries; no index is authoritative in SQL. */
  el.glucose,
  el.hba1c,
  el.tg,
  el.hdl,
  el.glucose_n,
  el.hba1c_n,
  el.tg_n,
  el.hdl_n,
  el.glucose_first_value,
  el.glucose_first_charttime,
  el.hba1c_first_value,
  el.hba1c_first_charttime,
  el.tg_first_value,
  el.tg_first_charttime,
  el.hdl_first_value,
  el.hdl_first_charttime,
  el.glucose_closest_value,
  el.glucose_closest_charttime,
  el.hba1c_closest_value,
  el.hba1c_closest_charttime,
  el.tg_closest_value,
  el.tg_closest_charttime,
  el.hdl_closest_value,
  el.hdl_closest_charttime
FROM outcomes AS o
LEFT JOIN mimiciv_derived.first_day_sofa AS fs ON fs.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.first_day_gcs AS g ON g.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.first_day_vitalsign AS vs ON vs.stay_id = o.stay_id
LEFT JOIN routine_labs AS rl ON rl.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.first_day_height AS h ON h.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.first_day_weight AS w ON w.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.first_day_urine_output AS u ON u.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.first_day_rrt AS rr ON rr.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.apsiii AS ap ON ap.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.lods AS lo ON lo.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.oasis AS oa ON oa.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.sapsii AS sa ON sa.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.sirs AS si ON si.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.charlson AS ch ON ch.hadm_id = o.hadm_id
LEFT JOIN diagnosis_flags AS df ON df.hadm_id = o.hadm_id
LEFT JOIN ventilation_24h AS vn ON vn.stay_id = o.stay_id
LEFT JOIN vasopressor_24h AS va ON va.stay_id = o.stay_id
LEFT JOIN prescription_flags AS pf ON pf.stay_id = o.stay_id
LEFT JOIN exposure_labs AS el ON el.stay_id = o.stay_id
ORDER BY o.stay_id;
