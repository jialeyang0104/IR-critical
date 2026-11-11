\set ON_ERROR_STOP on
\timing on

-- 本地MIMIC-III宽表聚合在PostgreSQL默认work_mem下会产生大量临时分批文件。
-- 仅调整当前会话，不修改数据库配置或持久对象。
SET work_mem = '256MB';

DROP TABLE IF EXISTS pg_temp.mimiciii_critical_ill_patients_export;
DROP TABLE IF EXISTS pg_temp.mimiciii_base_export;

CREATE TEMP TABLE mimiciii_base_export AS
WITH patient_order AS (
    SELECT icustay_id, ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY intime, icustay_id) AS rn
    FROM mimiciii.icustays
), base0 AS (
    SELECT
        'mimiciii'::text AS source_database,
        i.subject_id,
        i.hadm_id,
        i.icustay_id,
        i.gender,
        i.ethnicity,
        i.ethnicity_grouped,
        i.admittime,
        i.dischtime,
        i.intime AS icu_intime,
        i.outtime AS icu_outtime,
        a.deathtime,
        p.dod,
        CASE WHEN i.admission_age >= 300 THEN 90 ELSE i.admission_age END AS admission_age,
        (i.admission_age >= 300)::int AS age_deidentified,
        i.admission_age AS age_original,
        i.los_hospital,
        i.los_icu,
        i.hospstay_seq,
        i.icustay_seq,
        i.first_hosp_stay,
        i.first_icu_stay,
        i.hospital_expire_flag,
        CASE WHEN i.hospital_expire_flag = 1 THEN 1 ELSE 0 END AS hospital_outcome,
        CASE
            WHEN a.deathtime IS NOT NULL
             AND a.deathtime >= i.intime
             AND a.deathtime <= i.outtime THEN 1
            ELSE 0
        END AS icu_outcome,
        CASE
            WHEN a.deathtime IS NOT NULL
             AND a.deathtime >= i.intime
             AND a.deathtime <= i.intime + INTERVAL '28 days' THEN 1
            WHEN a.deathtime IS NULL
             AND p.dod IS NOT NULL
             AND p.dod::date >= i.intime::date
             AND p.dod::date <= (i.intime + INTERVAL '28 days')::date THEN 1
            ELSE 0
        END AS outcome_28d,
        CASE
            WHEN a.deathtime IS NOT NULL AND a.deathtime >= i.intime THEN 'admissions.deathtime'
            WHEN p.dod IS NOT NULL AND p.dod::date >= i.intime::date THEN 'patients.dod_date'
            ELSE 'none'
        END AS death_time_source
    FROM mimiciii.icustay_detail i
    INNER JOIN mimiciii.admissions a
      ON a.subject_id = i.subject_id AND a.hadm_id = i.hadm_id
    INNER JOIN mimiciii.patients p
      ON p.subject_id = i.subject_id
    INNER JOIN patient_order r ON r.icustay_id=i.icustay_id
    WHERE r.rn=1
      AND i.admission_age >= 18
      -- 入组定义为自ICU入科起，ICU和住院观察时间均至少24小时。
      -- 不使用los_icu/los_hospital进行入组判断，避免整数化或舍入造成边界误分。
      AND i.outtime >= i.intime + INTERVAL '24 hours'
      AND i.dischtime >= i.intime + INTERVAL '24 hours'
)
SELECT
    b.*,
    CASE
      WHEN b.hospital_outcome = 1 THEN COALESCE(b.deathtime, b.dischtime)
      ELSE b.dischtime
    END AS hospital_event_or_censor_time,
    CASE WHEN b.icu_outcome = 1 THEN b.deathtime ELSE b.icu_outtime END
      AS icu_event_or_censor_time,
    CASE
      WHEN b.outcome_28d = 1 AND b.deathtime IS NOT NULL THEN b.deathtime
      WHEN b.outcome_28d = 1 AND b.dod IS NOT NULL THEN b.dod::date::timestamp
      ELSE b.icu_intime + INTERVAL '28 days'
    END AS day28_event_or_censor_time
FROM base0 b;

CREATE INDEX mimiciii_base_export_hadm_idx
  ON mimiciii_base_export (hadm_id, icu_intime);
CREATE UNIQUE INDEX mimiciii_base_export_stay_idx
  ON mimiciii_base_export (icustay_id);
ANALYZE mimiciii_base_export;

DROP TABLE IF EXISTS pg_temp.mimiciii_novel_labs_export;
CREATE TEMP TABLE mimiciii_novel_labs_export AS
SELECT
    b.subject_id,
    b.hadm_id,
    b.icustay_id,
    AVG(le.valuenum) FILTER (WHERE le.itemid IN (50809, 50931)) AS glucose_mean,
    AVG(le.valuenum) FILTER (WHERE le.itemid = 50852) AS hba1c_mean,
    AVG(le.valuenum) FILTER (WHERE le.itemid = 51000) AS triglycerides_mean,
    AVG(le.valuenum) FILTER (WHERE le.itemid = 50904) AS hdl_mean,
    AVG(le.valuenum) FILTER (WHERE le.itemid IN (50905, 50906)) AS ldl_mean,
    AVG(le.valuenum) FILTER (WHERE le.itemid = 50907) AS total_cholesterol_mean,
    COUNT(*) FILTER (WHERE le.itemid IN (50809, 50931)) AS glucose_n,
    COUNT(*) FILTER (WHERE le.itemid = 50852) AS hba1c_n,
    COUNT(*) FILTER (WHERE le.itemid = 51000) AS triglycerides_n,
    COUNT(*) FILTER (WHERE le.itemid = 50904) AS hdl_n
FROM mimiciii_base_export b
LEFT JOIN mimiciii.labevents le
  ON le.hadm_id = b.hadm_id
 AND le.charttime >= b.icu_intime - INTERVAL '6 hours'
 AND le.charttime <= b.icu_intime + INTERVAL '24 hours'
 AND le.itemid IN (50809, 50931, 50852, 51000, 50904, 50905, 50906, 50907)
 AND le.valuenum IS NOT NULL
 AND le.valuenum > 0
GROUP BY b.subject_id, b.hadm_id, b.icustay_id;

CREATE UNIQUE INDEX mimiciii_novel_labs_export_stay_idx
  ON mimiciii_novel_labs_export (icustay_id);
ANALYZE mimiciii_novel_labs_export;

CREATE TEMP TABLE mimiciii_critical_ill_patients_export AS
WITH sofa_final AS (
    SELECT
        subject_id,
        hadm_id,
        icustay_id,
        sofa,
        respiration,
        coagulation,
        liver,
        cardiovascular,
        cns,
        renal
    FROM mimiciii.sofa
),
diagnosis_flags AS (
    SELECT
        subject_id,
        hadm_id,
        MAX(CASE WHEN REPLACE(icd9_code, '.', '') ~ '^40[1-5]' THEN 1 ELSE 0 END) AS icd_hypertension,
        MAX(CASE WHEN REPLACE(icd9_code, '.', '') ~ '^250' THEN 1 ELSE 0 END) AS icd_diabetes,
        MAX(CASE WHEN REPLACE(icd9_code, '.', '') ~ '^(570|571|572|573)' THEN 1 ELSE 0 END) AS icd_liver_disease,
        MAX(CASE WHEN REPLACE(icd9_code, '.', '') ~ '^58[4-6]' THEN 1 ELSE 0 END) AS icd_renal_failure,
        MAX(CASE WHEN REPLACE(icd9_code, '.', '') ~ '^(200|201|202|203|204|205|206|207|208)' THEN 1 ELSE 0 END) AS icd_hematologic_malignancy,
        MAX(CASE WHEN REPLACE(icd9_code, '.', '') ~ '^(14[0-9]|15[0-9]|16[0-9]|17[0-6]|18[0-9]|19[0-9]|20[0-9]|V10)' THEN 1 ELSE 0 END) AS icd_malignancy,
        MAX(CASE WHEN REPLACE(icd9_code, '.', '') ~ '^(V42|9968)' THEN 1 ELSE 0 END) AS icd_transplant
    FROM mimiciii.diagnoses_icd
    GROUP BY subject_id, hadm_id
),
charlson_diag AS (
    SELECT
        b.icustay_id,
        REGEXP_REPLACE(UPPER(di.icd9_code), '[^A-Z0-9]', '', 'g') AS clean_icd
    FROM mimiciii_base_export b
    INNER JOIN mimiciii.diagnoses_icd di
      ON di.hadm_id = b.hadm_id
),
charlson_components AS (
    SELECT
        icustay_id,
        MAX(CASE WHEN SUBSTR(clean_icd, 1, 3) IN ('410', '412') THEN 1 ELSE 0 END) AS myocardial_infarct,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) = '428'
              OR SUBSTR(clean_icd, 1, 5) IN ('39891', '40201', '40211', '40291', '40401', '40403', '40411', '40413', '40491', '40493')
              OR SUBSTR(clean_icd, 1, 4) BETWEEN '4254' AND '4259'
            THEN 1 ELSE 0 END) AS congestive_heart_failure,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) IN ('440', '441')
              OR SUBSTR(clean_icd, 1, 4) IN ('0930', '4373', '4471', '5571', '5579', 'V434')
              OR SUBSTR(clean_icd, 1, 4) BETWEEN '4431' AND '4439'
            THEN 1 ELSE 0 END) AS peripheral_vascular_disease,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) BETWEEN '430' AND '438'
              OR SUBSTR(clean_icd, 1, 5) = '36234'
            THEN 1 ELSE 0 END) AS cerebrovascular_disease,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) = '290'
              OR SUBSTR(clean_icd, 1, 4) IN ('2941', '3312')
            THEN 1 ELSE 0 END) AS dementia,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) BETWEEN '490' AND '505'
              OR SUBSTR(clean_icd, 1, 4) IN ('4168', '4169', '5064', '5081', '5088')
            THEN 1 ELSE 0 END) AS chronic_pulmonary_disease,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) = '725'
              OR SUBSTR(clean_icd, 1, 4) IN ('4465', '7100', '7101', '7102', '7103', '7104', '7140', '7141', '7142', '7148')
            THEN 1 ELSE 0 END) AS rheumatic_disease,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) IN ('531', '532', '533', '534')
            THEN 1 ELSE 0 END) AS peptic_ulcer_disease,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) IN ('570', '571')
              OR SUBSTR(clean_icd, 1, 4) IN ('0706', '0709', '5733', '5734', '5738', '5739', 'V427')
              OR SUBSTR(clean_icd, 1, 5) IN ('07022', '07023', '07032', '07033', '07044', '07054')
            THEN 1 ELSE 0 END) AS mild_liver_disease,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 4) IN ('2500', '2501', '2502', '2503', '2508', '2509')
            THEN 1 ELSE 0 END) AS diabetes_without_cc,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 4) IN ('2504', '2505', '2506', '2507')
            THEN 1 ELSE 0 END) AS diabetes_with_cc,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) IN ('342', '343')
              OR SUBSTR(clean_icd, 1, 4) IN ('3341', '3440', '3441', '3442', '3443', '3444', '3445', '3446', '3449')
            THEN 1 ELSE 0 END) AS paraplegia,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) IN ('582', '585', '586', 'V56')
              OR SUBSTR(clean_icd, 1, 4) IN ('5880', 'V420', 'V451')
              OR SUBSTR(clean_icd, 1, 4) BETWEEN '5830' AND '5837'
              OR SUBSTR(clean_icd, 1, 5) IN ('40301', '40311', '40391', '40402', '40403', '40412', '40413', '40492', '40493')
            THEN 1 ELSE 0 END) AS renal_disease,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) BETWEEN '140' AND '172'
              OR SUBSTR(clean_icd, 1, 4) BETWEEN '1740' AND '1958'
              OR SUBSTR(clean_icd, 1, 3) BETWEEN '200' AND '208'
              OR SUBSTR(clean_icd, 1, 4) = '2386'
            THEN 1 ELSE 0 END) AS malignant_cancer,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 4) IN ('4560', '4561', '4562')
              OR SUBSTR(clean_icd, 1, 4) BETWEEN '5722' AND '5728'
            THEN 1 ELSE 0 END) AS severe_liver_disease,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) IN ('196', '197', '198', '199')
            THEN 1 ELSE 0 END) AS metastatic_solid_tumor,
        MAX(CASE
            WHEN SUBSTR(clean_icd, 1, 3) IN ('042', '043', '044')
            THEN 1 ELSE 0 END) AS aids
    FROM charlson_diag
    GROUP BY icustay_id
),
charlson_scores AS (
    SELECT
        b.icustay_id,
        (
            CASE
                WHEN b.admission_age IS NULL OR b.admission_age <= 50 THEN 0
                WHEN b.admission_age <= 60 THEN 1
                WHEN b.admission_age <= 70 THEN 2
                WHEN b.admission_age <= 80 THEN 3
                ELSE 4
            END
            + COALESCE(cc.myocardial_infarct, 0)
            + COALESCE(cc.congestive_heart_failure, 0)
            + COALESCE(cc.peripheral_vascular_disease, 0)
            + COALESCE(cc.cerebrovascular_disease, 0)
            + COALESCE(cc.dementia, 0)
            + COALESCE(cc.chronic_pulmonary_disease, 0)
            + COALESCE(cc.rheumatic_disease, 0)
            + COALESCE(cc.peptic_ulcer_disease, 0)
            + GREATEST(COALESCE(cc.mild_liver_disease, 0), 3 * COALESCE(cc.severe_liver_disease, 0))
            + GREATEST(COALESCE(cc.diabetes_without_cc, 0), 2 * COALESCE(cc.diabetes_with_cc, 0))
            + 2 * COALESCE(cc.paraplegia, 0)
            + 2 * COALESCE(cc.renal_disease, 0)
            + GREATEST(2 * COALESCE(cc.malignant_cancer, 0), 6 * COALESCE(cc.metastatic_solid_tumor, 0))
            + 6 * COALESCE(cc.aids, 0)
        )::integer AS charlson_comorbidity_index
    FROM mimiciii_base_export b
    LEFT JOIN charlson_components cc
      ON cc.icustay_id = b.icustay_id
),
vasopressor_agg AS (
    SELECT
        icustay_id,
        1 AS vasopressor_icu,
        SUM(duration_hours) AS vasopressor_duration_hours
    FROM (
        SELECT icustay_id, duration_hours FROM mimiciii.vasopressor_durations
        UNION ALL
        SELECT icustay_id, duration_hours FROM mimiciii.vasopressin_durations
    ) v
    GROUP BY icustay_id
),
-- 仅保留本研究分析、基线审计和预测性能需要的字段，避免全库无关宽表聚合。
analysis_ready AS (SELECT 1 AS ready)
SELECT
    b.source_database,
    b.subject_id,
    b.hadm_id,
    b.icustay_id,
    b.gender,
    b.ethnicity,
    b.ethnicity_grouped,
    b.admittime,
    b.dischtime,
    b.icu_intime,
    b.icu_outtime,
    b.deathtime,
    b.dod,
    b.admission_age,
    b.age_deidentified,
    b.age_original,
    b.los_hospital,
    b.los_icu,
    b.hospstay_seq,
    b.icustay_seq,
    b.first_hosp_stay,
    b.first_icu_stay,
    b.hospital_expire_flag,
    b.hospital_outcome,
    b.icu_outcome,
    b.outcome_28d,
    b.death_time_source,
    b.hospital_event_or_censor_time,
    b.icu_event_or_censor_time,
    b.day28_event_or_censor_time,
    GREATEST(EXTRACT(EPOCH FROM (b.hospital_event_or_censor_time - b.icu_intime))/86400.0, 0.001)
      AS hospital_time_days,
    GREATEST(EXTRACT(EPOCH FROM (b.icu_event_or_censor_time - b.icu_intime))/86400.0, 0.001)
      AS icu_time_days,
    CASE
      WHEN b.outcome_28d = 1 AND b.deathtime IS NOT NULL THEN
        GREATEST(EXTRACT(EPOCH FROM (b.deathtime - b.icu_intime))/86400.0, 0.001)
      WHEN b.outcome_28d = 1 AND b.dod IS NOT NULL THEN
        GREATEST((b.dod::date - b.icu_intime::date)::double precision, 0.001)
      ELSE 28.0
    END AS time_28d_days,
    sf.sofa,
    sf.respiration AS sofa_respiration,
    sf.coagulation AS sofa_coagulation,
    sf.liver AS sofa_liver,
    sf.cardiovascular AS sofa_cardiovascular,
    sf.cns AS sofa_cns,
    sf.renal AS sofa_renal,
    g.mingcs,
    g.gcsmotor,
    g.gcsverbal,
    g.gcseyes,
    g.endotrachflag,
    vent.vent AS ventilation_day1,
    vent.vent AS ventilation,
    COALESCE(va.vasopressor_icu, 0) AS vasopressor_icu,
    COALESCE(va.vasopressor_icu, 0) AS vasopressor,
    va.vasopressor_duration_hours,
    ap.apsiii,
    ap.apsiii_prob,
    oa.oasis,
    oa.oasis_prob,
    sa.sapsii,
    sa.sapsii_prob,
    sirs.sirs,
    cs.charlson_comorbidity_index,
    df.icd_hypertension,
    df.icd_diabetes,
    df.icd_liver_disease,
    df.icd_liver_disease AS liver_failure,
    df.icd_renal_failure,
    df.icd_renal_failure AS renal_failure,
    df.icd_hematologic_malignancy,
    df.icd_malignancy,
    df.icd_transplant,
    nl.glucose_mean AS first_day_glucose_mean,
    nl.hba1c_mean AS first_day_hba1c_mean,
    nl.triglycerides_mean AS first_day_triglycerides_mean,
    nl.hdl_mean AS first_day_hdl_mean,
    nl.ldl_mean AS first_day_ldl_mean,
    nl.total_cholesterol_mean AS first_day_total_cholesterol_mean,
    nl.glucose_n,
    nl.hba1c_n,
    nl.triglycerides_n,
    nl.hdl_n,
    CASE
        WHEN nl.glucose_mean > 0
         AND nl.hba1c_mean > 0
         AND (28.7 * nl.hba1c_mean - 46.7) > 0
            THEN nl.glucose_mean / (28.7 * nl.hba1c_mean - 46.7)
    END AS shr,
    CASE
        WHEN nl.glucose_mean > 0 AND nl.triglycerides_mean > 0
            THEN LN(nl.triglycerides_mean * nl.glucose_mean / 2.0)
    END AS tyg,
    CASE
        WHEN nl.triglycerides_mean > 0 AND nl.hdl_mean > 0
            THEN LN((nl.triglycerides_mean / 88.57) / (nl.hdl_mean / 38.67)) / LN(10.0)
    END AS aip
FROM mimiciii_base_export b
LEFT JOIN sofa_final sf
  ON sf.icustay_id = b.icustay_id
LEFT JOIN mimiciii.gcs_first_day g
  ON g.icustay_id = b.icustay_id
LEFT JOIN mimiciii.ventilation_first_day vent
  ON vent.icustay_id = b.icustay_id
LEFT JOIN vasopressor_agg va
  ON va.icustay_id = b.icustay_id
LEFT JOIN mimiciii.apsiii ap
  ON ap.icustay_id = b.icustay_id
LEFT JOIN mimiciii.oasis oa
  ON oa.icustay_id = b.icustay_id
LEFT JOIN mimiciii.sapsii sa
  ON sa.icustay_id = b.icustay_id
LEFT JOIN mimiciii.sirs sirs
  ON sirs.icustay_id = b.icustay_id
LEFT JOIN diagnosis_flags df
  ON df.subject_id = b.subject_id
 AND df.hadm_id = b.hadm_id
LEFT JOIN charlson_scores cs
  ON cs.icustay_id = b.icustay_id
LEFT JOIN mimiciii_novel_labs_export nl
  ON nl.icustay_id = b.icustay_id
ORDER BY b.subject_id, b.hadm_id, b.icustay_id;

-- 受限原始事件文件：仅保留ICU入科前6 h至入科后24 h的四项暴露实验室事件。
DROP TABLE IF EXISTS pg_temp.mimiciii_exposure_lab_events_export;
CREATE TEMP TABLE mimiciii_exposure_lab_events_export AS
SELECT
    e.subject_id,
    e.hadm_id,
    e.icustay_id,
    le.row_id AS lab_event_id,
    le.charttime,
    EXTRACT(EPOCH FROM (le.charttime - e.icu_intime))/3600.0 AS hours_from_icu_admission,
    CASE
      WHEN le.itemid IN (50809, 50931) THEN 'glucose'
      WHEN le.itemid = 50852 THEN 'hba1c'
      WHEN le.itemid = 51000 THEN 'triglycerides'
      WHEN le.itemid = 50904 THEN 'hdl'
    END AS analyte,
    le.itemid,
    le.valuenum,
    le.valueuom
FROM mimiciii_critical_ill_patients_export e
INNER JOIN mimiciii.labevents le
  ON le.hadm_id = e.hadm_id
 AND le.charttime >= e.icu_intime - INTERVAL '6 hours'
 AND le.charttime <= e.icu_intime + INTERVAL '24 hours'
 AND le.itemid IN (50809, 50931, 50852, 51000, 50904)
WHERE le.valuenum IS NOT NULL AND le.valuenum > 0;

DROP TABLE IF EXISTS pg_temp.mimiciii_exposure_lab_definitions_export;
CREATE TEMP TABLE mimiciii_exposure_lab_definitions_export AS
WITH ranked AS (
  SELECT
    x.*,
    ROW_NUMBER() OVER (
      PARTITION BY icustay_id, analyte ORDER BY charttime, lab_event_id
    ) AS rn_first,
    ROW_NUMBER() OVER (
      PARTITION BY icustay_id, analyte
      ORDER BY ABS(hours_from_icu_admission), charttime, lab_event_id
    ) AS rn_closest
  FROM mimiciii_exposure_lab_events_export x
)
SELECT
  subject_id,
  hadm_id,
  icustay_id,
  analyte,
  MAX(valuenum) FILTER (WHERE rn_first = 1) AS first_available_value,
  MAX(charttime) FILTER (WHERE rn_first = 1) AS first_available_charttime,
  MAX(hours_from_icu_admission) FILTER (WHERE rn_first = 1) AS first_available_hours,
  MAX(valuenum) FILTER (WHERE rn_closest = 1) AS closest_to_icu_value,
  MAX(charttime) FILTER (WHERE rn_closest = 1) AS closest_to_icu_charttime,
  MAX(hours_from_icu_admission) FILTER (WHERE rn_closest = 1) AS closest_to_icu_hours,
  AVG(valuenum) AS early_window_mean,
  COUNT(*) AS early_window_n
FROM ranked
GROUP BY subject_id, hadm_id, icustay_id, analyte;

-- 运行时通过 psql 变量传入受限输出路径：
-- psql -v mimiciii_output_csv='/authorized/path/mimiciii_external_validation.csv' -f this_script.sql
-- 使用\o和COPY TO STDOUT，确保psql变量作为输出文件名正确展开。
\o :mimiciii_output_csv
COPY mimiciii_critical_ill_patients_export TO STDOUT WITH (FORMAT CSV, HEADER true);
\o

\o :mimiciii_lab_events_csv
COPY mimiciii_exposure_lab_events_export TO STDOUT WITH (FORMAT CSV, HEADER true);
\o

\o :mimiciii_lab_definitions_csv
COPY mimiciii_exposure_lab_definitions_export TO STDOUT WITH (FORMAT CSV, HEADER true);
\o

SELECT
    COUNT(*) AS exported_rows,
    COUNT(*) FILTER (WHERE hospital_outcome = 1) AS hospital_deaths,
    COUNT(*) FILTER (WHERE icu_outcome = 1) AS icu_deaths,
    COUNT(*) FILTER (WHERE outcome_28d = 1) AS deaths_28d,
    COUNT(DISTINCT icustay_id) AS unique_icustays,
    COUNT(*) FILTER (WHERE deathtime IS NULL AND hospital_outcome = 1) AS hospital_deaths_without_deathtime,
    COUNT(*) FILTER (WHERE charlson_comorbidity_index IS NULL) AS charlson_missing,
    COUNT(*) FILTER (WHERE first_day_glucose_mean IS NOT NULL AND first_day_hba1c_mean IS NOT NULL) AS shr_complete,
    COUNT(*) FILTER (WHERE first_day_glucose_mean IS NOT NULL AND first_day_triglycerides_mean IS NOT NULL) AS tyg_complete,
    COUNT(*) FILTER (WHERE first_day_triglycerides_mean IS NOT NULL AND first_day_hdl_mean IS NOT NULL) AS aip_complete,
    COUNT(*) FILTER (
        WHERE first_day_glucose_mean IS NOT NULL
          AND first_day_hba1c_mean IS NOT NULL
          AND first_day_triglycerides_mean IS NOT NULL
          AND first_day_hdl_mean IS NOT NULL
    ) AS all_three_complete
FROM mimiciii_critical_ill_patients_export;
