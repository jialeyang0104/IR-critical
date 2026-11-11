\set ON_ERROR_STOP on

DROP TABLE IF EXISTS pg_temp.nwicu_external_validation_export;

CREATE TEMP TABLE nwicu_external_validation_export AS
WITH ranked_icu_stays AS (
    SELECT
        i.subject_id,
        i.stay_id,
        ROW_NUMBER() OVER (
            PARTITION BY i.subject_id
            ORDER BY i.intime, i.stay_id
        ) AS icu_stay_sequence
    FROM icu.icustays i
), eligible AS (
    SELECT
        c.subject_id,
        c.hadm_id,
        c.stay_id,
        c.icu_intime,
        c.icu_outtime,
        a.admittime,
        a.dischtime,
        a.deathtime,
        p.dod,
        a.hospital_expire_flag,
        c.sex,
        c.age,
        r.icu_stay_sequence,
        c.charlson_comorbidity_index,
        GREATEST(COALESCE(c.charlson_diabetes_without_cc, 0),
                 COALESCE(c.charlson_diabetes_with_cc, 0)) AS diabetes,
        -- NWICU uses locally adapted SOFA; its component substitutions must be reported.
        ms.sofa_adapted_24hours AS sofa,
        c.charlson_mild_liver_disease,
        c.charlson_severe_liver_disease,
        c.charlson_renal_disease,
        c.itd_used_ventilation_icu AS ventilation_icu,
        c.itd_used_vasopressor_icu AS vasopressor_icu
    FROM derived.cohort_all_icu_lab_final_select c
    INNER JOIN ranked_icu_stays r
      ON r.subject_id = c.subject_id
     AND r.stay_id = c.stay_id
     AND r.icu_stay_sequence = 1
    LEFT JOIN hosp.admissions a
      ON a.subject_id = c.subject_id AND a.hadm_id = c.hadm_id
    INNER JOIN hosp.patients p
      ON p.subject_id = c.subject_id
    LEFT JOIN derived.nwicu_sofa_day1_adapted ms
      ON ms.stay_id=c.stay_id
    WHERE c.age >= 18
      -- 仅保留患者数据库记录中的首次ICU住院；随后应用精确双24小时条件。
      -- ICU和住院观察时间均自ICU入科起至少24小时；避免派生LOS舍入造成边界误分。
      AND c.icu_outtime >= c.icu_intime + INTERVAL '24 hours'
      AND a.dischtime >= c.icu_intime + INTERVAL '24 hours'
), lab_agg AS (
    SELECT
      e.stay_id,
      AVG(le.valuenum) FILTER (WHERE le.itemid IN (100001, 100045)) AS glucose,
      AVG(le.valuenum) FILTER (WHERE le.itemid = 100079) AS hba1c,
      AVG(le.valuenum) FILTER (WHERE le.itemid = 100056) AS triglycerides,
      AVG(le.valuenum) FILTER (WHERE le.itemid = 100086) AS hdl,
      COUNT(*) FILTER (WHERE le.itemid IN (100001, 100045)) AS glucose_n,
      COUNT(*) FILTER (WHERE le.itemid = 100079) AS hba1c_n,
      COUNT(*) FILTER (WHERE le.itemid = 100056) AS triglycerides_n,
      COUNT(*) FILTER (WHERE le.itemid = 100086) AS hdl_n
    FROM eligible e
    LEFT JOIN hosp.labevents le
      ON le.hadm_id = e.hadm_id
     AND le.charttime >= e.icu_intime - INTERVAL '6 hours'
     AND le.charttime <= e.icu_intime + INTERVAL '24 hours'
     AND le.itemid IN (100001, 100045, 100079, 100056, 100086)
     AND le.valuenum IS NOT NULL
     AND le.valuenum > 0
    GROUP BY e.stay_id
), base AS (
SELECT
    e.*,
    l.glucose,
    l.hba1c,
    l.triglycerides,
    l.hdl,
    l.glucose_n,
    l.hba1c_n,
    l.triglycerides_n,
    l.hdl_n
FROM eligible e
LEFT JOIN lab_agg l ON l.stay_id = e.stay_id
), outcome_flags AS (
SELECT
    base.*,
    CASE
      WHEN COALESCE(hospital_expire_flag, 0) = 1
        OR (deathtime IS NOT NULL AND deathtime >= admittime AND deathtime <= dischtime)
      THEN 1 ELSE 0
    END AS hospital_outcome,
    CASE
      WHEN deathtime IS NOT NULL AND deathtime >= icu_intime AND deathtime <= icu_outtime
      THEN 1 ELSE 0
    END AS icu_outcome,
    CASE
      WHEN deathtime IS NOT NULL AND deathtime >= icu_intime
       AND deathtime <= icu_intime + INTERVAL '28 days'
      THEN 1 ELSE 0
    END AS outcome_28d_exact,
    CASE
      WHEN deathtime IS NOT NULL AND deathtime >= icu_intime
       AND deathtime <= icu_intime + INTERVAL '28 days' THEN 1
      WHEN deathtime IS NULL AND dod IS NOT NULL
       AND dod >= icu_intime::date
       AND dod <= (icu_intime + INTERVAL '28 days')::date THEN 1
      ELSE 0
    END AS outcome_28d,
    CASE
      WHEN deathtime IS NOT NULL AND deathtime >= icu_intime THEN 'admissions.deathtime'
      WHEN dod IS NOT NULL AND dod >= icu_intime::date THEN 'patients.dod_date'
      ELSE 'none'
    END AS death_time_source
FROM base
)
SELECT
    *,
    GREATEST(EXTRACT(EPOCH FROM (
      CASE
        WHEN hospital_outcome = 1 AND deathtime IS NOT NULL
          THEN deathtime
        ELSE COALESCE(dischtime, icu_outtime)
      END - icu_intime))/86400.0, 0.001)
      AS hospital_time_days,
    GREATEST(EXTRACT(EPOCH FROM (
      CASE
        WHEN icu_outcome = 1 AND deathtime IS NOT NULL
          THEN deathtime
        ELSE icu_outtime
      END - icu_intime))/86400.0, 0.001)
      AS icu_time_days,
    CASE
      WHEN outcome_28d = 1 AND deathtime IS NOT NULL
      THEN GREATEST(EXTRACT(EPOCH FROM (deathtime - icu_intime))/86400.0, 0.001)
      WHEN outcome_28d = 1 AND dod IS NOT NULL
      THEN GREATEST((dod - icu_intime::date)::double precision, 0.001)
      ELSE 28.0
    END AS time_28d_days
FROM outcome_flags
ORDER BY subject_id, hadm_id, stay_id;

-- 受限原始事件文件：ICU入科前6 h至入科后24 h的四项暴露实验室事件。
DROP TABLE IF EXISTS pg_temp.nwicu_exposure_lab_events_export;
CREATE TEMP TABLE nwicu_exposure_lab_events_export AS
SELECT
    e.subject_id,
    e.hadm_id,
    e.stay_id,
    le.labevent_id AS lab_event_id,
    le.charttime,
    EXTRACT(EPOCH FROM (le.charttime - e.icu_intime))/3600.0 AS hours_from_icu_admission,
    CASE
      WHEN le.itemid IN (100001, 100045) THEN 'glucose'
      WHEN le.itemid = 100079 THEN 'hba1c'
      WHEN le.itemid = 100056 THEN 'triglycerides'
      WHEN le.itemid = 100086 THEN 'hdl'
    END AS analyte,
    le.itemid,
    le.valuenum,
    le.valueuom
FROM nwicu_external_validation_export e
INNER JOIN hosp.labevents le
  ON le.hadm_id = e.hadm_id
 AND le.charttime >= e.icu_intime - INTERVAL '6 hours'
 AND le.charttime <= e.icu_intime + INTERVAL '24 hours'
 AND le.itemid IN (100001, 100045, 100079, 100056, 100086)
WHERE le.valuenum IS NOT NULL AND le.valuenum > 0;

DROP TABLE IF EXISTS pg_temp.nwicu_exposure_lab_definitions_export;
CREATE TEMP TABLE nwicu_exposure_lab_definitions_export AS
WITH ranked AS (
  SELECT
    x.*,
    ROW_NUMBER() OVER (
      PARTITION BY stay_id, analyte ORDER BY charttime, lab_event_id
    ) AS rn_first,
    ROW_NUMBER() OVER (
      PARTITION BY stay_id, analyte
      ORDER BY ABS(hours_from_icu_admission), charttime, lab_event_id
    ) AS rn_closest
  FROM nwicu_exposure_lab_events_export x
)
SELECT
  subject_id,
  hadm_id,
  stay_id,
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
GROUP BY subject_id, hadm_id, stay_id, analyte;

-- 运行时通过psql变量传入受限输出路径。
-- 使用\o和COPY TO STDOUT，避免\copy不展开文件名变量的问题。
\o :nwicu_output_csv
COPY nwicu_external_validation_export TO STDOUT WITH (FORMAT CSV, HEADER true);
\o

\o :nwicu_lab_events_csv
COPY nwicu_exposure_lab_events_export TO STDOUT WITH (FORMAT CSV, HEADER true);
\o

\o :nwicu_lab_definitions_csv
COPY nwicu_exposure_lab_definitions_export TO STDOUT WITH (FORMAT CSV, HEADER true);
\o

SELECT
    COUNT(*) AS exported_rows,
    COUNT(*) FILTER (WHERE hospital_outcome = 1) AS hospital_deaths,
    COUNT(*) FILTER (WHERE icu_outcome = 1) AS icu_deaths,
    COUNT(*) FILTER (WHERE outcome_28d = 1) AS deaths_28d,
    COUNT(DISTINCT stay_id) AS unique_stays,
    COUNT(sofa) AS sofa_available,
    COUNT(*) FILTER (WHERE deathtime IS NULL AND hospital_outcome = 1) AS hospital_deaths_without_deathtime,
    COUNT(*) FILTER (WHERE charlson_comorbidity_index IS NULL) AS charlson_missing,
    COUNT(*) FILTER (WHERE glucose IS NOT NULL) AS glucose_available,
    COUNT(*) FILTER (WHERE hba1c IS NOT NULL) AS hba1c_available,
    COUNT(*) FILTER (WHERE triglycerides IS NOT NULL) AS triglycerides_available,
    COUNT(*) FILTER (WHERE hdl IS NOT NULL) AS hdl_available,
    COUNT(*) FILTER (
      WHERE glucose > 0 AND hba1c > 0 AND triglycerides > 0 AND hdl > 0
    ) AS all_three_complete
FROM nwicu_external_validation_export;
