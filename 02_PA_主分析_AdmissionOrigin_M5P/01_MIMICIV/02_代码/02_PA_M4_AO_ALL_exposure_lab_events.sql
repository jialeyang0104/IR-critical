/* Restricted event-level source for alternative laboratory-timing analyses.
   This query creates no persistent database objects. */
WITH patient_order AS (
  SELECT stay_id, ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY intime, stay_id) AS rn
  FROM mimiciv_icu.icustays
), eligible AS (
  SELECT d.hadm_id, d.stay_id, d.icu_intime
  FROM mimiciv_derived.icustay_detail AS d
  JOIN patient_order r ON r.stay_id=d.stay_id
  WHERE r.rn=1
    AND d.admission_age >= 18
    AND d.icu_outtime >= d.icu_intime + INTERVAL '24 hours'
    AND d.dischtime >= d.icu_intime + INTERVAL '24 hours'
)
SELECT
  e.stay_id,
  CASE
    WHEN le.itemid IN (50809,50931,52027,52569) THEN 'glucose'
    WHEN le.itemid IN (50852,51631) THEN 'hba1c'
    WHEN le.itemid = 51000 THEN 'tg'
    WHEN le.itemid = 50904 THEN 'hdl'
  END AS analyte,
  le.valuenum,
  le.valueuom,
  le.charttime,
  EXTRACT(EPOCH FROM (le.charttime - e.icu_intime)) / 3600.0 AS hours_from_icu,
  ABS(EXTRACT(EPOCH FROM (le.charttime - e.icu_intime)) / 3600.0) AS abs_hours_from_icu
FROM eligible AS e
CROSS JOIN LATERAL (
  SELECT x.itemid, x.valuenum, x.valueuom, x.charttime, x.labevent_id
  FROM mimiciv_hosp.labevents AS x
  WHERE x.hadm_id = e.hadm_id
    AND x.charttime >= e.icu_intime - INTERVAL '6 hours'
    AND x.charttime <= e.icu_intime + INTERVAL '24 hours'
    AND x.itemid IN (50809,50931,52027,52569,50852,51631,51000,50904)
    AND x.valuenum IS NOT NULL
  OFFSET 0
) AS le
ORDER BY e.stay_id, analyte, le.charttime, le.labevent_id;
