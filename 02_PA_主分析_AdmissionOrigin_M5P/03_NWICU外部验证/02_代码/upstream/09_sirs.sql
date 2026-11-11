-- ====================================================================================================
-- sirs.sql (nwicu adaptation)
-- SIRS score computed on ICU Day 1 (0–24h from ICU admission)
-- Data sources: icu.chartevents (vitals), hosp.labevents (WBC, PaCO2)
--
-- chartevents itemids used:
--   320045 = PULSE (heart rate)
--   320210 = RESPIRATIONS (respiratory rate)
--   323761 = TEMPERATURE
-- labevents itemids:
--   100016 = White Blood Cells
--   100028 = Immature Granulocytes (bands %)
-- Blood Gas PaCO2 is not reliably available; resp_rate alone is used for respiratory criterion.
-- ====================================================================================================
DROP TABLE IF EXISTS derived.sirs;
CREATE TABLE derived.sirs AS
WITH vitals_d1 AS (
  -- Aggregate vital sign min/max over ICU Day 1 from chartevents
  SELECT
    ie.stay_id,
    MIN(CASE WHEN c.itemid = 323761 THEN c.valuenum END) AS temp_min,
    MAX(CASE WHEN c.itemid = 323761 THEN c.valuenum END) AS temp_max,
    MAX(CASE WHEN c.itemid IN (320045, 320050, 320051) THEN c.valuenum END) AS hr_max,
    MAX(CASE WHEN c.itemid = 320210 THEN c.valuenum END) AS rr_max
  FROM icu.icustays AS ie
  JOIN icu.chartevents AS c ON c.stay_id = ie.stay_id
  WHERE c.charttime BETWEEN ie.intime AND ie.intime + INTERVAL '24 hours'
    AND c.valuenum IS NOT NULL
    AND c.itemid IN (320045, 320050, 320051, 320210, 323761)
  GROUP BY ie.stay_id
),
labs_d1 AS (
  -- WBC and bands from labevents over ICU Day 1
  SELECT
    ie.stay_id,
    MIN(CASE WHEN l.itemid = 100016 THEN l.valuenum END) AS wbc_min,
    MAX(CASE WHEN l.itemid = 100016 THEN l.valuenum END) AS wbc_max,
    MAX(CASE WHEN l.itemid = 100028 THEN l.valuenum END) AS bands_max
  FROM icu.icustays AS ie
  JOIN hosp.labevents AS l ON l.subject_id = ie.subject_id
  WHERE l.charttime BETWEEN ie.intime AND ie.intime + INTERVAL '24 hours'
    AND l.valuenum IS NOT NULL
    AND l.itemid IN (100016, 100028)
  GROUP BY ie.stay_id
),
scorecomp AS (
  SELECT
    ie.stay_id,
    v.temp_min,
    v.temp_max,
    v.hr_max,
    v.rr_max,
    l.wbc_min,
    l.wbc_max,
    l.bands_max
  FROM icu.icustays AS ie
  LEFT JOIN vitals_d1 AS v ON ie.stay_id = v.stay_id
  LEFT JOIN labs_d1   AS l ON ie.stay_id = l.stay_id
),
scorecalc AS (
  SELECT
    stay_id,
    -- Temperature criterion: <36°C or >38°C
    CASE
      WHEN temp_min < 36.0 THEN 1
      WHEN temp_max > 38.0 THEN 1
      WHEN temp_min IS NULL THEN NULL
      ELSE 0
    END AS temp_score,
    -- Heart rate criterion: >90 bpm
    CASE
      WHEN hr_max > 90.0 THEN 1
      WHEN hr_max IS NULL THEN NULL
      ELSE 0
    END AS heart_rate_score,
    -- Respiratory criterion: RR >20
    CASE
      WHEN rr_max > 20.0 THEN 1
      WHEN rr_max IS NULL THEN NULL
      ELSE 0
    END AS resp_score,
    -- WBC criterion: <4 or >12 x10^9/L, or bands >10%
    CASE
      WHEN wbc_min < 4.0  THEN 1
      WHEN wbc_max > 12.0 THEN 1
      WHEN bands_max > 10 THEN 1
      WHEN COALESCE(wbc_min, bands_max) IS NULL THEN NULL
      ELSE 0
    END AS wbc_score
  FROM scorecomp
)
SELECT
  ie.subject_id,
  ie.hadm_id,
  ie.stay_id,
  COALESCE(temp_score, 0)
    + COALESCE(heart_rate_score, 0)
    + COALESCE(resp_score, 0)
    + COALESCE(wbc_score, 0)                    AS sirs,
  temp_score,
  heart_rate_score,
  resp_score,
  wbc_score
FROM icu.icustays AS ie
LEFT JOIN scorecalc AS s ON ie.stay_id = s.stay_id
ORDER BY ie.subject_id, ie.stay_id;
