-- ====================================================================================================
-- kdigo_stages.sql (nwicu adaptation — creatinine-only)
-- AKI staging per KDIGO 2012 guidelines, based on serum creatinine only.
-- NOTE: Urine output staging is omitted because outputevents does not exist in nwicu.
--
-- Creatinine KDIGO stages:
--   Stage 1: SCr >= baseline × 1.5  OR  SCr increase >= 0.3 within 48h
--   Stage 2: SCr >= baseline × 2.0
--   Stage 3: SCr >= baseline × 3.0  OR  SCr >= 4.0 with acute rise >= 0.3
-- ====================================================================================================
DROP TABLE IF EXISTS derived.kdigo_stages;
CREATE TABLE derived.kdigo_stages AS
WITH cr_raw AS (
  -- All creatinine measurements per ICU stay
  SELECT
    ie.stay_id,
    ie.subject_id,
    ie.hadm_id,
    l.charttime,
    l.valuenum AS creat
  FROM icu.icustays AS ie
  JOIN hosp.labevents AS l
    ON l.subject_id = ie.subject_id
    AND l.charttime BETWEEN (ie.intime - INTERVAL '7 days') AND ie.outtime
  WHERE l.itemid = 100002  -- Creatinine
    AND l.valuenum > 0
    AND l.valuenum IS NOT NULL
),
cr_with_baseline AS (
  SELECT
    cr.stay_id,
    cr.subject_id,
    cr.hadm_id,
    cr.charttime,
    cr.creat,
    -- Rolling minimum over the past 7 days
    MIN(cr.creat) OVER (
      PARTITION BY cr.stay_id
      ORDER BY cr.charttime
      RANGE BETWEEN INTERVAL '7 days' PRECEDING AND CURRENT ROW
    ) AS creat_low_past_7day,
    -- Rolling minimum over the past 48 hours
    MIN(cr.creat) OVER (
      PARTITION BY cr.stay_id
      ORDER BY cr.charttime
      RANGE BETWEEN INTERVAL '48 hours' PRECEDING AND CURRENT ROW
    ) AS creat_low_past_48hr
  FROM cr_raw AS cr
),
cr_stg AS (
  SELECT
    stay_id,
    subject_id,
    hadm_id,
    charttime,
    creat,
    creat_low_past_7day,
    creat_low_past_48hr,
    CASE
      -- KDIGO Stage 3: SCr >= 3x baseline
      WHEN creat >= creat_low_past_7day * 3.0 THEN 3
      -- KDIGO Stage 3: SCr >= 4.0 with acute rise
      WHEN creat >= 4.0
       AND (creat_low_past_48hr <= 3.7 OR creat >= creat_low_past_7day * 1.5) THEN 3
      -- KDIGO Stage 2
      WHEN creat >= creat_low_past_7day * 2.0 THEN 2
      -- KDIGO Stage 1: acute 0.3 rise within 48h
      WHEN creat >= creat_low_past_48hr + 0.3 THEN 1
      -- KDIGO Stage 1: 1.5x baseline
      WHEN creat >= creat_low_past_7day * 1.5 THEN 1
      ELSE 0
    END AS aki_stage_creat
  FROM cr_with_baseline
)
SELECT
  cr.stay_id,
  cr.subject_id,
  cr.hadm_id,
  cr.charttime,
  cr.creat,
  cr.creat_low_past_7day,
  cr.creat_low_past_48hr,
  cr.aki_stage_creat,
  -- Smoothed: carry forward highest KDIGO stage from the past 6 hours
  MAX(cr.aki_stage_creat) OVER (
    PARTITION BY cr.stay_id
    ORDER BY cr.charttime
    RANGE BETWEEN INTERVAL '6 hours' PRECEDING AND CURRENT ROW
  ) AS aki_stage_smoothed
FROM cr_stg AS cr
ORDER BY cr.stay_id, cr.charttime;
