\set ON_ERROR_STOP on
CREATE SCHEMA IF NOT EXISTS derived;
SET search_path TO derived, hosp, icu, public;

-- Build cohort_all_icu without sepsis filtration
DROP TABLE IF EXISTS derived.cohort_all_icu;
CREATE TABLE derived.cohort_all_icu AS
SELECT
  i.subject_id,
  i.hadm_id,
  i.stay_id,
  i.intime AS icu_intime,
  i.outtime AS icu_outtime,
  i.los AS icu_los_days,
  p.gender AS sex,
  p.anchor_age AS age
FROM icu.icustays i
INNER JOIN hosp.patients p
  ON i.subject_id = p.subject_id
INNER JOIN hosp.admissions a
  ON i.hadm_id = a.hadm_id
WHERE p.anchor_age >= 18
  AND i.los >= 1.0;

CREATE INDEX IF NOT EXISTS idx_nw_cohort_all_icu_stay
  ON derived.cohort_all_icu (stay_id);
CREATE INDEX IF NOT EXISTS idx_nw_cohort_all_icu_hadm
  ON derived.cohort_all_icu (hadm_id);
