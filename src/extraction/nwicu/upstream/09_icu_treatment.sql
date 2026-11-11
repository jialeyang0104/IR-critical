-- ====================================================================================================
-- icu_treatment.sql (nwicu)
-- ICU treatment flags for ALL ICU patients (one row per stay_id)
-- 
-- Data sources:
--   hosp.prescriptions → vasopressors, steroids
--   icu.procedureevents → mechanical ventilation, NIV, intubation, CRRT, HD, PD
-- ====================================================================================================
DROP TABLE IF EXISTS derived.icu_treatment;
CREATE TABLE derived.icu_treatment AS
WITH vasopressor AS (
  -- IV vasopressor/inotrope prescriptions that overlap the ICU stay
  SELECT DISTINCT ie.stay_id
  FROM icu.icustays AS ie
  JOIN hosp.prescriptions AS p
    ON p.subject_id = ie.subject_id
    AND p.starttime < ie.outtime
    AND p.stoptime  > ie.intime
  WHERE p.route ILIKE ANY(ARRAY['%IV%','%INTRAVENOUS%','%IVPB%'])
    AND LOWER(p.drug) SIMILAR TO
      '%(norepinephrine|dopamine|dobutamine|epinephrine|vasopressin|phenylephrine)%'
    -- Exclude topical / local / ophthalmic which contain epinephrine but are not vasopressors
    AND p.route NOT ILIKE ANY(ARRAY['%TOP%','%OPHT%','%NASAL%','%ORAL%','%PO%','%NG%','%INJ%'])
),
steroid AS (
  -- Systemic steroid prescriptions overlapping ICU stay
  SELECT DISTINCT ie.stay_id
  FROM icu.icustays AS ie
  JOIN hosp.prescriptions AS p
    ON p.subject_id = ie.subject_id
    AND p.starttime < ie.outtime
    AND p.stoptime  > ie.intime
  WHERE LOWER(p.drug) SIMILAR TO
      '%(hydrocortisone|methylprednisolone|dexamethasone|prednisolone|prednisone|fludrocortisone)%'
    -- Only systemic routes (not topical/ophthalmic)
    AND p.route NOT ILIKE ANY(ARRAY['%TOP%','%OPHT%','%OTIC%','%NASAL%','%CREAM%'])
),
vent AS (
  -- Mechanical ventilation (invasive) procedure events
  SELECT DISTINCT ie.stay_id
  FROM icu.icustays AS ie
  JOIN icu.procedureevents AS pe ON pe.stay_id = ie.stay_id
  WHERE pe.itemid = 787541  -- MECHANICAL VENTILATION
),
niv AS (
  -- Non-invasive ventilation
  SELECT DISTINCT ie.stay_id
  FROM icu.icustays AS ie
  JOIN icu.procedureevents AS pe ON pe.stay_id = ie.stay_id
  WHERE pe.itemid = 792843  -- NON-INVASIVE VENTILATION
),
intub AS (
  -- Intubation procedure
  SELECT DISTINCT ie.stay_id
  FROM icu.icustays AS ie
  JOIN icu.procedureevents AS pe ON pe.stay_id = ie.stay_id
  WHERE pe.itemid = 753461  -- INTUBATION
),
crrt AS (
  -- Continuous Renal Replacement Therapy
  SELECT DISTINCT ie.stay_id
  FROM icu.icustays AS ie
  JOIN icu.procedureevents AS pe ON pe.stay_id = ie.stay_id
  WHERE pe.itemid = 740477  -- CONTINUOUS RENAL REPLACEMENT THERAPY
),
hd AS (
  -- Intermittent hemodialysis
  SELECT DISTINCT ie.stay_id
  FROM icu.icustays AS ie
  JOIN icu.procedureevents AS pe ON pe.stay_id = ie.stay_id
  WHERE pe.itemid = 704890  -- HEMODIALYSIS
),
pd AS (
  -- Peritoneal dialysis
  SELECT DISTINCT ie.stay_id
  FROM icu.icustays AS ie
  JOIN icu.procedureevents AS pe ON pe.stay_id = ie.stay_id
  WHERE pe.itemid = 772042  -- PERITONEAL DIALYSIS
)
SELECT
  ie.subject_id,
  ie.hadm_id,
  ie.stay_id,
  ie.intime   AS icu_intime,
  ie.outtime  AS icu_outtime,
  ie.los      AS icu_los_days,
  -- Treatment flags (1 = used; 0 = not documented)
  CASE WHEN vp.stay_id IS NOT NULL THEN 1 ELSE 0 END AS vasopressor,
  CASE WHEN st.stay_id IS NOT NULL THEN 1 ELSE 0 END AS steroid,
  CASE WHEN vt.stay_id IS NOT NULL THEN 1 ELSE 0 END AS mechanical_ventilation,
  CASE WHEN ni.stay_id IS NOT NULL THEN 1 ELSE 0 END AS non_invasive_ventilation,
  CASE WHEN ib.stay_id IS NOT NULL THEN 1 ELSE 0 END AS intubation,
  CASE WHEN cr.stay_id IS NOT NULL THEN 1 ELSE 0 END AS crrt,
  CASE WHEN hd.stay_id IS NOT NULL THEN 1 ELSE 0 END AS hemodialysis,
  CASE WHEN pd.stay_id IS NOT NULL THEN 1 ELSE 0 END AS peritoneal_dialysis,
  -- Any renal replacement therapy
  CASE WHEN cr.stay_id IS NOT NULL OR hd.stay_id IS NOT NULL OR pd.stay_id IS NOT NULL
       THEN 1 ELSE 0 END AS any_rrt
FROM icu.icustays AS ie
LEFT JOIN vasopressor  AS vp ON ie.stay_id = vp.stay_id
LEFT JOIN steroid      AS st ON ie.stay_id = st.stay_id
LEFT JOIN vent         AS vt ON ie.stay_id = vt.stay_id
LEFT JOIN niv          AS ni ON ie.stay_id = ni.stay_id
LEFT JOIN intub        AS ib ON ie.stay_id = ib.stay_id
LEFT JOIN crrt         AS cr ON ie.stay_id = cr.stay_id
LEFT JOIN hd                 ON ie.stay_id = hd.stay_id
LEFT JOIN pd                 ON ie.stay_id = pd.stay_id
ORDER BY ie.subject_id, ie.stay_id;
