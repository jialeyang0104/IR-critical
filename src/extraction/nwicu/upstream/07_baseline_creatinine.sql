-- ====================================================================================================
-- creatinine_baseline.sql (nwicu adaptation)
-- Estimates baseline serum creatinine for each hospital admission.
-- 
-- Logic (from MIMIC-IV):
--   i.  If lowest creatinine during admission <= 1.1 → use that value
--   ii. If patient has CKD (N18.x ICD-10) → use lowest value during admission
--   iii. Otherwise → estimate using Simplified MDRD:
--          eGFR = 186 × SCr^(-1.154) × Age^(-0.203) × [0.742 if female]
--          rearranged: SCr = (75 / 186 / age^(-0.203) [/ 0.742 if female])^(1/-1.154)
-- ====================================================================================================
DROP TABLE IF EXISTS derived.creatinine_baseline;
CREATE TABLE derived.creatinine_baseline AS
WITH p AS (
  SELECT
    ag.subject_id,
    ag.hadm_id,
    ag.age,
    pt.gender,
    -- Simplified MDRD estimated baseline SCr assuming eGFR = 75 mL/min/1.73m²
    CASE
      WHEN pt.gender = 'F'
        THEN POWER(75.0 / 186.0 / POWER(ag.age, -0.203) / 0.742, 1.0 / -1.154)
      ELSE
        POWER(75.0 / 186.0 / POWER(ag.age, -0.203), 1.0 / -1.154)
    END AS mdrd_est
  FROM derived.age AS ag
  LEFT JOIN hosp.patients AS pt ON ag.subject_id = pt.subject_id
  WHERE ag.age >= 18
),
lab AS (
  -- Lowest creatinine during the entire hospital admission
  SELECT
    ad.hadm_id,
    MIN(l.valuenum) AS scr_min
  FROM hosp.admissions AS ad
  JOIN hosp.labevents AS l
    ON l.subject_id = ad.subject_id
    AND l.charttime BETWEEN ad.admittime AND ad.dischtime
  WHERE l.itemid = 100002  -- Creatinine (Blood, Chemistry)
    AND l.valuenum IS NOT NULL
    AND l.valuenum > 0
  GROUP BY ad.hadm_id
),
ckd AS (
  -- CKD flag based on ICD-10 N18.x codes
  SELECT
    hadm_id,
    1 AS ckd_flag
  FROM hosp.diagnoses_icd
  WHERE icd_version = 10
    AND SUBSTR(icd_code, 1, 3) = 'N18'
  GROUP BY hadm_id
)
SELECT
  p.subject_id,
  p.hadm_id,
  p.gender,
  p.age,
  lab.scr_min,
  COALESCE(ckd.ckd_flag, 0)                        AS ckd,
  ROUND(p.mdrd_est::numeric, 3)                     AS mdrd_est,
  -- Final baseline determination
  CASE
    WHEN lab.scr_min <= 1.1    THEN lab.scr_min
    WHEN ckd.ckd_flag = 1      THEN lab.scr_min
    ELSE ROUND(p.mdrd_est::numeric, 3)
  END                                               AS scr_baseline
FROM p
LEFT JOIN lab ON p.hadm_id = lab.hadm_id
LEFT JOIN ckd ON p.hadm_id = ckd.hadm_id
ORDER BY p.subject_id, p.hadm_id;
