-- ====================================================================================================
-- icu_fluid_daily.sql (nwicu)
-- Daily fluid balance approximation for all ICU patients.
--
-- LIMITATIONS (important — document in any publication):
--   - INPUT:  Approximated from hosp.prescriptions IV fluid orders.
--             Records whether an IV fluid order was active each ICU day.
--             Actual infused volumes are NOT recorded (no inputevents table).
--   - OUTPUT: Urine output and other fluid losses are NOT available (no outputevents).
--             fluid_out and net_balance columns will always be NULL.
--
-- Day definitions (relative to icu_intime):
--   d1: [−6h, +24h)   d2: [24h, 48h)  ... d7: [144h, 168h)
-- ====================================================================================================
DROP TABLE IF EXISTS derived.icu_fluid_daily;
CREATE TABLE derived.icu_fluid_daily AS
WITH day_grid AS (
  -- Generate one row per patient per ICU day (d1–d7)
  SELECT
    ie.subject_id,
    ie.hadm_id,
    ie.stay_id,
    ie.intime AS icu_intime,
    d.icu_day,
    ie.intime + (d.day_start * INTERVAL '1 hour') AS window_start,
    ie.intime + (d.day_end   * INTERVAL '1 hour') AS window_end
  FROM icu.icustays AS ie
  CROSS JOIN (VALUES
    ('d1',  -6,  24),
    ('d2',  24,  48),
    ('d3',  48,  72),
    ('d4',  72,  96),
    ('d5',  96, 120),
    ('d6', 120, 144),
    ('d7', 144, 168)
  ) AS d(icu_day, day_start, day_end)
  -- Only include days that overlap with actual ICU stay
  WHERE ie.intime + (d.day_start * INTERVAL '1 hour') < ie.outtime
),
iv_fluid AS (
  -- IV fluid orders (crystalloids/colloids) active in each day window
  -- classified from prescriptions: route is IV or IVPB, drug is a common fluid
  SELECT
    dg.stay_id,
    dg.icu_day,
    COUNT(DISTINCT p.pharmacy_id) AS iv_fluid_order_count
  FROM day_grid AS dg
  JOIN hosp.prescriptions AS p
    ON p.subject_id = dg.subject_id
    AND p.starttime < dg.window_end
    AND p.stoptime  > dg.window_start
  WHERE p.route ILIKE ANY(ARRAY['%IV%','%IVPB%','%INTRAVENOUS%'])
    AND LOWER(p.drug) SIMILAR TO
      '%(sodium chloride|normal saline|lactated ringer|dextrose|albumin|hetastarch|gelatin|plasmalyte|ringers|fluid|saline)%'
  GROUP BY dg.stay_id, dg.icu_day
),
vasopressor_daily AS (
  -- Track vasopressor use per day
  SELECT
    dg.stay_id,
    dg.icu_day,
    COUNT(DISTINCT p.pharmacy_id) AS vasopressor_order_count
  FROM day_grid AS dg
  JOIN hosp.prescriptions AS p
    ON p.subject_id = dg.subject_id
    AND p.starttime < dg.window_end
    AND p.stoptime  > dg.window_start
  WHERE p.route ILIKE ANY(ARRAY['%IV%','%INTRAVENOUS%'])
    AND LOWER(p.drug) SIMILAR TO
      '%(norepinephrine|dopamine|dobutamine|epinephrine|vasopressin|phenylephrine)%'
    AND p.route NOT ILIKE ANY(ARRAY['%TOP%','%OPHT%','%NASAL%','%ORAL%','%PO%'])
  GROUP BY dg.stay_id, dg.icu_day
)
SELECT
  dg.subject_id,
  dg.hadm_id,
  dg.stay_id,
  dg.icu_intime,
  dg.icu_day,
  dg.window_start,
  dg.window_end,
  -- IV fluid orders active during this day
  COALESCE(fl.iv_fluid_order_count, 0)          AS iv_fluid_order_count,
  CASE WHEN fl.iv_fluid_order_count > 0
       THEN 1 ELSE 0 END                        AS received_iv_fluid,
  -- Vasopressor orders active during this day
  COALESCE(vp.vasopressor_order_count, 0)       AS vasopressor_order_count,
  CASE WHEN vp.vasopressor_order_count > 0
       THEN 1 ELSE 0 END                        AS received_vasopressor,
  -- Output not available
  NULL::double precision                         AS fluid_out_ml,
  NULL::double precision                         AS net_fluid_balance_ml
FROM day_grid AS dg
LEFT JOIN iv_fluid        AS fl ON dg.stay_id = fl.stay_id AND dg.icu_day = fl.icu_day
LEFT JOIN vasopressor_daily AS vp ON dg.stay_id = vp.stay_id AND dg.icu_day = vp.icu_day
ORDER BY dg.subject_id, dg.stay_id, dg.icu_day;
