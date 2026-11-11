\set ON_ERROR_STOP on

/*
NWICU first-day adapted SOFA (mSOFA), version 1.0

This script intentionally DOES NOT overwrite derived.sofa. The NWICU source
does not contain a directly usable PaO2 field, complete GCS observations,
true urine output, or inputevents-style administered vasopressor rates.

Adaptations:
  1. Respiratory: nearest SpO2/FiO2 ratio within 2 hours; FiO2 is itemid 100029.
  2. Cardiovascular: MAP from paired SBP/DBP plus EMAR-confirmed prescription rate.
  3. Renal: creatinine only; urine-output contribution is unavailable.
  4. Neurologic: uses derived.gcs when present; otherwise remains missing.
  5. Missing component scores are assumed normal only in sofa_adapted_24hours.

The table exposes provenance and missingness. It must not be labelled as a
fully standard SOFA or silently substituted into a cross-database main model.
*/

BEGIN;

CREATE SCHEMA IF NOT EXISTS derived;
DROP TABLE IF EXISTS derived.nwicu_sofa_day1_adapted;

CREATE TABLE derived.nwicu_sofa_day1_adapted AS
WITH cohort AS (
    SELECT
        i.subject_id,
        i.hadm_id,
        i.stay_id,
        i.intime,
        i.outtime,
        i.intime - INTERVAL '6 hours' AS lab_start,
        i.intime + INTERVAL '24 hours' AS day1_end
    FROM icu.icustays i
),
weight_ranked AS (
    SELECT
        c.stay_id,
        ce.valuenum AS weight_kg,
        ce.charttime,
        ROW_NUMBER() OVER (
            PARTITION BY c.stay_id
            ORDER BY ABS(EXTRACT(EPOCH FROM (ce.charttime-c.intime))), ce.charttime
        ) AS rn
    FROM cohort c
    JOIN icu.chartevents ce
      ON ce.stay_id=c.stay_id
     AND ce.itemid=326531
     AND ce.charttime>=c.lab_start
     AND ce.charttime<c.day1_end
     AND ce.valuenum BETWEEN 20 AND 350
),
weight AS (
    SELECT stay_id, weight_kg, charttime AS weight_charttime
    FROM weight_ranked
    WHERE rn=1
),
labs AS (
    SELECT
        c.stay_id,
        MIN(le.valuenum) FILTER (
            WHERE le.itemid=100014 AND le.valuenum BETWEEN 1 AND 2000
        ) AS platelet_min,
        MAX(le.valuenum) FILTER (
            WHERE le.itemid=100020 AND le.valuenum BETWEEN 0 AND 60
        ) AS bilirubin_max,
        MAX(le.valuenum) FILTER (
            WHERE le.itemid=100002 AND le.valuenum BETWEEN 0.1 AND 30
        ) AS creatinine_max
    FROM cohort c
    LEFT JOIN hosp.labevents le
      ON le.subject_id=c.subject_id
     AND le.hadm_id=c.hadm_id
     AND le.charttime>=c.lab_start
     AND le.charttime<c.day1_end
     AND le.itemid IN (100014,100020,100002)
     AND le.valuenum IS NOT NULL
     AND le.valuenum<>9999999
    GROUP BY c.stay_id
),
fio2_events AS (
    SELECT
        c.stay_id,
        le.charttime AS fio2_charttime,
        AVG(le.valuenum)/100.0 AS fio2_fraction,
        STRING_AGG(DISTINCT COALESCE(le.valueuom,'NULL'), '|') AS fio2_units
    FROM cohort c
    JOIN hosp.labevents le
      ON le.subject_id=c.subject_id
     AND le.hadm_id=c.hadm_id
     AND le.charttime>=c.intime
     AND le.charttime<c.day1_end
     AND le.itemid=100029
     AND le.valuenum BETWEEN 21 AND 100
     AND (le.valueuom IS NULL OR LOWER(le.valueuom) IN ('%','percent'))
    GROUP BY c.stay_id, le.charttime
),
resp_pairs_ranked AS (
    SELECT
        f.stay_id,
        f.fio2_charttime,
        f.fio2_fraction,
        f.fio2_units,
        s.charttime AS spo2_charttime,
        s.valuenum AS spo2,
        ABS(EXTRACT(EPOCH FROM (s.charttime-f.fio2_charttime)))/60.0 AS pair_gap_minutes,
        ROW_NUMBER() OVER (
            PARTITION BY f.stay_id, f.fio2_charttime
            ORDER BY ABS(EXTRACT(EPOCH FROM (s.charttime-f.fio2_charttime))), s.charttime
        ) AS rn
    FROM fio2_events f
    JOIN icu.chartevents s
      ON s.stay_id=f.stay_id
     AND s.itemid=320277
     AND s.charttime>=f.fio2_charttime-INTERVAL '2 hours'
     AND s.charttime<=f.fio2_charttime+INTERVAL '2 hours'
     AND s.valuenum BETWEEN 50 AND 100
),
resp_pairs AS (
    SELECT
        stay_id,
        fio2_charttime,
        fio2_fraction,
        fio2_units,
        spo2_charttime,
        spo2,
        pair_gap_minutes,
        spo2/fio2_fraction AS sf_ratio
    FROM resp_pairs_ranked
    WHERE rn=1
),
resp_worst AS (
    SELECT DISTINCT ON (stay_id)
        stay_id,
        sf_ratio AS sf_ratio_min,
        spo2 AS sf_spo2,
        fio2_fraction AS sf_fio2,
        fio2_charttime,
        spo2_charttime,
        pair_gap_minutes,
        fio2_units
    FROM resp_pairs
    WHERE sf_ratio>0
    ORDER BY stay_id, sf_ratio, pair_gap_minutes
),
resp_support AS (
    SELECT
        c.stay_id,
        BOOL_OR(pe.itemid=787541) AS invasive_vent_day1,
        BOOL_OR(pe.itemid=792843) AS noninvasive_vent_day1
    FROM cohort c
    LEFT JOIN icu.procedureevents pe
      ON pe.stay_id=c.stay_id
     AND pe.itemid IN (787541,792843)
     AND pe.starttime<c.day1_end
     AND COALESCE(pe.endtime,pe.starttime)>=c.intime
    GROUP BY c.stay_id
),
bp_paired AS (
    SELECT
        c.stay_id,
        ce.charttime,
        MAX(ce.valuenum) FILTER (WHERE ce.itemid=320050) AS art_sys,
        MAX(ce.valuenum) FILTER (WHERE ce.itemid=320051) AS art_dia,
        MAX(ce.valuenum) FILTER (WHERE ce.itemid=320179) AS nibp_sys,
        MAX(ce.valuenum) FILTER (WHERE ce.itemid=320180) AS nibp_dia
    FROM cohort c
    JOIN icu.chartevents ce
      ON ce.stay_id=c.stay_id
     AND ce.charttime>=c.intime
     AND ce.charttime<c.day1_end
     AND ce.itemid IN (320050,320051,320179,320180)
    GROUP BY c.stay_id, ce.charttime
),
map_values AS (
    SELECT
        b.stay_id,
        b.charttime,
        v.map_value,
        v.map_source
    FROM bp_paired b
    CROSS JOIN LATERAL (
        VALUES
          (CASE WHEN b.art_sys BETWEEN 20 AND 300
                      AND b.art_dia BETWEEN 10 AND 200
                      AND b.art_sys>=b.art_dia
                 THEN (b.art_sys+2.0*b.art_dia)/3.0 END, 'arterial'),
          (CASE WHEN b.nibp_sys BETWEEN 20 AND 300
                      AND b.nibp_dia BETWEEN 10 AND 200
                      AND b.nibp_sys>=b.nibp_dia
                 THEN (b.nibp_sys+2.0*b.nibp_dia)/3.0 END, 'noninvasive')
    ) AS v(map_value,map_source)
    WHERE v.map_value BETWEEN 20 AND 200
),
map_worst AS (
    SELECT DISTINCT ON (stay_id)
        stay_id,
        map_value AS map_min,
        charttime AS map_charttime,
        map_source
    FROM map_values
    ORDER BY stay_id, map_value, charttime
),
gcs_day1 AS (
    SELECT
        c.stay_id,
        MIN(g.gcs) FILTER (WHERE g.gcs BETWEEN 3 AND 15) AS gcs_min
    FROM cohort c
    LEFT JOIN derived.gcs g
      ON g.stay_id=c.stay_id
     AND g.charttime>=c.intime
     AND g.charttime<c.day1_end
    GROUP BY c.stay_id
),
emar_vaso_positive AS (
    SELECT DISTINCT
        c.stay_id,
        e.hadm_id,
        e.pharmacy_id
    FROM cohort c
    JOIN hosp.emar e
      ON e.hadm_id=c.hadm_id
     AND e.charttime>=c.intime
     AND e.charttime<c.day1_end
    WHERE e.pharmacy_id IS NOT NULL
      AND (
        LOWER(e.medication) LIKE 'norepinephrine%'
        OR LOWER(e.medication) LIKE 'epinephrine%'
        OR LOWER(e.medication) LIKE 'dopamine%'
        OR LOWER(e.medication) LIKE 'dobutamine%'
      )
      AND NULLIF(TRIM(e.event_txt),'') IS NOT NULL
      AND LOWER(COALESCE(e.event_txt,'')) NOT IN
          ('not given','cancelled','canceled','held','stopped')
),
vaso_orders_raw AS (
    SELECT
        c.stay_id,
        p.pharmacy_id,
        CASE
          WHEN LOWER(p.drug) LIKE 'norepinephrine%' THEN 'norepinephrine'
          WHEN LOWER(p.drug) LIKE 'epinephrine%' THEN 'epinephrine'
          WHEN LOWER(p.drug) LIKE 'dopamine%' THEN 'dopamine'
          WHEN LOWER(p.drug) LIKE 'dobutamine%' THEN 'dobutamine'
        END AS vaso_drug,
        p.starttime,
        p.stoptime,
        p.dose_unit_rx,
        CASE
          WHEN TRIM(COALESCE(p.dose_val_rx,'')) ~ '^[0-9]+([.][0-9]+)?$'
          THEN TRIM(p.dose_val_rx)::double precision
        END AS ordered_rate,
        w.weight_kg,
        (ev.pharmacy_id IS NOT NULL) AS emar_positive
    FROM cohort c
    JOIN hosp.prescriptions p
      ON p.subject_id=c.subject_id
     AND p.hadm_id=c.hadm_id
     AND p.starttime<c.day1_end
     AND COALESCE(p.stoptime,p.starttime+INTERVAL '1 hour')>=c.intime
     AND LOWER(COALESCE(p.route,''))='intravenous'
     AND (
        LOWER(p.drug) LIKE 'norepinephrine%'
        OR LOWER(p.drug) LIKE 'epinephrine%'
        OR LOWER(p.drug) LIKE 'dopamine%'
        OR LOWER(p.drug) LIKE 'dobutamine%'
     )
    LEFT JOIN weight w USING (stay_id)
    LEFT JOIN emar_vaso_positive ev
      ON ev.stay_id=c.stay_id
     AND ev.hadm_id=p.hadm_id
     AND ev.pharmacy_id=p.pharmacy_id
),
vaso_orders AS (
    SELECT
        *,
        CASE
          WHEN ordered_rate IS NULL OR NOT emar_positive THEN NULL
          WHEN LOWER(dose_unit_rx)='mcg/kg/min' THEN ordered_rate
          WHEN LOWER(dose_unit_rx)='mcg/min' AND weight_kg>0 THEN ordered_rate/weight_kg
        END AS rate_mcg_kg_min
    FROM vaso_orders_raw
),
vaso_agg AS (
    SELECT
        stay_id,
        MAX(rate_mcg_kg_min) FILTER (WHERE vaso_drug='norepinephrine') AS norepinephrine_max,
        MAX(rate_mcg_kg_min) FILTER (WHERE vaso_drug='epinephrine') AS epinephrine_max,
        MAX(rate_mcg_kg_min) FILTER (WHERE vaso_drug='dopamine') AS dopamine_max,
        MAX(rate_mcg_kg_min) FILTER (WHERE vaso_drug='dobutamine') AS dobutamine_max,
        BOOL_OR(vaso_drug='dobutamine' AND emar_positive) AS dobutamine_emar_positive,
        COUNT(*) FILTER (WHERE emar_positive) AS emar_positive_orders,
        COUNT(*) FILTER (
          WHERE emar_positive
            AND rate_mcg_kg_min IS NULL
            AND vaso_drug IN ('norepinephrine','epinephrine','dopamine')
        ) AS unquantified_positive_orders
    FROM vaso_orders
    GROUP BY stay_id
),
inputs AS (
    SELECT
        c.*,
        w.weight_kg,
        l.platelet_min,
        l.bilirubin_max,
        l.creatinine_max,
        rw.sf_ratio_min,
        rw.sf_spo2,
        rw.sf_fio2,
        rw.fio2_charttime,
        rw.spo2_charttime,
        rw.pair_gap_minutes,
        COALESCE(rs.invasive_vent_day1,false) AS invasive_vent_day1,
        COALESCE(rs.noninvasive_vent_day1,false) AS noninvasive_vent_day1,
        mw.map_min,
        mw.map_charttime,
        mw.map_source,
        gd.gcs_min,
        va.norepinephrine_max,
        va.epinephrine_max,
        va.dopamine_max,
        va.dobutamine_max,
        COALESCE(va.dobutamine_emar_positive,false) AS dobutamine_emar_positive,
        COALESCE(va.emar_positive_orders,0) AS emar_positive_orders,
        COALESCE(va.unquantified_positive_orders,0) AS unquantified_positive_orders
    FROM cohort c
    LEFT JOIN weight w USING (stay_id)
    LEFT JOIN labs l USING (stay_id)
    LEFT JOIN resp_worst rw USING (stay_id)
    LEFT JOIN resp_support rs USING (stay_id)
    LEFT JOIN map_worst mw USING (stay_id)
    LEFT JOIN gcs_day1 gd USING (stay_id)
    LEFT JOIN vaso_agg va USING (stay_id)
),
scores AS (
    SELECT
        i.*,
        CASE
          WHEN sf_ratio_min IS NULL THEN NULL
          WHEN (invasive_vent_day1 OR noninvasive_vent_day1) AND sf_ratio_min<100 THEN 4
          WHEN (invasive_vent_day1 OR noninvasive_vent_day1) AND sf_ratio_min<200 THEN 3
          WHEN sf_ratio_min<300 THEN 2
          WHEN sf_ratio_min<400 THEN 1
          ELSE 0
        END AS respiration_score,
        CASE
          WHEN platelet_min IS NULL THEN NULL
          WHEN platelet_min<20 THEN 4
          WHEN platelet_min<50 THEN 3
          WHEN platelet_min<100 THEN 2
          WHEN platelet_min<150 THEN 1
          ELSE 0
        END AS coagulation_score,
        CASE
          WHEN bilirubin_max IS NULL THEN NULL
          WHEN bilirubin_max>=12 THEN 4
          WHEN bilirubin_max>=6 THEN 3
          WHEN bilirubin_max>=2 THEN 2
          WHEN bilirubin_max>=1.2 THEN 1
          ELSE 0
        END AS liver_score,
        CASE
          WHEN dopamine_max>15 OR epinephrine_max>0.1 OR norepinephrine_max>0.1 THEN 4
          WHEN unquantified_positive_orders>0 THEN NULL
          WHEN dopamine_max>5
            OR (epinephrine_max>0 AND epinephrine_max<=0.1)
            OR (norepinephrine_max>0 AND norepinephrine_max<=0.1) THEN 3
          WHEN dopamine_max>0 OR dobutamine_max>0 OR dobutamine_emar_positive THEN 2
          WHEN map_min<70 THEN 1
          WHEN map_min IS NOT NULL THEN 0
          ELSE NULL
        END AS cardiovascular_score,
        CASE
          WHEN gcs_min IS NULL THEN NULL
          WHEN gcs_min<6 THEN 4
          WHEN gcs_min<10 THEN 3
          WHEN gcs_min<13 THEN 2
          WHEN gcs_min<15 THEN 1
          ELSE 0
        END AS neurologic_score,
        CASE
          WHEN creatinine_max IS NULL THEN NULL
          WHEN creatinine_max>=5 THEN 4
          WHEN creatinine_max>=3.5 THEN 3
          WHEN creatinine_max>=2 THEN 2
          WHEN creatinine_max>=1.2 THEN 1
          ELSE 0
        END AS renal_score
    FROM inputs i
)
SELECT
    subject_id,
    hadm_id,
    stay_id,
    intime,
    outtime,
    respiration_score,
    coagulation_score,
    liver_score,
    cardiovascular_score,
    neurologic_score,
    renal_score,
    COALESCE(respiration_score,0)
      +COALESCE(coagulation_score,0)
      +COALESCE(liver_score,0)
      +COALESCE(cardiovascular_score,0)
      +COALESCE(neurologic_score,0)
      +COALESCE(renal_score,0) AS sofa_adapted_24hours,
    COALESCE(respiration_score,0)
      +COALESCE(coagulation_score,0)
      +COALESCE(liver_score,0)
      +COALESCE(cardiovascular_score,0)
      +COALESCE(renal_score,0) AS sofa_non_neuro_24hours,
    NULL::integer AS sofa_standard_24hours,
    (respiration_score IS NOT NULL)::integer
      +(coagulation_score IS NOT NULL)::integer
      +(liver_score IS NOT NULL)::integer
      +(cardiovascular_score IS NOT NULL)::integer
      +(neurologic_score IS NOT NULL)::integer
      +(renal_score IS NOT NULL)::integer AS observed_component_count,
    CASE
      WHEN respiration_score IS NOT NULL
       AND coagulation_score IS NOT NULL
       AND liver_score IS NOT NULL
       AND cardiovascular_score IS NOT NULL
       AND renal_score IS NOT NULL
       AND neurologic_score IS NOT NULL THEN 'all_components_observed_but_adapted_sources'
      WHEN respiration_score IS NOT NULL
       AND coagulation_score IS NOT NULL
       AND liver_score IS NOT NULL
       AND cardiovascular_score IS NOT NULL
       AND renal_score IS NOT NULL THEN 'five_non_neurologic_components_observed'
      WHEN ((respiration_score IS NOT NULL)::integer
          +(coagulation_score IS NOT NULL)::integer
          +(liver_score IS NOT NULL)::integer
          +(cardiovascular_score IS NOT NULL)::integer
          +(neurologic_score IS NOT NULL)::integer
          +(renal_score IS NOT NULL)::integer)>=4 THEN 'four_or_more_components_observed'
      ELSE 'fewer_than_four_components_observed'
    END AS quality_tier,
    false AS is_standard_sofa,
    'NWICU_mSOFA_v1_SF_creatinine_only_EMAR_prescription_rate'::text AS algorithm_version,
    sf_ratio_min,
    sf_spo2,
    sf_fio2,
    fio2_charttime,
    spo2_charttime,
    pair_gap_minutes,
    invasive_vent_day1,
    noninvasive_vent_day1,
    platelet_min,
    bilirubin_max,
    map_min,
    map_charttime,
    map_source,
    gcs_min,
    creatinine_max,
    weight_kg,
    norepinephrine_max,
    epinephrine_max,
    dopamine_max,
    dobutamine_max,
    emar_positive_orders,
    unquantified_positive_orders,
    (respiration_score IS NULL) AS respiration_missing,
    (coagulation_score IS NULL) AS coagulation_missing,
    (liver_score IS NULL) AS liver_missing,
    (cardiovascular_score IS NULL) AS cardiovascular_missing,
    (neurologic_score IS NULL) AS neurologic_missing,
    (renal_score IS NULL) AS renal_missing,
    true AS respiratory_uses_sf_proxy,
    true AS renal_lacks_urine_output,
    true AS vasopressor_uses_emar_confirmed_prescription_rate,
    'Not valid as an unqualified standard SOFA; use only as explicitly labelled NWICU-adapted sensitivity covariate.'::text AS use_note
FROM scores;

ALTER TABLE derived.nwicu_sofa_day1_adapted
  ADD CONSTRAINT nwicu_sofa_day1_adapted_pk PRIMARY KEY (stay_id);

CREATE INDEX nwicu_sofa_day1_adapted_hadm_idx
  ON derived.nwicu_sofa_day1_adapted (subject_id,hadm_id);

COMMENT ON TABLE derived.nwicu_sofa_day1_adapted IS
'NWICU first-day adapted SOFA. Uses S/F respiratory proxy, creatinine-only renal score, and EMAR-confirmed prescription rates. Not a fully standard SOFA.';

COMMENT ON COLUMN derived.nwicu_sofa_day1_adapted.sofa_adapted_24hours IS
'Sum of six component columns with missing components assumed normal (0); inspect observed_component_count and quality_tier.';

COMMENT ON COLUMN derived.nwicu_sofa_day1_adapted.sofa_standard_24hours IS
'Always NULL because the available NWICU source cannot establish a fully standard SOFA.';

COMMIT;

-- Required post-build QA.
SELECT
    COUNT(*) AS rows,
    COUNT(DISTINCT stay_id) AS distinct_stays,
    MIN(sofa_adapted_24hours) AS min_score,
    MAX(sofa_adapted_24hours) AS max_score,
    COUNT(*) FILTER (WHERE observed_component_count>=5) AS observed_5plus,
    COUNT(*) FILTER (WHERE respiration_score IS NOT NULL) AS respiration_n,
    COUNT(*) FILTER (WHERE coagulation_score IS NOT NULL) AS coagulation_n,
    COUNT(*) FILTER (WHERE liver_score IS NOT NULL) AS liver_n,
    COUNT(*) FILTER (WHERE cardiovascular_score IS NOT NULL) AS cardiovascular_n,
    COUNT(*) FILTER (WHERE neurologic_score IS NOT NULL) AS neurologic_n,
    COUNT(*) FILTER (WHERE renal_score IS NOT NULL) AS renal_n,
    COUNT(sofa_standard_24hours) AS standard_sofa_n
FROM derived.nwicu_sofa_day1_adapted;

SELECT quality_tier, COUNT(*) AS rows
FROM derived.nwicu_sofa_day1_adapted
GROUP BY quality_tier
ORDER BY rows DESC;

SELECT sofa_adapted_24hours, COUNT(*) AS rows
FROM derived.nwicu_sofa_day1_adapted
GROUP BY sofa_adapted_24hours
ORDER BY sofa_adapted_24hours;

