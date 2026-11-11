\set ON_ERROR_STOP on
CREATE SCHEMA IF NOT EXISTS derived;
SET search_path TO derived, hosp, icu, public;

/*
NWICU final extraction (aligned to zigong s3 part2 format/rules).

Dependencies (run before this script):
  - derived.cohort_all_icu
  - derived.all_icu_labs_d1
  - derived.all_icu_labs_daily
  - derived.charlson
  - derived.sirs
  - derived.gcs
  - derived.creatinine_baseline
  - derived.kdigo_stages
  - derived.icu_treatment
  - derived.icu_fluid_daily

Output:
  - derived.cohort_all_icu_lab_final_select

Notes:
  - Fluid output/net values are aggregated from derived.icu_fluid_daily when present.
  - KDIGO UO stage remains NULL when urine-output data are unavailable in source.
*/

DROP TABLE IF EXISTS derived.cohort_all_icu_lab_final_select;
CREATE TABLE derived.cohort_all_icu_lab_final_select AS
WITH labs_extra_d1 AS (
    SELECT
        c.subject_id,
        c.hadm_id,
        c.stay_id,
        AVG(CASE WHEN l.itemid IN (100027, 100063) THEN l.valuenum END) AS total_protein,
        AVG(
            CASE
                WHEN lower(COALESCE(dl.label, '')) LIKE '%pco2%'
                  OR lower(COALESCE(dl.label, '')) LIKE '%partial pressure of carbon dioxide%'
                    THEN l.valuenum
                ELSE NULL
            END
        ) AS bg_pco2
    FROM derived.cohort_all_icu c
    JOIN hosp.labevents l
      ON l.subject_id = c.subject_id
     AND l.hadm_id = c.hadm_id
     AND l.charttime >= c.icu_intime - INTERVAL '6 hours'
     AND l.charttime <  c.icu_intime + INTERVAL '24 hours'
     AND l.valuenum IS NOT NULL
    LEFT JOIN hosp.d_labitems dl
      ON l.itemid = dl.itemid
    WHERE l.itemid IN (100027, 100063)
       OR lower(COALESCE(dl.label, '')) LIKE '%pco2%'
       OR lower(COALESCE(dl.label, '')) LIKE '%partial pressure of carbon dioxide%'
    GROUP BY c.subject_id, c.hadm_id, c.stay_id
),
final_agg AS (
    SELECT
        b.subject_id AS agg_subject_id,
        b.hadm_id AS agg_hadm_id,
        b.stay_id AS agg_stay_id,
        d.wbc AS "fd_WBC_mean",
        d.rbc AS "fd_RBC_mean",
        d.hemoglobin AS "fd_HGB_mean",
        d.hematocrit AS "fd_HCT_mean",
        d.platelet AS "fd_PLT_mean",
        d.neutrophils_pct AS "fd_NEU_PCT_mean",
        d.abs_neutrophil AS "fd_NEU_COUNT_mean",
        d.lymphocytes_pct AS "fd_LYM_PCT_mean",
        d.abs_lymphocyte AS "fd_LYM_COUNT_mean",
        d.monocytes_pct AS "fd_MON_PCT_mean",
        d.abs_monocyte AS "fd_MON_COUNT_mean",
        d.eosinophils_pct AS "fd_EOS_PCT_mean",
        d.abs_eosinophil AS "fd_EOS_COUNT_mean",
        d.basophils_pct AS "fd_BAS_PCT_mean",
        d.abs_basophil AS "fd_BAS_COUNT_mean",
        d.rdw AS "fd_RDW_mean",
        d.mcv AS "fd_MCV_mean",
        d.mch AS "fd_MCH_mean",
        d.mchc AS "fd_MCHC_mean",
        d.glucose AS "fd_GLU_mean",
        d.hba1c AS "fd_HBA1C_mean",
        d.cholesterol_total AS "fd_TC_mean",
        d.triglycerides AS "fd_TG_mean",
        d.hdl AS "fd_HDL_mean",
        d.ldl AS "fd_LDL_mean",
        d.alt AS "fd_ALT_mean",
        d.ast AS "fd_AST_mean",
        d.alp AS "fd_ALP_mean",
        d.ggt AS "fd_GGT_mean",
        d.bilirubin_total AS "fd_TBIL_mean",
        d.bilirubin_direct AS "fd_DBIL_mean",
        d.albumin AS "fd_ALB_mean",
        le.total_protein AS "fd_TP_mean",
        d.creatinine AS "fd_CRE_mean",
        d.bun AS "fd_UREA_mean",
        cb.mdrd_est::double precision AS "fd_EGFR_mean",
        d.sodium AS "fd_NA_mean",
        d.potassium AS "fd_K_mean",
        d.chloride AS "fd_CL_mean",
        d.calcium_total AS "fd_CA_mean",
        d.magnesium AS "fd_MG_mean",
        d.phosphate AS "fd_P_mean",
        d.total_co2 AS "fd_CO2_mean",
        d.anion_gap AS "fd_ANIONGAP_mean",
        d.crp AS "fd_CRP_mean",
        d.ferritin AS "fd_FERRITIN_mean",
        d.ldh AS "fd_LDH_mean",
        d.d_dimer AS "fd_DDIMER_mean",
        d.pt AS "fd_PT_mean",
        d.inr AS "fd_INR_mean",
        d.ptt AS "fd_APTT_mean",
        d.fibrinogen AS "fd_FIB_mean",
        d.bg_ph AS "fd_PH_mean",
        d.bg_po2 AS "fd_PO2_mean",
        le.bg_pco2 AS "fd_PCO2_mean",
        d.bg_o2sat AS "fd_SO2_mean",
        d.bg_lactate AS "fd_LAC_mean",
        sw.ldh_m6_24h_mean, sw.ldh_24_48h_mean, sw.ldh_48_72h_mean, sw.ldh_72_96h_mean, sw.ldh_96_120h_mean, sw.ldh_120_144h_mean, sw.ldh_144_168h_mean,
        sw.creatinine_m6_24h_mean, sw.creatinine_24_48h_mean, sw.creatinine_48_72h_mean, sw.creatinine_72_96h_mean, sw.creatinine_96_120h_mean, sw.creatinine_120_144h_mean, sw.creatinine_144_168h_mean,
        sw.plt_m6_24h_mean, sw.plt_24_48h_mean, sw.plt_48_72h_mean, sw.plt_72_96h_mean, sw.plt_96_120h_mean, sw.plt_120_144h_mean, sw.plt_144_168h_mean,
        sw.mono_count_m6_24h_mean, sw.mono_count_24_48h_mean, sw.mono_count_48_72h_mean, sw.mono_count_72_96h_mean, sw.mono_count_96_120h_mean, sw.mono_count_120_144h_mean, sw.mono_count_144_168h_mean,
        sw.neu_count_m6_24h_mean, sw.neu_count_24_48h_mean, sw.neu_count_48_72h_mean, sw.neu_count_72_96h_mean, sw.neu_count_96_120h_mean, sw.neu_count_120_144h_mean, sw.neu_count_144_168h_mean,
        sw.lym_count_m6_24h_mean, sw.lym_count_24_48h_mean, sw.lym_count_48_72h_mean, sw.lym_count_72_96h_mean, sw.lym_count_96_120h_mean, sw.lym_count_120_144h_mean, sw.lym_count_144_168h_mean
    FROM derived.cohort_all_icu b
    LEFT JOIN derived.all_icu_labs_d1 d
      ON b.subject_id = d.subject_id
     AND b.hadm_id = d.hadm_id
     AND b.stay_id = d.stay_id
    LEFT JOIN derived.creatinine_baseline cb
      ON b.hadm_id = cb.hadm_id
    LEFT JOIN labs_extra_d1 le
      ON b.subject_id = le.subject_id
     AND b.hadm_id = le.hadm_id
     AND b.stay_id = le.stay_id
    LEFT JOIN (
        SELECT
            t.subject_id,
            t.hadm_id,
            t.stay_id,
            AVG(CASE WHEN t.icu_day = 'd1' THEN t.ldh END) AS ldh_m6_24h_mean,
            AVG(CASE WHEN t.icu_day = 'd2' THEN t.ldh END) AS ldh_24_48h_mean,
            AVG(CASE WHEN t.icu_day = 'd3' THEN t.ldh END) AS ldh_48_72h_mean,
            AVG(CASE WHEN t.icu_day = 'd4' THEN t.ldh END) AS ldh_72_96h_mean,
            AVG(CASE WHEN t.icu_day = 'd5' THEN t.ldh END) AS ldh_96_120h_mean,
            AVG(CASE WHEN t.icu_day = 'd6' THEN t.ldh END) AS ldh_120_144h_mean,
            AVG(CASE WHEN t.icu_day = 'd7' THEN t.ldh END) AS ldh_144_168h_mean,
            AVG(CASE WHEN t.icu_day = 'd1' THEN t.creatinine END) AS creatinine_m6_24h_mean,
            AVG(CASE WHEN t.icu_day = 'd2' THEN t.creatinine END) AS creatinine_24_48h_mean,
            AVG(CASE WHEN t.icu_day = 'd3' THEN t.creatinine END) AS creatinine_48_72h_mean,
            AVG(CASE WHEN t.icu_day = 'd4' THEN t.creatinine END) AS creatinine_72_96h_mean,
            AVG(CASE WHEN t.icu_day = 'd5' THEN t.creatinine END) AS creatinine_96_120h_mean,
            AVG(CASE WHEN t.icu_day = 'd6' THEN t.creatinine END) AS creatinine_120_144h_mean,
            AVG(CASE WHEN t.icu_day = 'd7' THEN t.creatinine END) AS creatinine_144_168h_mean,
            AVG(CASE WHEN t.icu_day = 'd1' THEN t.platelet END) AS plt_m6_24h_mean,
            AVG(CASE WHEN t.icu_day = 'd2' THEN t.platelet END) AS plt_24_48h_mean,
            AVG(CASE WHEN t.icu_day = 'd3' THEN t.platelet END) AS plt_48_72h_mean,
            AVG(CASE WHEN t.icu_day = 'd4' THEN t.platelet END) AS plt_72_96h_mean,
            AVG(CASE WHEN t.icu_day = 'd5' THEN t.platelet END) AS plt_96_120h_mean,
            AVG(CASE WHEN t.icu_day = 'd6' THEN t.platelet END) AS plt_120_144h_mean,
            AVG(CASE WHEN t.icu_day = 'd7' THEN t.platelet END) AS plt_144_168h_mean,
            AVG(CASE WHEN t.icu_day = 'd1' THEN t.abs_monocyte END) AS mono_count_m6_24h_mean,
            AVG(CASE WHEN t.icu_day = 'd2' THEN t.abs_monocyte END) AS mono_count_24_48h_mean,
            AVG(CASE WHEN t.icu_day = 'd3' THEN t.abs_monocyte END) AS mono_count_48_72h_mean,
            AVG(CASE WHEN t.icu_day = 'd4' THEN t.abs_monocyte END) AS mono_count_72_96h_mean,
            AVG(CASE WHEN t.icu_day = 'd5' THEN t.abs_monocyte END) AS mono_count_96_120h_mean,
            AVG(CASE WHEN t.icu_day = 'd6' THEN t.abs_monocyte END) AS mono_count_120_144h_mean,
            AVG(CASE WHEN t.icu_day = 'd7' THEN t.abs_monocyte END) AS mono_count_144_168h_mean,
            AVG(CASE WHEN t.icu_day = 'd1' THEN t.abs_neutrophil END) AS neu_count_m6_24h_mean,
            AVG(CASE WHEN t.icu_day = 'd2' THEN t.abs_neutrophil END) AS neu_count_24_48h_mean,
            AVG(CASE WHEN t.icu_day = 'd3' THEN t.abs_neutrophil END) AS neu_count_48_72h_mean,
            AVG(CASE WHEN t.icu_day = 'd4' THEN t.abs_neutrophil END) AS neu_count_72_96h_mean,
            AVG(CASE WHEN t.icu_day = 'd5' THEN t.abs_neutrophil END) AS neu_count_96_120h_mean,
            AVG(CASE WHEN t.icu_day = 'd6' THEN t.abs_neutrophil END) AS neu_count_120_144h_mean,
            AVG(CASE WHEN t.icu_day = 'd7' THEN t.abs_neutrophil END) AS neu_count_144_168h_mean,
            AVG(CASE WHEN t.icu_day = 'd1' THEN t.abs_lymphocyte END) AS lym_count_m6_24h_mean,
            AVG(CASE WHEN t.icu_day = 'd2' THEN t.abs_lymphocyte END) AS lym_count_24_48h_mean,
            AVG(CASE WHEN t.icu_day = 'd3' THEN t.abs_lymphocyte END) AS lym_count_48_72h_mean,
            AVG(CASE WHEN t.icu_day = 'd4' THEN t.abs_lymphocyte END) AS lym_count_72_96h_mean,
            AVG(CASE WHEN t.icu_day = 'd5' THEN t.abs_lymphocyte END) AS lym_count_96_120h_mean,
            AVG(CASE WHEN t.icu_day = 'd6' THEN t.abs_lymphocyte END) AS lym_count_120_144h_mean,
            AVG(CASE WHEN t.icu_day = 'd7' THEN t.abs_lymphocyte END) AS lym_count_144_168h_mean
        FROM derived.all_icu_labs_daily t
        GROUP BY t.subject_id, t.hadm_id, t.stay_id
    ) sw
      ON b.subject_id = sw.subject_id
     AND b.hadm_id = sw.hadm_id
     AND b.stay_id = sw.stay_id
),
day1_flags AS (
    SELECT
        ie.stay_id,
        MAX(
            CASE
                WHEN LOWER(p.drug) SIMILAR TO '%(hydrocortisone|methylprednisolone|dexamethasone|prednisolone|prednisone|fludrocortisone)%'
                  THEN 1 ELSE 0
            END
        ) AS steroid_day1,
        MAX(
            CASE
                WHEN p.route ILIKE ANY(ARRAY['%IV%','%INTRAVENOUS%','%IVPB%'])
                 AND LOWER(p.drug) SIMILAR TO '%(norepinephrine|dopamine|dobutamine|epinephrine|vasopressin|phenylephrine)%'
                  THEN 1 ELSE 0
            END
        ) AS vasopressor_day1,
        MAX(
            CASE
                WHEN pe.itemid IN (787541, 792843, 753461) THEN 1 ELSE 0
            END
        ) AS ventilation_day1,
        MAX(
            CASE
                WHEN pe.itemid IN (740477, 704890, 772042) THEN 1 ELSE 0
            END
        ) AS crrt_day1
    FROM icu.icustays ie
    LEFT JOIN hosp.prescriptions p
      ON p.subject_id = ie.subject_id
     AND p.starttime < ie.intime + INTERVAL '24 hours'
     AND p.stoptime > ie.intime - INTERVAL '6 hours'
    LEFT JOIN icu.procedureevents pe
      ON pe.stay_id = ie.stay_id
     AND pe.starttime < ie.intime + INTERVAL '24 hours'
     AND COALESCE(pe.endtime, pe.starttime + INTERVAL '1 hour') > ie.intime - INTERVAL '6 hours'
    GROUP BY ie.stay_id
),
icu_treat_agg AS (
    SELECT
        t.subject_id,
        t.hadm_id,
        t.stay_id,
        1::int AS itd_icu_episode_count,
        GREATEST(1, CEIL(COALESCE(t.icu_los_days, EXTRACT(EPOCH FROM (t.icu_outtime - t.icu_intime)) / 86400.0)))::int AS itd_icu_days_total,
        EXTRACT(EPOCH FROM (t.icu_intime - ad.admittime)) / 3600.0 AS itd_first_icu_intime,
        EXTRACT(EPOCH FROM (t.icu_outtime - ad.admittime)) / 3600.0 AS itd_last_icu_outtime,
        COALESCE(t.steroid, 0) AS itd_used_steroid_icu,
        GREATEST(COALESCE(t.mechanical_ventilation, 0), COALESCE(t.non_invasive_ventilation, 0), COALESCE(t.intubation, 0)) AS itd_used_ventilation_icu,
        COALESCE(t.any_rrt, 0) AS itd_used_crrt_icu,
        COALESCE(t.vasopressor, 0) AS itd_used_vasopressor_icu,
        COALESCE(d1.steroid_day1, 0) AS itd_steroid_day1,
        COALESCE(d1.ventilation_day1, 0) AS itd_ventilation_day1,
        COALESCE(d1.crrt_day1, 0) AS itd_crrt_day1,
        COALESCE(d1.vasopressor_day1, 0) AS itd_vasopressor_day1,
        COUNT(fd.icu_day)::int AS itd_icu_day_rows,
        SUM(COALESCE(fd.iv_fluid_order_count::double precision, 0)) AS itd_intake_total_proxy_sum,
        SUM(fd.fluid_out_ml) AS itd_output_fluid_total_ml_proxy_sum,
        SUM(fd.net_fluid_balance_ml) AS itd_net_fluid_balance_proxy_sum,
        AVG(fd.iv_fluid_order_count::double precision) AS itd_intake_total_proxy_mean_daily,
        AVG(fd.fluid_out_ml) AS itd_output_fluid_total_ml_proxy_mean_daily,
        AVG(fd.net_fluid_balance_ml) AS itd_net_fluid_balance_proxy_mean_daily,
        MIN(fd.net_fluid_balance_ml) AS itd_net_fluid_balance_proxy_min_daily,
        MAX(fd.net_fluid_balance_ml) AS itd_net_fluid_balance_proxy_max_daily
    FROM derived.icu_treatment t
    LEFT JOIN hosp.admissions ad
      ON t.hadm_id = ad.hadm_id
    LEFT JOIN derived.icu_fluid_daily fd
      ON t.stay_id = fd.stay_id
    LEFT JOIN day1_flags d1
      ON t.stay_id = d1.stay_id
    GROUP BY
        t.subject_id, t.hadm_id, t.stay_id, t.icu_intime, t.icu_outtime, t.icu_los_days,
        t.steroid, t.mechanical_ventilation, t.non_invasive_ventilation, t.intubation, t.any_rrt, t.vasopressor,
        ad.admittime,
        d1.steroid_day1, d1.ventilation_day1, d1.crrt_day1, d1.vasopressor_day1
),
charlson_agg AS (
    SELECT
        ch.subject_id,
        ch.hadm_id,
        si.stay_id,
        MAX(ch.age) AS charlson_age,
        MAX(ch.age_score) AS charlson_age_score,
        MAX(ch.myocardial_infarct) AS charlson_myocardial_infarct,
        MAX(ch.congestive_heart_failure) AS charlson_congestive_heart_failure,
        MAX(ch.peripheral_vascular_disease) AS charlson_peripheral_vascular_disease,
        MAX(ch.cerebrovascular_disease) AS charlson_cerebrovascular_disease,
        MAX(ch.dementia) AS charlson_dementia,
        MAX(ch.chronic_pulmonary_disease) AS charlson_chronic_pulmonary_disease,
        MAX(ch.rheumatic_disease) AS charlson_rheumatic_disease,
        MAX(ch.peptic_ulcer_disease) AS charlson_peptic_ulcer_disease,
        MAX(ch.mild_liver_disease) AS charlson_mild_liver_disease,
        MAX(ch.diabetes_without_cc) AS charlson_diabetes_without_cc,
        MAX(ch.diabetes_with_cc) AS charlson_diabetes_with_cc,
        MAX(ch.paraplegia) AS charlson_paraplegia,
        MAX(ch.renal_disease) AS charlson_renal_disease,
        MAX(ch.malignant_cancer) AS charlson_malignant_cancer,
        MAX(ch.severe_liver_disease) AS charlson_severe_liver_disease,
        MAX(ch.metastatic_solid_tumor) AS charlson_metastatic_solid_tumor,
        MAX(ch.aids) AS charlson_aids,
        MAX(ch.charlson_comorbidity_index) AS charlson_comorbidity_index
    FROM derived.charlson ch
    INNER JOIN derived.cohort_all_icu si
      ON ch.hadm_id = si.hadm_id
    GROUP BY ch.subject_id, ch.hadm_id, si.stay_id
),
sirs_raw AS (
    SELECT
        c.subject_id,
        c.hadm_id,
        c.stay_id,
        MIN(CASE WHEN ce.itemid = 323761 THEN ce.valuenum END) AS temperature_min,
        MAX(CASE WHEN ce.itemid = 323761 THEN ce.valuenum END) AS temperature_max,
        MAX(CASE WHEN ce.itemid IN (320045, 320050, 320051) THEN ce.valuenum END) AS heart_rate_max,
        MAX(CASE WHEN ce.itemid = 320210 THEN ce.valuenum END) AS resp_rate_max,
        MIN(
            CASE
                WHEN lower(COALESCE(dl.label, '')) LIKE '%pco2%'
                  OR lower(COALESCE(dl.label, '')) LIKE '%partial pressure of carbon dioxide%'
                    THEN le.valuenum
                ELSE NULL
            END
        ) AS paco2_min,
        MIN(CASE WHEN le.itemid = 100016 THEN le.valuenum END) AS wbc_min,
        MAX(CASE WHEN le.itemid = 100016 THEN le.valuenum END) AS wbc_max,
        MAX(CASE WHEN le.itemid = 100028 THEN le.valuenum END) AS bands_max
    FROM derived.cohort_all_icu c
    LEFT JOIN icu.chartevents ce
      ON c.stay_id = ce.stay_id
     AND ce.charttime >= c.icu_intime
     AND ce.charttime <  c.icu_intime + INTERVAL '24 hours'
     AND ce.valuenum IS NOT NULL
     AND ce.itemid IN (323761, 320045, 320050, 320051, 320210)
    LEFT JOIN hosp.labevents le
      ON c.subject_id = le.subject_id
     AND c.hadm_id = le.hadm_id
     AND le.charttime >= c.icu_intime
     AND le.charttime <  c.icu_intime + INTERVAL '24 hours'
     AND le.valuenum IS NOT NULL
     AND (
          le.itemid IN (100016, 100028)
          OR le.itemid IN (
               SELECT itemid
               FROM hosp.d_labitems
               WHERE lower(label) LIKE '%pco2%'
                  OR lower(label) LIKE '%partial pressure of carbon dioxide%'
          )
     )
    LEFT JOIN hosp.d_labitems dl
      ON le.itemid = dl.itemid
    GROUP BY c.subject_id, c.hadm_id, c.stay_id
),
sirs_agg AS (
    SELECT
        s.subject_id,
        s.hadm_id,
        s.stay_id,
        MAX(s.sirs) AS sirs_score,
        MAX(s.temp_score) AS sirs_temp_score,
        MAX(s.heart_rate_score) AS sirs_heart_rate_score,
        MAX(s.resp_score) AS sirs_resp_score,
        MAX(s.wbc_score) AS sirs_wbc_score,
        MAX(sr.temperature_min) AS sirs_temperature_min,
        MAX(sr.temperature_max) AS sirs_temperature_max,
        MAX(sr.heart_rate_max) AS sirs_heart_rate_max,
        MAX(sr.resp_rate_max) AS sirs_resp_rate_max,
        MIN(sr.paco2_min) AS sirs_paco2_min,
        MIN(sr.wbc_min) AS sirs_wbc_min,
        MAX(sr.wbc_max) AS sirs_wbc_max,
        MAX(sr.bands_max) AS sirs_bands_max
    FROM derived.sirs s
    LEFT JOIN sirs_raw sr
      ON s.subject_id = sr.subject_id
     AND s.hadm_id = sr.hadm_id
     AND s.stay_id = sr.stay_id
    GROUP BY s.subject_id, s.hadm_id, s.stay_id
),
gcs_agg AS (
    SELECT
        c.subject_id,
        c.hadm_id,
        c.stay_id,
        MIN(g.gcs) AS gcs_min_icu,
        MAX(g.gcs) AS gcs_max_icu,
        AVG(g.gcs::double precision) AS gcs_mean_icu,
        MIN(g.gcs) FILTER (
            WHERE g.charttime >= c.icu_intime - INTERVAL '6 hours'
              AND g.charttime <= c.icu_intime + INTERVAL '24 hours'
        ) AS gcs_min_d1,
        MAX(g.gcs) FILTER (
            WHERE g.charttime >= c.icu_intime - INTERVAL '6 hours'
              AND g.charttime <= c.icu_intime + INTERVAL '24 hours'
        ) AS gcs_max_d1,
        AVG(g.gcs::double precision) FILTER (
            WHERE g.charttime >= c.icu_intime - INTERVAL '6 hours'
              AND g.charttime <= c.icu_intime + INTERVAL '24 hours'
        ) AS gcs_mean_d1
    FROM derived.cohort_all_icu c
    LEFT JOIN derived.gcs g
      ON c.stay_id = g.stay_id
    GROUP BY c.subject_id, c.hadm_id, c.stay_id, c.icu_intime
),
kdigo_agg AS (
    SELECT
        c.subject_id,
        c.hadm_id,
        c.stay_id,
        MAX(cb.scr_baseline) AS kdigo_scr_baseline,
        GREATEST(
            COALESCE(MAX(k.aki_stage_creat), 0),
            COALESCE(
                MAX(
                    CASE
                        WHEN t.stay_id IS NULL THEN NULL
                        WHEN COALESCE(t.crrt, 0) = 1 OR COALESCE(t.any_rrt, 0) = 1 THEN 3
                        ELSE 0
                    END
                ),
                0
            )
        ) AS kdigo_aki_stage_max,
        MAX(k.aki_stage_smoothed) AS kdigo_aki_stage_smoothed_max,
        MAX(k.aki_stage_creat) AS kdigo_aki_stage_creat_max,
        NULL::int AS kdigo_aki_stage_uo_max,
        MAX(
            CASE
                WHEN t.stay_id IS NULL THEN NULL
                WHEN COALESCE(t.crrt, 0) = 1 OR COALESCE(t.any_rrt, 0) = 1 THEN 3
                ELSE 0
            END
        ) AS kdigo_aki_stage_crrt_max,
        MAX(k.creat) AS kdigo_creat_max,
        MIN(k.creat) AS kdigo_creat_min,
        AVG(k.creat) AS kdigo_creat_mean,
        MAX(k.aki_stage_creat) FILTER (
            WHERE k.charttime >= c.icu_intime - INTERVAL '6 hours'
              AND k.charttime <= c.icu_intime + INTERVAL '24 hours'
        ) AS kdigo_aki_stage_d1_max
    FROM derived.cohort_all_icu c
    LEFT JOIN derived.kdigo_stages k
      ON c.stay_id = k.stay_id
    LEFT JOIN derived.creatinine_baseline cb
      ON c.hadm_id = cb.hadm_id
    LEFT JOIN derived.icu_treatment t
      ON c.stay_id = t.stay_id
    GROUP BY c.subject_id, c.hadm_id, c.stay_id, c.icu_intime
)
SELECT
    c.*,
    fa.*,
    ita.itd_icu_episode_count,
    ita.itd_icu_days_total,
    ita.itd_first_icu_intime,
    ita.itd_last_icu_outtime,
    ita.itd_used_steroid_icu,
    ita.itd_used_ventilation_icu,
    ita.itd_used_crrt_icu,
    ita.itd_used_vasopressor_icu,
    ita.itd_steroid_day1,
    ita.itd_ventilation_day1,
    ita.itd_crrt_day1,
    ita.itd_vasopressor_day1,
    ita.itd_icu_day_rows,
    ita.itd_intake_total_proxy_sum,
    ita.itd_output_fluid_total_ml_proxy_sum,
    ita.itd_net_fluid_balance_proxy_sum,
    ita.itd_intake_total_proxy_mean_daily,
    ita.itd_output_fluid_total_ml_proxy_mean_daily,
    ita.itd_net_fluid_balance_proxy_mean_daily,
    ita.itd_net_fluid_balance_proxy_min_daily,
    ita.itd_net_fluid_balance_proxy_max_daily,
    cha.charlson_age,
    cha.charlson_age_score,
    cha.charlson_myocardial_infarct,
    cha.charlson_congestive_heart_failure,
    cha.charlson_peripheral_vascular_disease,
    cha.charlson_cerebrovascular_disease,
    cha.charlson_dementia,
    cha.charlson_chronic_pulmonary_disease,
    cha.charlson_rheumatic_disease,
    cha.charlson_peptic_ulcer_disease,
    cha.charlson_mild_liver_disease,
    cha.charlson_diabetes_without_cc,
    cha.charlson_diabetes_with_cc,
    cha.charlson_paraplegia,
    cha.charlson_renal_disease,
    cha.charlson_malignant_cancer,
    cha.charlson_severe_liver_disease,
    cha.charlson_metastatic_solid_tumor,
    cha.charlson_aids,
    cha.charlson_comorbidity_index,
    sia.sirs_score,
    sia.sirs_temp_score,
    sia.sirs_heart_rate_score,
    sia.sirs_resp_score,
    sia.sirs_wbc_score,
    sia.sirs_temperature_min,
    sia.sirs_temperature_max,
    sia.sirs_heart_rate_max,
    sia.sirs_resp_rate_max,
    sia.sirs_paco2_min,
    sia.sirs_wbc_min,
    sia.sirs_wbc_max,
    sia.sirs_bands_max,
    gca.gcs_min_icu,
    gca.gcs_max_icu,
    gca.gcs_mean_icu,
    gca.gcs_min_d1,
    gca.gcs_max_d1,
    gca.gcs_mean_d1,
    kda.kdigo_scr_baseline,
    kda.kdigo_aki_stage_max,
    kda.kdigo_aki_stage_smoothed_max,
    kda.kdigo_aki_stage_creat_max,
    kda.kdigo_aki_stage_uo_max,
    kda.kdigo_aki_stage_crrt_max,
    kda.kdigo_creat_max,
    kda.kdigo_creat_min,
    kda.kdigo_creat_mean,
    kda.kdigo_aki_stage_d1_max
FROM derived.cohort_all_icu c
LEFT JOIN final_agg fa
  ON c.subject_id = fa.agg_subject_id
 AND c.hadm_id = fa.agg_hadm_id
 AND c.stay_id = fa.agg_stay_id
LEFT JOIN icu_treat_agg ita
  ON c.subject_id = ita.subject_id
 AND c.hadm_id = ita.hadm_id
 AND c.stay_id = ita.stay_id
LEFT JOIN charlson_agg cha
  ON c.subject_id = cha.subject_id
 AND c.hadm_id = cha.hadm_id
 AND c.stay_id = cha.stay_id
LEFT JOIN sirs_agg sia
  ON c.subject_id = sia.subject_id
 AND c.hadm_id = sia.hadm_id
 AND c.stay_id = sia.stay_id
LEFT JOIN gcs_agg gca
  ON c.subject_id = gca.subject_id
 AND c.hadm_id = gca.hadm_id
 AND c.stay_id = gca.stay_id
LEFT JOIN kdigo_agg kda
  ON c.subject_id = kda.subject_id
 AND c.hadm_id = kda.hadm_id
 AND c.stay_id = kda.stay_id
ORDER BY c.subject_id, c.hadm_id, c.stay_id;

CREATE INDEX IF NOT EXISTS idx_nwicu_final_select_stay
  ON derived.cohort_all_icu_lab_final_select (stay_id);
CREATE INDEX IF NOT EXISTS idx_nwicu_final_select_hadm
  ON derived.cohort_all_icu_lab_final_select (hadm_id);
CREATE INDEX IF NOT EXISTS idx_nwicu_final_select_subject
  ON derived.cohort_all_icu_lab_final_select (subject_id);
