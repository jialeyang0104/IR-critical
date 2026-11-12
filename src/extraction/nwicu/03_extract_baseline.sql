-- Read-only supplement for baseline tables; existing formal SOFA and cohort remain unchanged.
WITH ranked AS (
    SELECT i.*,row_number() OVER (PARTITION BY subject_id ORDER BY intime,stay_id) AS rn
    FROM icu.icustays i
), eligible AS MATERIALIZED (
    SELECT r.* FROM ranked r
    JOIN derived.cohort_all_icu_lab_final_select c USING (stay_id)
    JOIN hosp.admissions a ON a.hadm_id=r.hadm_id
    WHERE r.rn=1 AND c.age>=18
      AND r.outtime>=r.intime+interval '24 hours'
      AND a.dischtime>=r.intime+interval '24 hours'
), vital_times AS (
    SELECT e.stay_id,ce.charttime,
        avg(ce.valuenum) FILTER (WHERE ce.itemid=320045 AND ce.valuenum>0 AND ce.valuenum<300) AS hr,
        avg(ce.valuenum) FILTER (WHERE ce.itemid=320210 AND ce.valuenum>0 AND ce.valuenum<100) AS rr,
        avg((ce.valuenum-32)*5.0/9.0) FILTER (WHERE ce.itemid=323761
            AND (ce.valuenum-32)*5.0/9.0 BETWEEN 25 AND 45) AS temp_c,
        avg(ce.valuenum) FILTER (WHERE ce.itemid=320179 AND ce.valuenum>0 AND ce.valuenum<400) AS sbp,
        avg(ce.valuenum) FILTER (WHERE ce.itemid=320180 AND ce.valuenum>0 AND ce.valuenum<300) AS dbp
    FROM eligible e JOIN icu.chartevents ce ON ce.stay_id=e.stay_id
      AND ce.charttime BETWEEN e.intime AND e.intime+interval '24 hours'
      AND ce.itemid IN (320045,320210,323761,320179,320180) AND ce.valuenum IS NOT NULL
    GROUP BY e.stay_id,ce.charttime
), vitals AS (
    SELECT stay_id,avg(hr) AS heart_rate_mean,avg(rr) AS resp_rate_mean,
        avg(temp_c) AS temperature_mean,
        avg((sbp+2*dbp)/3) FILTER (WHERE sbp>=dbp) AS mbp_mean,
        min(temp_c) AS temp_min,max(temp_c) AS temp_max,max(hr) AS hr_max,max(rr) AS rr_max
    FROM vital_times GROUP BY stay_id
), diagnoses AS (
    SELECT e.stay_id,
        max(((d.icd_version=9 AND d.icd_code ~ '^40[1-5]')
            OR (d.icd_version=10 AND substr(d.icd_code,1,3) IN ('I10','I11','I15','I16')))::int) AS hypertension,
        max(((d.icd_version=9 AND d.icd_code='V5867')
            OR (d.icd_version=10 AND d.icd_code='Z794'))::int) AS longterm_insulin,
        max((d.icd_version=10 AND d.icd_code='Z7984')::int) AS longterm_antidiabetics
    FROM eligible e LEFT JOIN hosp.diagnoses_icd d ON d.hadm_id=e.hadm_id GROUP BY e.stay_id
), steroids AS (
    SELECT e.stay_id,max((lower(p.drug) ~
        '(dexamethasone|hydrocortisone|methylprednisolone|prednisone|prednisolone|cortisone|betamethasone|triamcinolone|fludrocortisone)')::int)
        AS glucocorticoids_1d
    FROM eligible e LEFT JOIN hosp.prescriptions p ON p.hadm_id=e.hadm_id
      AND p.starttime<=e.intime+interval '24 hours'
      AND coalesce(p.stoptime,p.starttime)>=e.intime-interval '6 hours'
    GROUP BY e.stay_id
)
SELECT e.subject_id,e.hadm_id,e.stay_id,
    v.heart_rate_mean,v.resp_rate_mean,v.temperature_mean,v.mbp_mean,
    c."fd_WBC_mean" AS wbc_mean,c."fd_PLT_mean" AS platelets_mean,
    c."fd_HGB_mean" AS hemoglobin_mean,c."fd_UREA_mean" AS bun_mean,
    c."fd_CRE_mean" AS creatinine_mean,c."fd_NA_mean" AS sodium_mean,
    c."fd_K_mean" AS potassium_mean,c."fd_INR_mean" AS inr_mean,
    c."fd_PT_mean" AS pt_mean,c."fd_ALT_mean" AS alt_mean,c."fd_AST_mean" AS ast_mean,
    c.gcs_min_d1 AS gcs_min,
    -- Temperature is Fahrenheit in NWICU. Recalculate the temperature and pulse
    -- components from validated day-1 vital signs; arterial BP is not heart rate.
    coalesce(((v.temp_min<36 OR v.temp_max>38)::int),0)
      +coalesce((v.hr_max>90)::int,0)+coalesce((v.rr_max>20)::int,0)
      +coalesce(c.sirs_wbc_score,0) AS sirs,
    c.charlson_myocardial_infarct AS myocardial_infarct,
    c.charlson_congestive_heart_failure AS congestive_heart_failure,
    c.charlson_chronic_pulmonary_disease AS chronic_pulmonary_disease,
    dx.hypertension,dx.longterm_insulin,dx.longterm_antidiabetics,
    coalesce(st.glucocorticoids_1d,0) AS glucocorticoids_1d,
    c.itd_ventilation_day1 AS ventilation_day1,c.itd_vasopressor_day1 AS vasopressor_day1
FROM eligible e
JOIN derived.cohort_all_icu_lab_final_select c USING (stay_id)
LEFT JOIN vitals v USING (stay_id)
LEFT JOIN diagnoses dx USING (stay_id)
LEFT JOIN steroids st USING (stay_id)
ORDER BY e.stay_id;
