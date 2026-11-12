# Restricted inputs and source definitions

Place authorized CSVs under `$IR_REVISION_ROOT/data/restricted/`. The directory is not part of this distribution. Raw extracts are `mimiciv_analysis.csv`, `mimiciii_external_validation.csv`, and `nwicu_external_validation.csv`. Script 00a preserves these files and writes `mimiciv_patient_level.csv`, `mimiciii_patient_level.csv`, and `nwicu_patient_level.csv`. It joins each raw stay to its rank in the complete ICU table and retains rank=1 before applying downstream laboratory selection. Patient uniqueness is asserted. Original source row and original age are retained only in the restricted derivatives.

|Fields (preferred name; aliases in core code)|Format / meaning|
|---|---|
|subject_id, hadm_id, stay_id or icustay_id|Restricted linkage identifiers; never distributed|
|age or admission_age; age_original; age_deidentified|Years; original value and explicit flag; MIMIC-III codes ≥300 become 90; 91.4 is retained only in historical analysis records|
|gender or sex|Sex category|
|sofa; charlson_comorbidity_index; diabetes|Formal adjustment; NWICU uses the adapted SOFA described below|
|glucose; hba1c; tg/triglycerides; hdl|Window means: mg/dL; %; mg/dL; mg/dL|
|icu_intime; icu_outtime; dischtime; deathtime; dod|Source timestamps or day-level death dates; retained under source access restrictions|
|hospital_outcome; icu_outcome|Death in the index hospitalization or first ICU; 0/1|
|outcome_28d; time_28d_days; mortality_28d_time_source|Verified day-28 event, follow-up duration and date/timestamp provenance|
|admission_type; first_careunit; diagnosis_category; dbsource|Case-mix/selection and MIMIC-III CareVue fields, queried by 00a|
|liver_failure; renal_failure; longterm_insulin; longterm_antidiabetics|Optional exclusion fields; missing availability is reported, not fabricated|
|coded_sepsis; coded_head_trauma|Discharge-code descriptions from 00a; no onset-time adjudication|
|covid_u071; coronavirus_b9729; diagnosis_available; steroid_24h|NWICU diagnosis availability and documented systemic EMAR administration in [ICU admission,+24h)|

Raw eligibility uses age ≥18, ICU duration ≥24h and discharge time ≥ICU admission+24h. Complete components are required after patient selection. Cleaning ranges: glucose 10–2,000, HbA1c 2–25, TG 1–10,000 and HDL 1–500. SHR = glucose/(28.7×HbA1c−46.7); TyG = ln(TG×glucose/2); AIP = log10[(TG/88.57)/(HDL/38.67)]. Formal covariates are complete in the current final cohorts. Auxiliary prediction filling is learned in training folds.

The MIMIC-IV baseline/prediction extract also supplies vital signs, blood count, chemistry, SAPS II, APS III, OASIS, GCS, SIRS, comorbidities and treatment indicators. Exact fields and fixed Basic Model variable sets are declared in shared script 21 and its aggregate variable-set CSV. All analysis-specific candidate aliases are explicit in shared core 18.

## Observed baseline and selection tables

The baseline-only read-only export is `src/extraction/01_extract_baseline.py`, using `mimiciii/03_extract_baseline.sql` and `nwicu/03_extract_baseline.sql`. It saves `mimiciii_baseline_supplement.csv` and `nwicu_baseline_supplement.csv` under the ignored restricted directory. Their identifiers must match the current eligible first-stay inputs exactly. The main regression input files are not modified. `src/analysis/03_descriptives.R` combines these supplements with the existing three cohort inputs.

Table 1 contains the same 42 characteristics for MIMIC-IV, MIMIC-III and NWICU; Table S2 reports their observed missing counts in each final common-index cohort. Table S3 uses 48 characteristics (the 42 baseline rows plus SAPS II, APS III, OASIS and the three mortality outcomes) in three cohort-specific included/excluded panels. Percentages use the nonmissing denominator for each variable. Laboratory and index summaries in excluded patients use the available measurements, even if another required component is missing.

Table S4 stratifies the MIMIC-IV final cohort by index-hospitalization vital status: 2,735 survivors and 321 non-survivors. Its 42 variables and overall summaries match Table 1. Continuous variables are observed medians (interquartile ranges), and categorical values are n (%) of nonmissing observations. Two-sided P values use the Mann–Whitney U test for continuous variables and Pearson's chi-square or Fisher's exact test for categorical variables (Fisher's test when any expected cell count is below 5). The source `00_CR_M4_survivor_nonsurvivor_baseline.csv` retains total and variable-specific nonmissing denominators for both groups.

Hospital length of stay is `(dischtime - admittime)` in exact days; ICU length of stay is `(icu_outtime - icu_intime)` in exact days. These descriptive values do not replace the post-ICU-admission durations used for eligibility or regression. Vitals use 0–24 h; routine laboratory means use −6 to +24 h. NWICU temperature item 323761 is Fahrenheit and is converted to Celsius. NWICU mean blood pressure is estimated from paired cuff systolic/diastolic readings as `(SBP + 2*DBP)/3`; SIRS temperature and pulse criteria use the first-day temperature and pulse items, not arterial blood pressure values. Glucocorticoid flags use overlapping early-window prescriptions; MIMIC-III prescription dates have day-level precision. Ventilation and vasopressor flags use early treatment records from the available cohort-specific source tables.

MIMIC-III lacks a comparable dedicated ICD-9 code for long-term oral antidiabetic use, so that definition is marked unavailable rather than assigned zero or 100% patient-level missingness. GCS and INR have no observed values in the current final NWICU sample; their observed missingness is 100%. SAPS II, APS III and OASIS are not available in the NWICU extraction and are marked NA in Table S3. No baseline values are imputed. Full field mappings and availability are in `00_CR_ALL_baseline_source_fields.csv` and `00_CR_ALL_baseline_field_availability.csv`.

## Alternative laboratory definitions

Provide `mimiciv_lab_timing_definitions.csv`, `mimiciii_lab_timing_definitions.csv` and `nwicu_lab_timing_definitions.csv`. Each stay–analyte pair is unique, with `stay_id, analyte, first_available_value, closest_to_icu_value, early_window_mean`. Analyte is glucose/hba1c/tg/hdl. First means earliest valid value in −6 to +24h; closest minimizes absolute time from ICU admission; mean is the arithmetic window mean. The MIMIC-IV extract may alternatively contain embedded `*_first_value` and `*_closest_value` fields. Cohort extraction SQL includes event/definition queries; retain outputs in the restricted directory. Timing comparisons use patients with all indices under all definitions.

## Database prerequisites and derived SOFA

The explicit local specifications are `mimiciv` on localhost:5432 (mimiciv_icu/hosp/derived), `mimiciii` on localhost:5433 (mimiciii schema) and `nwicu` on localhost:5432 (icu/hosp/derived). Adapt those specifications to a licensed local installation. Prepare MIMIC upstream derived tables referenced in each SQL (ICU-stay detail, first-day severity and Charlson) using the official MIT-LCP mimic-code concepts. The included NWICU upstream scripts build the local cohort/covariate table in the order stated in docs/RUN_ORDER.md. The package does not include database dumps.

NWICU uses `derived.nwicu_sofa_day1_adapted`. The included construction SQL documents respiratory SpO2/FiO2 substitution, creatinine-only renal scoring, MAP and EMAR-confirmed prescription rates, use of derived GCS when available, and normal scoring for missing components in the total score. Component availability/provenance are retained. This adapted score is explicitly identified in the manuscript.

State 0 is right censoring, 1 target death, 2 alive departure for stay-specific competing-risk outcomes. MIMIC 28-day follow-up extends past discharge to recorded death or day 28. NWICU 28-day deaths have hospital-system ascertainment. Events with absent exact death times use departure-time proxies for stay-specific endpoints; date-only 28-day precision and origin-specific exclusions are reported in aggregate audit CSVs.
