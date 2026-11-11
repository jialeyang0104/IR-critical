-- ====================================================================================================
-- charlson.sql (nwicu adaptation)
-- Charlson Comorbidity Index based on ICD-10 codes (nwicu has ICD-10 only)
-- Reference: Quan H, et al. Med Care. 2005; 43(11):1130-9.
-- ====================================================================================================
DROP TABLE IF EXISTS derived.charlson;
CREATE TABLE derived.charlson AS
WITH diag AS (
  SELECT
    hadm_id,
    icd_code AS icd10_code
  FROM hosp.diagnoses_icd
  WHERE icd_version = 10
),
com AS (
  SELECT
    ad.hadm_id,
    -- Myocardial infarction
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('I21','I22') OR SUBSTR(icd10_code,1,4)='I252' THEN 1 ELSE 0 END) AS myocardial_infarct,
    -- Congestive heart failure
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('I43','I50')
              OR SUBSTR(icd10_code,1,4) IN ('I099','I110','I130','I132','I255','I420','I425','I426','I427','I428','I429','P290')
             THEN 1 ELSE 0 END) AS congestive_heart_failure,
    -- Peripheral vascular disease
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('I70','I71')
              OR SUBSTR(icd10_code,1,4) IN ('I731','I738','I739','I771','I790','I792','K551','K558','K559','Z958','Z959')
             THEN 1 ELSE 0 END) AS peripheral_vascular_disease,
    -- Cerebrovascular disease
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('G45','G46')
              OR SUBSTR(icd10_code,1,3) BETWEEN 'I60' AND 'I69'
              OR SUBSTR(icd10_code,1,4)='H340'
             THEN 1 ELSE 0 END) AS cerebrovascular_disease,
    -- Dementia
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('F00','F01','F02','F03','G30')
              OR SUBSTR(icd10_code,1,4) IN ('F051','G311')
             THEN 1 ELSE 0 END) AS dementia,
    -- Chronic pulmonary disease
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) BETWEEN 'J40' AND 'J47'
              OR SUBSTR(icd10_code,1,3) BETWEEN 'J60' AND 'J67'
              OR SUBSTR(icd10_code,1,4) IN ('I278','I279','J684','J701','J703')
             THEN 1 ELSE 0 END) AS chronic_pulmonary_disease,
    -- Rheumatic disease
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('M05','M06','M32','M33','M34')
              OR SUBSTR(icd10_code,1,4) IN ('M315','M351','M353','M360')
             THEN 1 ELSE 0 END) AS rheumatic_disease,
    -- Peptic ulcer disease
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('K25','K26','K27','K28') THEN 1 ELSE 0 END) AS peptic_ulcer_disease,
    -- Mild liver disease
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('B18','K73','K74')
              OR SUBSTR(icd10_code,1,4) IN ('K700','K701','K702','K703','K709','K713','K714','K715','K717','K760','K762','K763','K764','K768','K769','Z944')
             THEN 1 ELSE 0 END) AS mild_liver_disease,
    -- Diabetes without chronic complication
    MAX(CASE WHEN SUBSTR(icd10_code,1,4) IN ('E100','E101','E106','E108','E109','E110','E111','E116','E118','E119',
                                              'E120','E121','E126','E128','E129','E130','E131','E136','E138','E139',
                                              'E140','E141','E146','E148','E149')
             THEN 1 ELSE 0 END) AS diabetes_without_cc,
    -- Diabetes with chronic complication
    MAX(CASE WHEN SUBSTR(icd10_code,1,4) IN ('E102','E103','E104','E105','E107','E112','E113','E114','E115','E117',
                                              'E122','E123','E124','E125','E127','E132','E133','E134','E135','E137',
                                              'E142','E143','E144','E145','E147')
             THEN 1 ELSE 0 END) AS diabetes_with_cc,
    -- Hemiplegia / paraplegia
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('G81','G82')
              OR SUBSTR(icd10_code,1,4) IN ('G041','G114','G801','G802','G830','G831','G832','G833','G834','G839')
             THEN 1 ELSE 0 END) AS paraplegia,
    -- Renal disease
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('N18','N19')
              OR SUBSTR(icd10_code,1,4) IN ('I120','I131','N032','N033','N034','N035','N036','N037',
                                             'N052','N053','N054','N055','N056','N057','N250','Z490','Z491','Z492','Z940','Z992')
             THEN 1 ELSE 0 END) AS renal_disease,
    -- Malignant cancer (excluding skin)
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('C43','C88')
              OR SUBSTR(icd10_code,1,3) BETWEEN 'C00' AND 'C26'
              OR SUBSTR(icd10_code,1,3) BETWEEN 'C30' AND 'C34'
              OR SUBSTR(icd10_code,1,3) BETWEEN 'C37' AND 'C41'
              OR SUBSTR(icd10_code,1,3) BETWEEN 'C45' AND 'C58'
              OR SUBSTR(icd10_code,1,3) BETWEEN 'C60' AND 'C76'
              OR SUBSTR(icd10_code,1,3) BETWEEN 'C81' AND 'C85'
              OR SUBSTR(icd10_code,1,3) BETWEEN 'C90' AND 'C97'
             THEN 1 ELSE 0 END) AS malignant_cancer,
    -- Moderate / severe liver disease
    MAX(CASE WHEN SUBSTR(icd10_code,1,4) IN ('I850','I859','I864','I982','K704','K711','K721','K729','K765','K766','K767')
             THEN 1 ELSE 0 END) AS severe_liver_disease,
    -- Metastatic solid tumor
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('C77','C78','C79','C80') THEN 1 ELSE 0 END) AS metastatic_solid_tumor,
    -- AIDS/HIV
    MAX(CASE WHEN SUBSTR(icd10_code,1,3) IN ('B20','B21','B22','B24') THEN 1 ELSE 0 END) AS aids
  FROM hosp.admissions AS ad
  LEFT JOIN diag ON ad.hadm_id = diag.hadm_id
  GROUP BY ad.hadm_id
),
ag AS (
  SELECT
    hadm_id,
    age,
    CASE
      WHEN age <= 50 THEN 0
      WHEN age <= 60 THEN 1
      WHEN age <= 70 THEN 2
      WHEN age <= 80 THEN 3
      ELSE 4
    END AS age_score
  FROM derived.age
)
SELECT
  ad.subject_id,
  ad.hadm_id,
  ag.age,
  ag.age_score,
  myocardial_infarct,
  congestive_heart_failure,
  peripheral_vascular_disease,
  cerebrovascular_disease,
  dementia,
  chronic_pulmonary_disease,
  rheumatic_disease,
  peptic_ulcer_disease,
  mild_liver_disease,
  diabetes_without_cc,
  diabetes_with_cc,
  paraplegia,
  renal_disease,
  malignant_cancer,
  severe_liver_disease,
  metastatic_solid_tumor,
  aids,
  -- Charlson score (original Charlson 1987 weights)
  COALESCE(ag.age_score, 0)
  + COALESCE(myocardial_infarct, 0)
  + COALESCE(congestive_heart_failure, 0)
  + COALESCE(peripheral_vascular_disease, 0)
  + COALESCE(cerebrovascular_disease, 0)
  + COALESCE(dementia, 0)
  + COALESCE(chronic_pulmonary_disease, 0)
  + COALESCE(rheumatic_disease, 0)
  + COALESCE(peptic_ulcer_disease, 0)
  + GREATEST(COALESCE(mild_liver_disease,0), 3 * COALESCE(severe_liver_disease,0))
  + GREATEST(2 * COALESCE(diabetes_with_cc,0), COALESCE(diabetes_without_cc,0))
  + GREATEST(2 * COALESCE(malignant_cancer,0), 6 * COALESCE(metastatic_solid_tumor,0))
  + 2 * COALESCE(paraplegia, 0)
  + 2 * COALESCE(renal_disease, 0)
  + 6 * COALESCE(aids, 0)
  AS charlson_comorbidity_index
FROM hosp.admissions AS ad
LEFT JOIN com   ON ad.hadm_id = com.hadm_id
LEFT JOIN ag    ON ad.hadm_id = ag.hadm_id
ORDER BY ad.subject_id, ad.hadm_id;
