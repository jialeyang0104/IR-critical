-- ====================================================================================================
-- 07_sepsis_labs.sql
-- 目的：为 derived.cohort_all_icu 患者提取实验室指标
-- 
-- Part A：ICU 入住 -6h ~ +24h（Day 1）内的全套指标宽表（一患者一行）
-- Part B：ICU 入住 -6h ~ +168h 按天（d1~d7）提取
--           动态指标：LDH、Creatinine、Platelet、Abs Mono/Neu/Lym（平均值）
-- ====================================================================================================

-- ============================================================
-- PART A: ICU Day 1 全套实验室指标宽表 (-6h ~ +24h)
-- ============================================================
DROP TABLE IF EXISTS derived.all_icu_labs_d1;
CREATE TABLE derived.all_icu_labs_d1 AS
WITH lab_window AS (
  SELECT
    c.subject_id,
    c.hadm_id,
    c.stay_id,
    c.icu_intime,
    l.itemid,
    l.valuenum,
    l.charttime,
    -- 距 ICU 入时的小时数
    EXTRACT(EPOCH FROM (l.charttime - c.icu_intime)) / 3600.0 AS hours_from_icu
  FROM derived.cohort_all_icu AS c
  JOIN hosp.labevents AS l ON l.subject_id = c.subject_id
  WHERE
    -- 时间窗口：ICU 入住前 6h 至入住后 24h
    l.charttime BETWEEN (c.icu_intime - INTERVAL '6 hours')
                    AND (c.icu_intime + INTERVAL '24 hours')
    AND l.valuenum IS NOT NULL
    AND l.valuenum <> 9999999
    AND l.itemid IN (
      -- === Blood Gas ===
      100339,  -- pH
      100064,  -- pH
      100262,  -- pH
      100029,  -- pO2 (Oxygen)
      100109,  -- Oxygen Saturation
      100031,  -- Lactate (Blood Gas)
      100044,  -- Free Calcium (iCa)
      100150,  -- Carboxyhemoglobin
      100145,  -- Methemoglobin
      -- === Chemistry: 肝功能 ===
      100042,  -- ALT
      100032,  -- AST
      100033,  -- Alkaline Phosphatase
      100020,  -- Bilirubin, Total
      100049,  -- Bilirubin, Direct
      100021,  -- Albumin
      100022,  -- Total Protein (via label, see note)
      100161,  -- GGT (Gamma Glutamyltransferase)
      -- === Chemistry: 肾功能 ===
      100002,  -- Creatinine
      100004,  -- Urea Nitrogen (BUN)
      -- === Chemistry: 电解质 ===
      100010,  -- Sodium
      100011,  -- Potassium
      100012,  -- Chloride
      100003,  -- Calcium, Total
      100009,  -- Magnesium
      100015,  -- Phosphate
      100013,  -- Total CO2 (Bicarbonate)
      100040,  -- Anion Gap
      100041,  -- Anion Gap (alt)
      -- === Chemistry: 炎症 & 生物标志物 ===
      100053,  -- CRP
      100251,  -- High-Sensitivity CRP
      100052,  -- Ferritin
      100055,  -- LDH (Lactate Dehydrogenase)
      100075,  -- D-Dimer
      100057,  -- Troponin T
      100059,  -- Troponin I
      100071,  -- NTproBNP
      -- === Chemistry: 血糖与血脂 ===
      100001,  -- Glucose (Chemistry)
      100045,  -- Glucose (Blood Gas)
      100088,  -- Cholesterol, Total
      100086,  -- Cholesterol, HDL
      100090,  -- Cholesterol, LDL, Calculated
      100248,  -- Cholesterol, LDL, Measured
      100056,  -- Triglycerides
      100079,  -- Hemoglobin A1c
      -- === Hematology: CBC ===
      100016,  -- WBC
      100019,  -- RBC
      100007,  -- Hemoglobin
      100006,  -- Hematocrit
      100014,  -- Platelet Count
      100017,  -- MCV
      100005,  -- MCH
      100018,  -- MCHC
      100008,  -- RDW
      100028,  -- Immature Granulocytes
      -- === Hematology: 白细胞分类 (%) ===
      100022,  -- Neutrophils %
      100023,  -- Lymphocytes %
      100026,  -- Monocytes %
      100024,  -- Eosinophils %
      100025,  -- Basophils %
      -- === Hematology: 白细胞绝对计数 ===
      100038,  -- Absolute Neutrophil Count
      100037,  -- Absolute Lymphocyte Count
      100039,  -- Absolute Monocyte Count
      100036,  -- Absolute Eosinophil Count
      100035,  -- Absolute Basophil Count
      -- === Hematology: 凝血 ===
      100034,  -- PT
      100030,  -- INR(PT)
      100046,  -- PTT
      100048   -- Fibrinogen
    )
)
SELECT
  subject_id,
  hadm_id,
  stay_id,
  icu_intime,
  -- === 血气 ===
  AVG(CASE WHEN itemid IN (100339, 100064, 100262) THEN valuenum END) AS bg_ph,
  AVG(CASE WHEN itemid = 100029 THEN valuenum END) AS bg_po2,
  AVG(CASE WHEN itemid = 100109 THEN valuenum END) AS bg_o2sat,
  AVG(CASE WHEN itemid = 100031 THEN valuenum END) AS bg_lactate,
  AVG(CASE WHEN itemid = 100044 THEN valuenum END) AS bg_free_calcium,
  -- === 肝功能 ===
  AVG(CASE WHEN itemid = 100042 THEN valuenum END) AS alt,
  AVG(CASE WHEN itemid = 100032 THEN valuenum END) AS ast,
  AVG(CASE WHEN itemid = 100033 THEN valuenum END) AS alp,
  AVG(CASE WHEN itemid = 100020 THEN valuenum END) AS bilirubin_total,
  AVG(CASE WHEN itemid = 100049 THEN valuenum END) AS bilirubin_direct,
  AVG(CASE WHEN itemid = 100021 THEN valuenum END) AS albumin,
  AVG(CASE WHEN itemid = 100161 THEN valuenum END) AS ggt,
  -- === 肾功能 ===
  AVG(CASE WHEN itemid = 100002 THEN valuenum END) AS creatinine,
  AVG(CASE WHEN itemid = 100004 THEN valuenum END) AS bun,
  -- === 电解质 ===
  AVG(CASE WHEN itemid = 100010 THEN valuenum END) AS sodium,
  AVG(CASE WHEN itemid = 100011 THEN valuenum END) AS potassium,
  AVG(CASE WHEN itemid = 100012 THEN valuenum END) AS chloride,
  AVG(CASE WHEN itemid = 100003 THEN valuenum END) AS calcium_total,
  AVG(CASE WHEN itemid = 100009 THEN valuenum END) AS magnesium,
  AVG(CASE WHEN itemid = 100015 THEN valuenum END) AS phosphate,
  AVG(CASE WHEN itemid = 100013 THEN valuenum END) AS total_co2,
  AVG(CASE WHEN itemid IN (100040, 100041) THEN valuenum END) AS anion_gap,
  -- === 炎症与生物标志物 ===
  AVG(CASE WHEN itemid IN (100053, 100251) THEN valuenum END) AS crp,
  AVG(CASE WHEN itemid = 100052 THEN valuenum END) AS ferritin,
  AVG(CASE WHEN itemid = 100055 THEN valuenum END) AS ldh,
  AVG(CASE WHEN itemid = 100075 THEN valuenum END) AS d_dimer,
  AVG(CASE WHEN itemid = 100057 THEN valuenum END) AS troponin_t,
  AVG(CASE WHEN itemid = 100059 THEN valuenum END) AS troponin_i,
  AVG(CASE WHEN itemid = 100071 THEN valuenum END) AS ntprobnp,
  -- === 血糖与血脂 ===
  AVG(CASE WHEN itemid IN (100001, 100045) THEN valuenum END) AS glucose,
  AVG(CASE WHEN itemid = 100088 THEN valuenum END) AS cholesterol_total,
  AVG(CASE WHEN itemid = 100086 THEN valuenum END) AS hdl,
  AVG(CASE WHEN itemid IN (100090, 100248) THEN valuenum END) AS ldl,
  AVG(CASE WHEN itemid = 100056 THEN valuenum END) AS triglycerides,
  AVG(CASE WHEN itemid = 100079 THEN valuenum END) AS hba1c,
  -- === CBC 常规 ===
  AVG(CASE WHEN itemid = 100016 THEN valuenum END) AS wbc,
  AVG(CASE WHEN itemid = 100019 THEN valuenum END) AS rbc,
  AVG(CASE WHEN itemid = 100007 THEN valuenum END) AS hemoglobin,
  AVG(CASE WHEN itemid = 100006 THEN valuenum END) AS hematocrit,
  AVG(CASE WHEN itemid = 100014 THEN valuenum END) AS platelet,
  AVG(CASE WHEN itemid = 100017 THEN valuenum END) AS mcv,
  AVG(CASE WHEN itemid = 100005 THEN valuenum END) AS mch,
  AVG(CASE WHEN itemid = 100018 THEN valuenum END) AS mchc,
  AVG(CASE WHEN itemid = 100008 THEN valuenum END) AS rdw,
  AVG(CASE WHEN itemid = 100028 THEN valuenum END) AS immature_granulocytes,
  -- === 白细胞分类 (%) ===
  AVG(CASE WHEN itemid = 100022 THEN valuenum END) AS neutrophils_pct,
  AVG(CASE WHEN itemid = 100023 THEN valuenum END) AS lymphocytes_pct,
  AVG(CASE WHEN itemid = 100026 THEN valuenum END) AS monocytes_pct,
  AVG(CASE WHEN itemid = 100024 THEN valuenum END) AS eosinophils_pct,
  AVG(CASE WHEN itemid = 100025 THEN valuenum END) AS basophils_pct,
  -- === 白细胞绝对计数 ===
  AVG(CASE WHEN itemid = 100038 THEN valuenum END) AS abs_neutrophil,
  AVG(CASE WHEN itemid = 100037 THEN valuenum END) AS abs_lymphocyte,
  AVG(CASE WHEN itemid = 100039 THEN valuenum END) AS abs_monocyte,
  AVG(CASE WHEN itemid = 100036 THEN valuenum END) AS abs_eosinophil,
  AVG(CASE WHEN itemid = 100035 THEN valuenum END) AS abs_basophil,
  -- === 凝血 ===
  AVG(CASE WHEN itemid = 100034 THEN valuenum END) AS pt,
  AVG(CASE WHEN itemid = 100030 THEN valuenum END) AS inr,
  AVG(CASE WHEN itemid = 100046 THEN valuenum END) AS ptt,
  AVG(CASE WHEN itemid = 100048 THEN valuenum END) AS fibrinogen
FROM lab_window
GROUP BY subject_id, hadm_id, stay_id, icu_intime
ORDER BY subject_id, stay_id;


-- ============================================================
-- PART B: 动态指标 Day 1~7 (-6h ~ +168h)，每天平均值（长表）
-- ============================================================
DROP TABLE IF EXISTS derived.all_icu_labs_daily;
CREATE TABLE derived.all_icu_labs_daily AS
WITH lab_timed AS (
  SELECT
    c.subject_id,
    c.hadm_id,
    c.stay_id,
    c.icu_intime,
    l.itemid,
    l.valuenum,
    l.charttime,
    EXTRACT(EPOCH FROM (l.charttime - c.icu_intime)) / 3600.0 AS hours_from_icu
  FROM derived.cohort_all_icu AS c
  JOIN hosp.labevents AS l ON l.subject_id = c.subject_id
  WHERE
    l.charttime BETWEEN (c.icu_intime - INTERVAL '6 hours')
                    AND (c.icu_intime + INTERVAL '168 hours')
    AND l.valuenum IS NOT NULL
    AND l.valuenum <> 9999999
    AND l.itemid IN (
      100055,  -- LDH
      100002,  -- Creatinine
      100014,  -- Platelet Count
      100039,  -- Absolute Monocyte Count
      100038,  -- Absolute Neutrophil Count
      100037   -- Absolute Lymphocyte Count
    )
),
lab_with_day AS (
  SELECT
    subject_id,
    hadm_id,
    stay_id,
    icu_intime,
    itemid,
    valuenum,
    hours_from_icu,
    -- 划分 ICU 天数（d1: -6~24h, d2: 24~48h, ... d7: 144~168h）
    CASE
      WHEN hours_from_icu >= -6   AND hours_from_icu < 24  THEN 'd1'
      WHEN hours_from_icu >= 24   AND hours_from_icu < 48  THEN 'd2'
      WHEN hours_from_icu >= 48   AND hours_from_icu < 72  THEN 'd3'
      WHEN hours_from_icu >= 72   AND hours_from_icu < 96  THEN 'd4'
      WHEN hours_from_icu >= 96   AND hours_from_icu < 120 THEN 'd5'
      WHEN hours_from_icu >= 120  AND hours_from_icu < 144 THEN 'd6'
      WHEN hours_from_icu >= 144  AND hours_from_icu < 168 THEN 'd7'
    END AS icu_day
  FROM lab_timed
)
SELECT
  subject_id,
  hadm_id,
  stay_id,
  icu_intime,
  icu_day,
  -- 动态指标各天平均值
  AVG(CASE WHEN itemid = 100055 THEN valuenum END) AS ldh,
  AVG(CASE WHEN itemid = 100002 THEN valuenum END) AS creatinine,
  AVG(CASE WHEN itemid = 100014 THEN valuenum END) AS platelet,
  AVG(CASE WHEN itemid = 100039 THEN valuenum END) AS abs_monocyte,
  AVG(CASE WHEN itemid = 100038 THEN valuenum END) AS abs_neutrophil,
  AVG(CASE WHEN itemid = 100037 THEN valuenum END) AS abs_lymphocyte
FROM lab_with_day
WHERE icu_day IS NOT NULL
GROUP BY subject_id, hadm_id, stay_id, icu_intime, icu_day
ORDER BY subject_id, stay_id, icu_day;
