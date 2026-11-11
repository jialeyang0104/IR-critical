# Metabolic indices and ICU mortality: patient-level competing-risk analysis

This package contains the code, aggregate results and figure sources for the current manuscript revision. The unit is one patient per database: select the earliest recorded ICU stay from the complete ICU table **before** applying eligibility criteria. The resulting common-index cohorts contain 3,056 MIMIC-IV, 1,133 MIMIC-III and 1,556 NWICU patients. MIMIC-III deidentified ages ≥300 are coded as 90 years (65 affected final-cohort patients).

In-hospital and first-ICU mortality use Fine–Gray models with alive departure as the competing event; MIMIC 28-day mortality uses Cox regression. NWICU recorded 28-day deaths are exploratory because ascertainment is limited to Northwestern Medicine hospital data. All five formal covariates are observed. Compatibility filenames retain M5/M20; m=0 and maxit=0 in result columns identify the actual complete-data analysis.

## Contents and evidence

|Path|Purpose|
|---|---|
|`08_CR_竞争风险修订/02_代码/`|Current statistical analysis and figure-generation code|
|`08_CR_竞争风险修订/03_汇总结果/`|Aggregate estimates, flows, diagnostics, bootstrap contrasts, table-source mapping|
|`08_CR_竞争风险修订/01_最终投稿文件/figures/`|16 current figure groups (6 main, 10 supplementary) (PDF/SVG/PNG), source CSVs, export manifest|
|`00_总说明与运行顺序/`|Shared preparation, descriptive, timing and prediction functions|
|`02_PA_主分析_AdmissionOrigin_M5P/`|Current PostgreSQL extraction scripts and adapted NWICU SOFA construction|
|`03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/`|Current ROC/DeLong, calibration, DCA and prediction outputs|
|`REPRODUCIBILITY/`|Executed software versions, input hashes and verified citation metadata|

`INPUT_DICTIONARY.md` defines restricted inputs and database prerequisites. No patient-level CSV, model object or credential is distributed. Code and result filenames are retained to preserve executable relative paths; this README and the input dictionary are the English navigation authorities.

## Run

Use R 4.5.3 with cmprsk 2.2-12, survival 3.8-12, Hmisc 5.3-0, pROC 1.19.1 and mice 3.19.0 (shared function dependency; not used for formal covariate imputation). Python requires pandas, numpy, matplotlib, Pillow, lxml and python-docx. Install a PostgreSQL client. Set `IR_REVISION_ROOT` to the extracted package root and `IR_CR_STAGE` to an empty temporary output directory. `IR_PSQL` optionally selects the psql executable. Set connection names/ports in the explicitly declared database specifications when reproducing on another host.

1. Obtain credentialed source data under the respective agreements and prepare the raw eligible extracts using the cohort SQL. Follow `RUN_ORDER.md` for exact extraction commands and the ordered NWICU derived-table build. MIMIC derived tables follow the official MIT-LCP concept definitions. NWICU prerequisite construction and the adapted SOFA SQL are included. These prerequisite scripts replace named derived objects; execute them in an authorized working database.
2. Run Python `08_CR_竞争风险修订/02_代码/00a_CR_PREPARE_PATIENT_INPUTS.py` and `00b_CR_DURATION_AUDIT.py`. These query complete source tables read-only, rank patients before eligibility, enrich metadata, and write restricted derivative inputs plus aggregate audits. The original raw extracts are preserved.
3. Run R scripts `00_CR_EXPORT_DESCRIPTIVES.R`, `02_CR_RUN_ANALYSES.R`, `03_CR_RUN_SENSITIVITIES.R`, `03b_CR_RUN_SUPPORTING_MODELS.R`, `03c_CR_EVENT_TIME_AUDIT.R` and `03d_CR_REVIEW_ANALYSES.R`, in that order, under `08_CR_竞争风险修订/02_代码/`.
4. Run R `00_总说明与运行顺序/21_RUN_PREDICTION_AND_SUBGROUPS.R`, then `08_CR_竞争风险修订/02_代码/03e_CR_PREDICTION_BOOTSTRAP.R`. The latter refits paired 10×5-fold CV in 200 patient bootstrap samples and supplies the current increment CIs and paired contrasts. The current subgroup/interaction entry is `02b_CR_RUN_SUBGROUPS.R`, also invoked by the main CR analysis. It fits 54 common-cohort-SD subgroup effects and 27 index-by-age/sex/diabetes interactions with two-sided Wald tests and Benjamini–Hochberg adjustment. Shared scripts 19/20 are function sources for the current pipeline; their standalone historical Cox workflows are not the manuscript's primary analysis.
5. Run Python `08_CR_竞争风险修订/02_代码/04_CR_BUILD_FIGURES.py`. For a subgroup-only update, run `02b_CR_RUN_SUBGROUPS.R` then `04a_CR_BUILD_SUBGROUP_FIGURE.py`. Every figure, including flow, ROC and incremental AUC, is rebuilt from current aggregates under `$IR_CR_STAGE/figures`. Table mapping is `03_汇总结果/10_CR_table_result_figure_mapping.csv`.
6. Run `python VALIDATE_PUBLIC.py` for matrix, denominator, figure, file-field and checksum verification. Three-cohort defaults are required for a complete package; single-cohort execution replaces combined summaries.

The statistical and figure workflow is self-contained apart from the credentialed database inputs and installed software. Manuscript-generation, reviewer-response, private Zotero mapping, and Word integration files are kept in the local author workspace. The distributed reference metadata contains bibliographic fields only.

## Definitions that affect interpretation

The ≥24-hour eligibility restriction selects patients who remain in ICU through the exposure window. The 24-hour-origin sensitivity is structurally redundant for stay-specific event ordering. Two MIMIC-IV day-level 28-day deaths are excluded by its positive-follow-up rule; risk sets are explicitly exported. Hospital discharge before recorded ICU departure explains the extra hospital-duration exclusions, reported in the duration audit.

Aggregate records retain previously computed analyses for reproducibility. The current manuscript omits the standalone standardized-effect comparisons, weighting figures, and cause-specific support figures requested for removal during author review. The manuscript retains four sensitivity categories: alternative laboratory summaries, a 24-hour landmark origin, selected comorbidity exclusions, and subgroup analysis. CareVue, COVID-19, glucocorticoid, and age-recoding sensitivities are retained only as historical computation records. The uncited proportionality-diagnostic and event-time tables are excluded from the supplement. The current table-source mapping identifies the 15 included supplementary tables, each cited in the main text and ordered by its first citation. Table 1 uses groupwise Shapiro–Wilk normality assessment to select mean (SD)/Welch t-test or median (IQR)/Mann–Whitney U. All 29 current continuous variables depart from normality in at least one group; the aggregate distribution audit records the actual decisions.

Prediction is internally evaluated in MIMIC-IV. Current prediction figures compare SOFA, SAPS II, APS III and OASIS alone with each index extension. Figure 1 uses three parallel patient-selection lanes with 15 sequential side exclusions. Figure 2 shows three-cohort associations; the MIMIC-IV three-model plot is Figure S1. There are six main and ten supplementary figures. Figure S10 and Table S15 report subgroup estimates and interaction tests using the full MIMIC-IV common-index cohort SD. TyG-by-diabetes interactions for in-hospital and 28-day mortality remain significant after Benjamini–Hochberg adjustment; no SHR or AIP interaction does. Calibration, DCA and Basic Model aggregates are retained as historical computation records and are not current manuscript tables or figures.

## Data resources and preservation

MIMIC-IV v3.1: https://doi.org/10.13026/kpb9-mt58; MIMIC-III v1.4: https://doi.org/10.13026/C2XW26; NWICU v0.1.0: https://doi.org/10.13026/s84w-1829. Follow resource-specific credentialing and use agreements. Local source-table counts are in the flow CSV and may differ from other snapshots.

`FILE_MANIFEST.csv` and `CHECKSUMS.sha256` describe this exact package. Restricted inputs stay outside distribution. Repository: https://github.com/jialeyang0104/IR-critical. No repository DOI has been assigned.
