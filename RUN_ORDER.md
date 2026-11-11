# Reproduce the analysis and figures

## 1. Prepare licensed source databases

Use credentialed MIMIC-IV 3.1, MIMIC-III 1.4 and NWICU 0.1.0 sources. MIMIC concept-table prerequisites are referenced by the extraction SQL and available from https://github.com/MIT-LCP/mimic-code. Database schemas and local ports are defined in `00a_CR_PREPARE_PATIENT_INPUTS.py` (5432 for MIMIC-IV/NWICU; 5433 for MIMIC-III). Adapt these explicit local settings to your installation. Use your normal PostgreSQL credential mechanism; do not put credentials in code.

For an NWICU working database, execute the following included SQL in order. These commands rebuild named `derived` objects and should be run in a dedicated authorized analysis database. The `mock_*` tables are empty compatibility structures; no synthetic events are inserted.

```bash
export IR_REVISION_ROOT="$(pwd)"
export IR_CR_STAGE="$(mktemp -d)"
export NW_SQL="$IR_REVISION_ROOT/02_PA_主分析_AdmissionOrigin_M5P/03_NWICU外部验证/02_代码"
for sql in 01_create_derived_schema.sql 02_build_sofa_dependencies.sql   critical_nwicu_part1_cohort_all_icu.sql 07_all_icu_labs.sql   08_charlson.sql 09_sirs.sql 10_creatinine_baseline.sql   11_kdigo_stages.sql 12_icu_treatment.sql 13_icu_fluid_daily.sql   critical_nwicu_part2_lab_extract_all_icu.sql; do
  psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d nwicu -f "$NW_SQL/upstream/$sql"
done
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d nwicu -f "$NW_SQL/00_NW_adapted_SOFA.sql"
```

## 2. Export restricted analysis inputs

Keep all patient records in the ignored `.restricted_input` directory. Extraction uses exact ICU and post-ICU-admission hospital duration of at least 24 hours. Ranking the complete ICU table precedes eligibility filters.

```bash
umask 077
mkdir -p "$IR_REVISION_ROOT/.restricted_input"
export INPUT="$IR_REVISION_ROOT/.restricted_input"
export PA="$IR_REVISION_ROOT/02_PA_主分析_AdmissionOrigin_M5P"
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d mimiciv --csv  -f "$PA/01_MIMICIV/02_代码/critical mimiciv3.1_local .sql" > "$INPUT/mimiciv_analysis.csv"
# MIMIC-IV exports first/closest values inside mimiciv_analysis.csv.
# Optional event-level export for provenance:
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d mimiciv --csv  -f "$PA/01_MIMICIV/02_代码/02_PA_M4_AO_ALL_exposure_lab_events.sql" > "$INPUT/mimiciv_exposure_lab_events.csv"
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5433 -d mimiciii  -v mimiciii_output_csv="$INPUT/mimiciii_external_validation.csv"  -v mimiciii_lab_events_csv="$INPUT/mimiciii_exposure_lab_events.csv"  -v mimiciii_lab_definitions_csv="$INPUT/mimiciii_lab_timing_definitions.csv"  -f "$PA/02_MIMICIII外部验证/02_代码/01_PA_M3_AO_M5P_ALL_extract.sql"
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d nwicu  -v nwicu_output_csv="$INPUT/nwicu_external_validation.csv"  -v nwicu_lab_events_csv="$INPUT/nwicu_exposure_lab_events.csv"  -v nwicu_lab_definitions_csv="$INPUT/nwicu_lab_timing_definitions.csv"  -f "$NW_SQL/01_PA_NW_AO_M5P_ALL_extract.sql"
```

## 3. Run the current statistical workflow

Use the software versions recorded in `REPRODUCIBILITY/*session.txt`. Python requires pandas, numpy, matplotlib and python-docx; R dependencies are listed in the README. Output aggregate paths are fixed under the package root; figure exports go under `$IR_CR_STAGE/figures`.

```bash
export CR="$IR_REVISION_ROOT/08_CR_竞争风险修订/02_代码"
python "$CR/00a_CR_PREPARE_PATIENT_INPUTS.py"
python "$CR/00b_CR_DURATION_AUDIT.py"
for script in 00_CR_EXPORT_DESCRIPTIVES.R 02_CR_RUN_ANALYSES.R   03_CR_RUN_SENSITIVITIES.R 03b_CR_RUN_SUPPORTING_MODELS.R   03c_CR_EVENT_TIME_AUDIT.R 03d_CR_REVIEW_ANALYSES.R; do
  Rscript "$CR/$script"
done
Rscript "$IR_REVISION_ROOT/00_总说明与运行顺序/21_RUN_PREDICTION_AND_SUBGROUPS.R"
Rscript "$CR/03e_CR_PREDICTION_BOOTSTRAP.R"
python "$CR/04_CR_BUILD_FIGURES.py"
```

For a flowchart-only update use `python "$CR/04b_CR_BUILD_FLOWCHART.py"`; it checks all 15 sequential count differences and writes the figure plus both source CSVs. For a subgroup-only update use `Rscript "$CR/02b_CR_RUN_SUBGROUPS.R"` followed by `python "$CR/04a_CR_BUILD_SUBGROUP_FIGURE.py"`.

## 4. Validate the distributed snapshot

Run `python VALIDATE_PUBLIC.py` before recomputation. The checksum manifest identifies the distributed snapshot; a fresh run naturally changes aggregate output files. Validation checks complete regression matrices, denominators, effect labels, bootstrap completion, figure inventory, aggregate-only file fields and transfer integrity. Patient-level reruns require licensed data access and are not performed by this validator.
