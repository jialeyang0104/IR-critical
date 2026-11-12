# Reproduce the analysis and figures

## 1. Prepare licensed source databases

Use credentialed MIMIC-IV 3.1, MIMIC-III 1.4 and NWICU 0.1.0 sources. MIMIC concept-table prerequisites are referenced by the extraction SQL and available from https://github.com/MIT-LCP/mimic-code. Database schemas and local ports are defined in `01_prepare_inputs.py` (5432 for MIMIC-IV/NWICU; 5433 for MIMIC-III). Adapt these explicit local settings to your installation. Use your normal PostgreSQL credential mechanism; do not put credentials in code.

For an NWICU working database, execute the following included SQL in order. These commands rebuild named `derived` objects and should be run in a dedicated authorized analysis database. The `mock_*` tables are empty compatibility structures; no synthetic events are inserted.

```bash
export IR_REVISION_ROOT="$(pwd)"
export IR_CR_STAGE="$(mktemp -d)"
export NW_SQL="$IR_REVISION_ROOT/src/extraction/nwicu"
for sql in 01_create_schema.sql 02_sofa_dependencies.sql   03_cohort.sql 04_laboratory_events.sql   05_charlson.sql 06_sirs.sql 07_baseline_creatinine.sql   08_kdigo.sql 09_icu_treatment.sql 10_fluid_balance.sql   11_laboratory_summary.sql; do
  psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d nwicu -f "$NW_SQL/upstream/$sql"
done
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d nwicu -f "$NW_SQL/01_build_sofa.sql"
```

## 2. Export restricted analysis inputs

Keep all patient records in the ignored `data/restricted` directory. Extraction uses exact ICU and post-ICU-admission hospital duration of at least 24 hours. Ranking the complete ICU table precedes eligibility filters.

```bash
umask 077
mkdir -p "$IR_REVISION_ROOT/data/restricted"
export INPUT="$IR_REVISION_ROOT/data/restricted"
export PA="$IR_REVISION_ROOT/src/extraction"
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d mimiciv --csv  -f "$PA/mimiciv/01_extract_cohort.sql" > "$INPUT/mimiciv_analysis.csv"
# MIMIC-IV exports first/closest values inside mimiciv_analysis.csv.
# Optional event-level export for provenance:
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d mimiciv --csv  -f "$PA/mimiciv/02_extract_lab_events.sql" > "$INPUT/mimiciv_exposure_lab_events.csv"
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5433 -d mimiciii  -v mimiciii_output_csv="$INPUT/mimiciii_external_validation.csv"  -v mimiciii_lab_events_csv="$INPUT/mimiciii_exposure_lab_events.csv"  -v mimiciii_lab_definitions_csv="$INPUT/mimiciii_lab_timing_definitions.csv"  -f "$PA/mimiciii/01_extract_cohort.sql"
psql -X -v ON_ERROR_STOP=1 -h localhost -p 5432 -d nwicu  -v nwicu_output_csv="$INPUT/nwicu_external_validation.csv"  -v nwicu_lab_events_csv="$INPUT/nwicu_exposure_lab_events.csv"  -v nwicu_lab_definitions_csv="$INPUT/nwicu_lab_timing_definitions.csv"  -f "$NW_SQL/02_extract_cohort.sql"
```

## 3. Run the current statistical workflow

Use the software versions recorded in `reports/reproducibility/*session.txt`. Python requires pandas, numpy, matplotlib and python-docx; R dependencies are listed in the README. Output aggregate paths are fixed under the package root; figure exports go under `$IR_CR_STAGE/figures`.

```bash
export CR="$IR_REVISION_ROOT/src/analysis"
python "$CR/01_prepare_inputs.py"
python "$CR/02_audit_durations.py"
python "$IR_REVISION_ROOT/src/extraction/01_extract_baseline.py"
for script in 03_descriptives.R 04_primary_models.R   06_sensitivities.R 07_supporting_models.R   08_event_time_audit.R 09_review_analyses.R; do
  Rscript "$CR/$script"
done
Rscript "$CR/10_prediction.R"
Rscript "$CR/11_prediction_bootstrap.R"
python "$CR/12_build_figures.py"
```

For a flowchart-only update use `python "$CR/13_build_flowchart.py"`; it checks all 15 sequential count differences and writes the figure plus both source CSVs. For a subgroup-only update use `Rscript "$CR/05_subgroups.R"` followed by `python "$CR/14_build_subgroup_figure.py"`.

For the NWICU recorded 28-day mortality forest only, run `python "$CR/12_build_figures.py" --only nwicu-recorded-28d`. It produces Figure S2 with the display labels “Admission origin” and “24h landmark” from the unchanged aggregate estimates.

For the three-cohort tertile forest only, run `python "$CR/12_build_figures.py" --only primary-tertiles`. Figure S1 contains 144 T2/T1 and T3/T1 estimates across Models 1–3, distinct from the 72 per-unit estimates in Figure 2. Table S6 uses the same three primary regression CSVs, including per-unit, tertile contrasts, and separate ordinal trend models (288 rows). The exploratory NWICU recorded 28-day endpoint remains separate in Figure S2 and Table S7.

For a baseline-only update, run `python src/extraction/01_extract_baseline.py`, then `Rscript src/analysis/03_descriptives.R` and `python src/validation/03_validate_baseline_tables.py`. Reuse the restricted baseline supplements when the inputs and field definitions have not changed. The read-only extraction adds MIMIC-III and NWICU baseline supplements under `data/restricted/` and leaves the formal regression inputs unchanged. Table 1 compares the three final cohorts, Table S2 uses the same 42 variables, and Table S3 compares included and laboratory-excluded patients within each cohort. Table S4 compares MIMIC-IV hospital survivors (n=2,735) with non-survivors (n=321), including the overall cohort, all 42 baseline variables, and P values. Current sources are `00_CR_ALL_table1_observed_baseline.csv`, `00_CR_ALL_Table1_variable_missingness.csv`, `00_CR_ALL_included_excluded_comparison.csv`, and `00_CR_M4_survivor_nonsurvivor_baseline.csv`. The older `00_CR_M4_table1_observed_baseline.csv` is historical and is no longer a manuscript-table input.

## 4. Validate the distributed snapshot

Run `python src/validation/01_validate_public.py` before recomputation. The checksum manifest identifies the distributed snapshot; a fresh run naturally changes aggregate output files. Validation checks complete regression matrices, denominators, effect labels, bootstrap completion, figure inventory, aggregate-only file fields and transfer integrity. Patient-level reruns require licensed data access and are not performed by this validator.

After editing public files, run `python src/validation/02_update_manifest.py`, then `python src/validation/01_validate_public.py`. The manifest writer selects the Git public file set and excludes ignored private inputs.

Python packages are listed in `requirements.txt`. Original R versions are recorded in the session files; use R 4.5.3 with cmprsk, survival, Hmisc, mice and pROC. Shared helper files are imported selectively by current scripts; do not execute their historical standalone workflows.
