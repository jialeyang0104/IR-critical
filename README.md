# IR-critical

## Manuscript

**Prognostic performance of the stress hyperglycemia ratio, triglyceride-glucose index, and atherogenic index of plasma in critically ill patients: a multidatabase cohort study**

This study compares three laboratory-derived metabolic indices—the stress hyperglycemia ratio (SHR), triglyceride-glucose (TyG) index, and atherogenic index of plasma (AIP)—for mortality associations and prognostic discrimination in critically ill adults. MIMIC-IV is the primary cohort; MIMIC-III and the Northwestern ICU (NWICU) database provide replication cohorts.

| Database | Version | Common-index cohort | Role |
|---|---|---:|---|
| MIMIC-IV | 3.1 | 3,056 | Primary analysis and internal prediction evaluation |
| MIMIC-III | 1.4 | 1,133 | Historical, same-institution replication |
| NWICU | 0.1.0 | 1,556 | Geographically distinct replication |

The first recorded ICU stay is selected for each patient before applying eligibility criteria. Eligible patients are adults with ICU duration and post-ICU-admission hospital follow-up of at least 24 hours and valid glucose, HbA1c, triglyceride, and HDL-C measurements within −6 to +24 hours of ICU admission. All three indices are compared within a common patient sample.

### Analysis overview

- **Primary outcome:** in-hospital mortality. Fine–Gray regression treats discharge alive as the competing event.
- **Secondary outcomes:** first-ICU mortality, using Fine–Gray regression with ICU departure alive as the competing event; and 28-day mortality, using Cox regression in MIMIC. NWICU 28-day recorded deaths are exploratory because ascertainment is limited to the hospital system.
- **Adjustment:** Model 1 is unadjusted; Model 2 adjusts for age and sex; Model 3 additionally adjusts for SOFA, Charlson comorbidity index, and diabetes.
- **Supporting analyses:** restricted cubic splines, cumulative incidence and Kaplan–Meier curves, ROC curves and paired DeLong comparisons, severity-score extensions, repeated cross-validation, and patient-bootstrap uncertainty.
- **Sensitivity analyses:** alternative laboratory summaries, a 24-hour landmark origin, selected comorbidity exclusions, and exploratory subgroup analyses.

Formal adjustment variables are complete in the final cohorts, so the current regression analyses do not use multiple imputation.

Table 1 presents the observed baseline characteristics of all three final cohorts, with between-cohort P values. Supplementary Table S2 reports missingness for exactly those 42 variables, and Table S3 provides three panels of included-versus-laboratory-excluded comparisons with P values. Table S4 reports the same 42 variables for the overall MIMIC-IV cohort and by in-hospital survival status (2,735 survivors and 321 non-survivors), with between-group P values. These tables use observed values without imputation. The three-cohort sources have the `00_CR_ALL_` prefix in [`results/analysis/`](results/analysis/); the survival-status source is [`00_CR_M4_survivor_nonsurvivor_baseline.csv`](results/analysis/00_CR_M4_survivor_nonsurvivor_baseline.csv). Source availability and field mappings are supplied alongside them.

SHR was positively associated with mortality across the three cohorts. In MIMIC-IV, it showed higher standalone mortality discrimination than TyG and AIP and increased cross-validated AUC when added to established severity scores. Exact estimates and uncertainty are available in the [aggregate results](results/analysis/) and [figure source data](deliverables/figures/source_data/).

## Repository contents

| Path | Contents |
|---|---|
| [`src/extraction/`](src/extraction/) | Cohort extraction SQL and NWICU prerequisite construction |
| [`src/analysis/`](src/analysis/) | Ordered preparation, regression, sensitivity, prediction, and figure entrypoints |
| [`src/shared/`](src/shared/) | Shared statistical and prediction functions |
| [`src/figures/`](src/figures/) | Shared plotting utilities |
| [`src/validation/`](src/validation/) | Aggregate-result and public-file validation |
| [`results/analysis/`](results/analysis/) | Regression estimates, cohort flows, risk sets, and supporting analyses |
| [`results/prediction/`](results/prediction/) | ROC, DeLong, cross-validation, and other prediction summaries |
| [`deliverables/figures/`](deliverables/figures/) | Six main and ten supplementary figure groups, exported as PDF, SVG, and PNG |
| [`reports/reproducibility/`](reports/reproducibility/) | Software sessions, public-file inventory, checksums, and path migration |
| [`docs/RUN_ORDER.md`](docs/RUN_ORDER.md) | Full extraction and execution commands |
| [`docs/INPUT_DICTIONARY.md`](docs/INPUT_DICTIONARY.md) | Restricted input fields, units, endpoint definitions, and database prerequisites |

The [table–result–figure mapping](results/analysis/10_CR_table_result_figure_mapping.csv) identifies the sources for the 16 supplementary tables. Some additional aggregates document historical or supporting calculations; the mapping distinguishes the current manuscript tables.

Figure 2 presents per-unit associations across the three cohorts and Models 1–3. Supplementary Figure S1 presents cohort-specific tertile contrasts (T2/T1 and T3/T1), rather than duplicating those continuous-index estimates. Table S6 reports both scales and separate ordinal trend models. Baseline continuous-variable comparisons use Kruskal–Wallis tests across three cohorts or Mann–Whitney U tests when only two cohorts have observed data; categorical comparisons use Pearson chi-square or Fisher exact tests when any expected cell count is below five.

Supplementary figures follow their first main-text citation: S1, three-cohort tertile associations; S2, exploratory NWICU recorded 28-day mortality; S3, 24-hour landmark; S4–S6, alternative laboratory timing; S7–S9, comorbidity exclusions; S10, subgroup associations.

## Reproduce the article

### 1. Obtain the repository and software

```bash
git clone https://github.com/jialeyang0104/IR-critical.git
cd IR-critical
python -m pip install -r requirements.txt
```

The recorded statistical environment used **R 4.5.3** with `cmprsk`, `survival`, `Hmisc`, `mice`, and `pROC`; exact versions are in [`reports/reproducibility/`](reports/reproducibility/). `mice` remains a shared-code dependency but is not used to impute the formal regression covariates. Database extraction requires the PostgreSQL client and licensed source databases.

### 2. Rebuild figures from the supplied aggregate results

This route requires no patient-level data:

```bash
export IR_REVISION_ROOT="$(pwd)"
export IR_CR_STAGE="$(mktemp -d)"
python src/validation/01_validate_public.py
python src/analysis/12_build_figures.py
printf 'Generated figures: %s/figures\n' "$IR_CR_STAGE"
```

The figure builder produces all 16 figure groups and their source CSVs in the temporary output directory. The committed figure exports remain unchanged. For a flowchart-only update, run `python src/analysis/13_build_flowchart.py` with the same environment variables.

### 3. Reproduce the patient-level analyses

Obtain authorized access to the source resources:

- [MIMIC-IV 3.1](https://doi.org/10.13026/kpb9-mt58)
- [MIMIC-III 1.4](https://doi.org/10.13026/C2XW26)
- [NWICU 0.1.0](https://doi.org/10.13026/s84w-1829)

Follow [RUN_ORDER.md](docs/RUN_ORDER.md), sections 1–2, to construct prerequisites and export the required inputs into `data/restricted/`. The [input dictionary](docs/INPUT_DICTIONARY.md) specifies fields, units, cleaning ranges, laboratory windows, and outcome coding. Configure database names and ports for your installation; use PostgreSQL's normal credential mechanism.

After extraction, run from the repository root:

```bash
export IR_REVISION_ROOT="$(pwd)"
export IR_CR_STAGE="$(mktemp -d)"
export CR="$IR_REVISION_ROOT/src/analysis"

python "$CR/01_prepare_inputs.py"
python "$CR/02_audit_durations.py"
python "$IR_REVISION_ROOT/src/extraction/01_extract_baseline.py"

for script in \
  03_descriptives.R \
  04_primary_models.R \
  06_sensitivities.R \
  07_supporting_models.R \
  08_event_time_audit.R \
  09_review_analyses.R \
  10_prediction.R \
  11_prediction_bootstrap.R; do
  Rscript "$CR/$script"
done

python "$CR/12_build_figures.py"
```

The primary-model entrypoint also invokes `05_subgroups.R`. The prediction bootstrap uses 200 patient resamples with repeated 10×5-fold cross-validation. Statistical scripts update aggregate files under `results/`; generated figures are written to `$IR_CR_STAGE/figures`. Use a separate clone when recomputing results while retaining the distributed snapshot.

### 4. Check the outputs

The public validator checks regression matrices, outcome-specific denominators, effect labels, subgroup and bootstrap completeness, the 16-figure inventory, excluded patient identifiers, and snapshot checksums:

```bash
python src/validation/01_validate_public.py
```

Run this command before recomputation to verify the distributed snapshot. After intentionally updating public outputs, refresh the inventory and checksums, then validate again:

```bash
python src/validation/02_update_manifest.py
python src/validation/01_validate_public.py
```

## Data and manuscript availability

The public repository distributes code, non-identifying aggregate results, and figure sources. Patient-level inputs remain subject to the respective database access agreements and are excluded from Git.
