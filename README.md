
# CLOVERS HTE / Biomarker CATE Analysis

## What this project does

This repository contains the  pipeline for evaluating whether 8 differnt 
biomarkers found in plasma improve detection of heterogeneous treatment effects (HTE) in the
CLOVERS sepsis trial patient population, beyond clinical covariates alone.

The pipeline mirrors the analytic structure of a prior HTE analysis on the
ADRENAL trial (V. Talisa et al.), adapting it to CLOVERS-specific data and
extending it to incorporate 8 plasma biomarkers (IL-1, Angiopoietin-1/2,
TNFR-1, IL-6, sTREM-1, KIM-1, sRAGE) measured at trial baseline.

---

## Central aims

1. **Aim 1 (descriptive).** Characterize baseline clinical and biomarker
   profiles of CLOVERS patients with available baseline (V1) biomarker draws,
   compared to the full CLOVERS cohort.

2. **Aim 2 (HTE detection — the headline question).** Determine whether
   adding plasma biomarkers to a clinical-covariate-based heterogeneous
   treatment effect (HTE) model meaningfully improves detection of HTE,
   versus a model using clinical covariates alone.

3. **Aim 3 (subgroup characterization).** If Aim 2 is positive, characterize
   the patient subgroup identified by the biomarker-augmented model as
   benefiting more (or less) from the experimental treatment arm, in
   clinical terms.

---

## Background / literature

Please see attached relevant literature folder. Reference the iARD and HTE work completed by Dr. Victor Talisa and Dr. Faraaz Shah. 

  - The original CLOVERS trial paper (Shapiro NI et al., NEJM 2023): aim,
    primary result, lack of treatment effect on the overall cohort.

  - The CLOVERS protocol and biomarker substudy protocol if separately
    published, for the rationale behind which biomarkers were drawn.

  - Prior HTE work in sepsis using causal forests / similar methods
    (e.g., Seymour et al. JAMA 2019 phenotypes, Sinha et al. clustering
    work, Wong et al. PERSEVERE).

  - Methodological references for the CATE estimators used here: Athey &
    Wager 2018 (grf), Hahn et al. 2020 (bcf), Tian et al. 2014 (modified
    covariate), and the AUTOC/AUQINI metrics from Yadlowsky et al. 2021.

  - The published ADRENAL HTE analysis (V. Talisa) 

---

## What was done, and what was changed from the ADRENAL pipeline

This pipeline follows the same five-stage structure as the ADRENAL HTE
analysis but adapts each stage to CLOVERS data and adds biomarker handling.

### Pipeline structure

1. Build a flat analytic file. I built a structured dataset (no imputation, no modeling)
2. Random 50/50 derivation/validation split, stratified by treatment arm (fluid/limited fluid)
3. missForest-based imputation (fit on derivation, applied to validation)
4. Risk modeling on imputed data
5. CATE modeling: 5 methods × cross-validation × TOC/AUTOC comparison

### Changes from ADRENAL

| Change | Reason |
|---|---|
| Two outcomes carried (`inhosp90`, `sofa_diff`) instead of single `y` | CLOVERS has both 90-day in-hospital mortality and a SOFA-change outcome |
| 8 baseline biomarkers added as covariates (log-transformed) | Project's central aim is to test biomarker contribution |
| CATE comparison run **twice** (clinical-only vs +biomarker) | Direct comparison is the headline result |
| Risk-modeling suite expanded from 3 methods to 5 | Adds elastic net (`enet.risk`) and BART (`bart.risk`) to Victor's existing logistic regression + 2 random forest variants |
| Charlson age-80 patch applied | Original `xvars_egdt.csv` has a bug giving 0 age-points to patients aged exactly 80 (18 patients affected); a corrected version is computed |
| All paths driven by `config/config.R` | Portability across machines (the ADRENAL scripts hardcoded Windows paths) |
| Each script writes to its own subfolder under `outputs/` | Output traceability — no two scripts can silently overwrite each other |

### Known issues within the analysis pipeline


- The Charlson scores for pts age-80+ : 18 patients get 0 age-points instead of 4. tried to implement a fix. need senior collaborator input on resolution. 
- Charlson omits tumor / liver / diabetes / kidney severity grades, even
  though the CRF collects them. 
- Two TNFR-1 measurements exist in the biomarker file (`tnfr_uplex` and
  `tnfr_calibrated_rplex`); we use the calibrated rPlex version. 

- `missforest_train.R` runs `randomForest` in regression mode for binary
  covariates, leaving 23 patients with fractional imputed values (e.g.,
  `kidney = 0.27`). Decision was to leave as-is rather than round, since
  causal forest can accept continuous indicators.
- `ADRENAL_internal_impute.R` writes `dsi_tr.RDS` + `dsi_te.RDS` but
  `ADRENAL_internal_derivation_v2.R` reads a single `dsi0.RDS` — there's
  a missing rbind step in Victor's chain. We build the combined object
  explicitly here in `04_impute.R`.

---

## Statistical methodology

### Cohort

The biomarker-CATE analytic cohort is the intersection of (a) randomized
CLOVERS participants with complete baseline clinical covariates and (b)
participants with a baseline (V1) biomarker draw. 

In total, n=1,340 from the 1,563 originally randomized patients
n=223 are excluded for absence of a V1 (visit 1) biomarker sample.
Treatment arm balance is preserved (672/668 in the
biomarker cohort vs 782/781 in the full cohort).

### Outcomes

Primary outcome: 90-day in-hospital mortality (`inhosp90`, binary).
Secondary outcome: SOFA score change from baseline to Day 3 (`sofa_diff`,
continuous), with 18.3% missing in the biomarker cohort due to absent Day-3
assessments (death, discharge, or missed visit).

### Covariates

35 clinical covariates from the CLOVERS curated datasets (`yw.csv`,
`xvars_egdt.csv`, `xvars_other.csv`), and 8 plasma biomarkers from
`Share.4.29.26.csv` (V1 draw only). Continuous skewed variables (bilirubin,
glucose, lactate, platelets, WBC, BMI, all biomarkers) are log- or sqrt-
transformed. Charlson Comorbidity Index is age-bracket-adjusted; the
upstream-data version has a known bug at age 80 that is corrected here.
Detailed variable-source mapping in `outputs/02b_build_variable_source_table/`.

### Missing data handling

Missing covariates are imputed via missForest (Stekhoven & Bühlmann 2012)
using the `missforest_train.R` implementation from the shared CATE repo.
Imputation forests are fit on the derivation set only; the same fitted
forests are applied to the validation set, preventing leakage of validation
information.

### Derivation / validation split

A 50/50 random split, stratified by treatment arm (`w`), with a fixed
random seed (`set.seed(03202026)`) for reproducibility. Yields 670/670
patients. Statistical equivalence of the two halves on every covariate
is verified in `outputs/03b_compare_der_val/`.

### Risk modeling

Five risk-prediction methods are compared on derivation, using out-of-sample
AUC across 50 repeated 50/50 cross-validation splits:

| Method | Implementation |
|---|---|
| Logistic regression | `glm(family=binomial)` |
| Random forest (regression forest) | `grf::regression_forest` |
| Random forest (Foster et al. 2011) | `randomForest` |
| Elastic net penalized logistic regression | `glmnet::cv.glmnet` (alpha=0.5) |
| BART (Bayesian additive regression trees) | `stochtree::bart` |

All five fit on control-arm patients only, predicting baseline (untreated)
risk, matching the convention in the shared CATE repo. (XGBoost was
intended as a sixth method but was unavailable for the installed R version;
omitted from the final comparison.)

### CATE / HTE estimation

Five CATE-scoring methods are compared, identical to those used in the
ADRENAL HTE analysis:

| Method | Implementation |
|---|---|
| Causal forest | `grf::causal_forest` via `cf.CATE` |
| Local linear random forest (R-learner) | `llrf.CATE.rlearn` |
| Bayesian causal forest | `bcf` via `bcf.CATE` |
| Probit BART (S-learner) | `stochtree::bart` via `pbart.slearn.CATE` |
| Modified-covariate elastic net (Tian) | `glmnet` via `ModCovElasticNet.CATE` |

Each method runs through 100 repeated 50/50 cross-validation splits on
derivation. Per-fold treatment effect estimates are post-processed via
`post.cv.fun` to compute AUTOC (area under the targeting operator
characteristic), AUQINI (area under the Qini curve), Kendall's tau, and
adjusted AUQINI metrics. The full comparison is run twice — once with
clinical covariates only, once with clinical + 8 biomarker covariates —
producing a side-by-side AUTOC table that is the project's primary result.

### Confirmatory analysis

[TODO: describe `confirm.analysis.R` usage on the validation set once it's
written.]

---

## How to reproduce

### One-time setup

```bash
# 1. Clone the repo
git clone <REPO_URL>
cd clover_hte_biomarkers

# 2. Set PROJECT_ROOT in config/config.R to the absolute path of this folder

# 3. Place the raw CLOVERS data under data/ (see below)

# 4. Place the shared CATE-derivation-and-evaluation-main repo under cate-repo/

# 5. Install R package dependencies (or use renv::restore() if renv.lock exists):
Rscript preprocessing/00_setup_renv.R
```

### Data layout expected

```
data/
├── Data/
│   ├── Curated datasets/
│   │   ├── yw.csv
│   │   ├── xvars_egdt.csv
│   │   └── xvars_other.csv
│   └── data/csv/
│       ├── DATASET.csv
│       └── DERIVED.csv
└── Share.4.29.26.csv     <- the biomarker file, directly in data/
```

### Run order

```bash
cd preprocessing
Rscript 01_check_data_availability.R    # diagnostic: do all files exist + load?
Rscript 02_build_flat_file.R            # stage 1: the flat file
Rscript 02b_build_variable_source_table.R  # variable → source CSV
Rscript 02c_check_missingness_and_qc.R  # missingness summary
Rscript 03_split_derivation_validation.R   # stage 2: 50/50 split
Rscript 03b_compare_der_val.R           # stage 2b: equivalence check
Rscript 04_impute.R                     # stage 3: missForest imputation
Rscript 05_risk_modeling.R              # stage 4: risk modeling comparison
Rscript 06_participant_flow_extract.R   # CONSORT-style counts + Mermaid
Rscript 07_cate_modeling.R              # stage 5: CATE comparison (LONG: hours)
Rscript 00_sofa_charlson_audit.R        # standalone: SOFA/Charlson audit
```

Each script writes to its own subfolder under `outputs/`. Total runtime
for the full pipeline is dominated by `07_cate_modeling.R` (several hours,
possibly overnight). All other scripts together run in well under one hour.

---

## Repo structure

```
clover_hte_biomarkers/
├── config/config.R              <- all paths, edit PROJECT_ROOT once
├── data/                        <- raw inputs
├── cate-repo/                   <- shared CATE-derivation-and-evaluation-main
├── preprocessing/               <- numbered analytic scripts
├── outputs/                     <- each script writes to its own subfolder
└── docs/                        <- README, supplementary tables, methods
```

---

## Contact

Contact Jaspreet Singh, Victor Talisa, or Faraaz Shah with questions on this analysis. 



# R Studio

Installation of new R packages 
Rscript -e 'install.packages("randomForest", repos="https://cloud.r-project.org")' 2>&1 | tail -20

Use of packages, including Random Forest requires installation of fortran (as random forest is written in fortran based code similiarly gboost has a Fortran/C build)
Rscript -e 'install.packages(c("glmnet","xgboost"), repos="https://cloud.r-project.org")'
Rscript -e 'library(stochtree); cat("stochtree OK\n")'
Rscript -e 'install.packages("grf", repos="https://cloud.r-project.org")'