
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
## Average treatment effects may mask true treatment differences

The CLOVERS trial (Shapiro et al., NEJM 2023) randomized 1,563 patients with sepsis-induced hypotension to either a restrictive fluid strategy (prioritize early vasopressors) or a liberal fluid strategy (prioritize crystalloid boluses).

Although the primary result was null: no statistically significant difference in 90-day in-hospital mortality between the two treatment across the full trial population, there was a more nuanced finding. 

A null average treatment effect (ATE) does NOT mean the treatment had no effect on anyone. It means the average effect across all patients was approximately zero. 

This average can hide real heterogeneity: some patients could have benefited from the liberal strategy, while others were actually adversed affected by it.

These competing effects may be masked when aggregated. 

Thus, this is the primary aim of this analysis: to explore the implications of heterogeneous treatment effects (HTE). 

Specifically, we wish to examined individualized absolute risk difference (iARD). For each patient, the iARD is the difference between their probability of the outcome (90-day in-hospital mortality) under treatment A versus treatment B. 

In a trial with a binary outcome, a patient's iARD represents how much their personal risk changes depending on which treatment they receive. The ATE reported in the CLOVERS paper is simply the average of all patients' iARDs. When the ATE is null, the question becomes: are the individual iARDs also all near zero (true homogeneity), or do some patients have strongly positive iARDs and others strongly negative iARDs that average to zero (hidden heterogeneity)?

 Data-driven methods for estimating individual-level treatment effects such as conditional average treatment effects (CATEs) address limitations by using the full covariate profile of each patient to estimate their personal iARD, without pre-specifying which covariates drive the heterogeneity.

For context, Kiernan et al. (AJRCCM 2025) demonstrated that this hidden heterogeneity exists in CLOVERS. Using latent class analysis (LCA) with 3 plasma biomarkers (angiopoietin-1, angiopoietin-2, and sTNFR-1), they classified 1,289 patients into two molecular subphenotypes (SP1 and SP2). SP2 patients, characterized by markers of endothelial injury and inflammation, had significantly higher 28-day mortality with the liberal fluid strategy compared to the restrictive strategy (41% vs 27%), while SP1 patients showed no difference (9% vs 9%). The interaction p-value was 0.02.

The pipeline in this repository is looking at the next step which is if CATE methods detect heterogeneity without pre-specifying subphenotypes. It also takes a look at MORE biomarkers (expanding from the original 3 markers to 8 markers). How do these additions improve the detection of the individualized differences in critically ill patients? 

---

## This pipeline's cohort 
Our analytic cohort (n=1,340) is defined as CLOVERS patients with any Visit 1 biomarker draw. 

This differs from Kiernan et al.'s cohort (n=1,289), which required complete measurements of all three LCA biomarkers (Ang-1, Ang-2, sTNFR-1) plus serum creatinine. Our broader inclusion is defensible because causal forest methods can handle missing individual biomarkers via imputation, whereas LCA requires all clustering variables to be complete.

From 1,563 randomized → 223 excluded (no V1 biomarker) → 1,340 analytic cohort (672 restrictive / 668 liberal so 186 events).

---
# Glossary

ATE (Average Treatment Effect): The mean difference in outcome between treatment arms across the entire study population. A null ATE does not rule out real heterogeneity.

HTE (Heterogeneous Treatment Effects): Non-random, explainable variability in the magnitude or direction of treatment effects across individuals within a population.

iARD (individualized Absolute Risk Difference): The estimated difference in outcome probability for a specific patient under treatment A versus treatment B. The ATE is the average of all patients' iARDs.

CATE (Conditional Average Treatment Effect): The expected treatment effect for patients with a given set of covariate values. CATE models estimate iARDs using each patient's full baseline profile — demographics, labs, vitals, comorbidities, and (in this project) biomarkers.

TOC (Targeting Operator Characteristic): A curve that plots the observed absolute risk difference as a function of the proportion of the population treated, starting from those predicted to benefit most. If a CATE model has detected real heterogeneity, the curve will separate from the ATE line.

AUTOC (Area Under the TOC Curve): Summary metric for how well a CATE model separates patients with different treatment effects. Higher AUTOC = more detected heterogeneity.

AUQINI (Area Under the Qini Curve): A linear-weighted version of AUTOC; emphasizes separation across the full score distribution.

Kendall's tau: Rank correlation between the CATE model's predicted treatment effect ordering and the observed ordering. Values near 1 mean the model correctly ranks who benefits most.

AUC (Area Under the ROC Curve): In the risk-modeling stage, measures how well a baseline risk model discriminates patients who will have the outcome from those who will not. AUC 0.5 = chance; AUC 1.0 = perfect discrimination.

LCA (Latent Class Analysis): An unsupervised clustering method that assigns patients to latent subgroups based on observed variables. Used by Kiernan et al. (2025) to define SP1/SP2 subphenotypes. Unlike CATE, LCA requires pre-specifying which variables define the clusters and produces hard (not continuous) assignments.

Causal forest: A random-forest-based method that estimates CATEs by adaptively splitting on covariates to find regions of differing treatment effects. The primary CATE method in this analysis.

BCF (Bayesian Causal Forest): Bayesian nonparametric CATE method that separately models baseline risk and the treatment effect modifier, using BART priors.

BART: A Bayesian nonparametric regression method that sums many small decision trees, each regularized by a prior. Used here both as a risk model (bart.risk) and as a CATE S-learner (pbart.slearn.CATE).

missForest:  A random-forest-based imputation method that iteratively predicts each missing variable using all other variables. Used here with a custom implementation that allows training on one dataset and applying to another.



## Background / literature

Please reference the iARD and HTE work completed by Dr. Victor Talisa and Dr. Faraaz Shah. 

Papers referenced for this pipeline include the titles: 
- Early Restrictive or Liberal Fluid Management for Sepsis-Induced Hypotension

- Early Restrictive Versus Liberal Fluid Management for Sepsis-induced HypotensionOnline Supplement 

- Molecular Phenotyping of Sepsis and Differential Responseto Fluid Resuscitation (Elizabeth Kiernan et al)

  - The original CLOVERS trial paper (Shapiro NI et al., NEJM 2023)

  - Prior HTE work in sepsis using causal forests / similar methods
    (e.g., Seymour et al. JAMA 2019  Sinha et al.
    work, Wong et al.).


  - Published ADRENAL HTE analysis (V. Talisa) 

---

## What was done, and what was changed from the ADRENAL pipeline
This pipeline was built to mirror the HTE analysis that Dr. Victor Talisa developed for the ADRENAL trial. 

Victor's framework is a five-stage pipeline for discovering and evaluating heterogeneous treatment effects in randomized trial data. The shared codebase can be found in 
/ cate-repo/ (CATE-derivation-and-evaluation-main) and consists of 9 R scripts. 

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

Missing covariates are imputed via missForest using the `missforest_train.R` implementation from the shared CATE repo.
Imputation forests are fit on the derivation set only. So, the same fitted
forests are applied to the validation set, preventing leakage of validation
information.

### Derivation / validation split

A 50/50 random split, stratified by treatment arm (`w`), with a fixed
random seed (`set.seed(03202026)`) for reproducibility. Y

ields 670/670
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


### Results 

---

## Results

### Risk modeling (stage 4)

Five baseline risk-prediction methods were compared on the derivation set
using 50-fold cross-validation, fitting on control-arm patients only.

| Method | Mean AUC | SD | Median AUC |
|---|---|---|---|
| Elastic net (alpha=0.5) | 0.862 | 0.017 | 0.864 |
| Random forest (randomForest) | 0.860 | 0.015 | 0.861 |
| BART (stochtree) | 0.855 | 0.020 | 0.856 |
| Random forest (grf) | 0.851 | 0.018 | 0.851 |
| Logistic regression | 0.815 | 0.026 | 0.816 |

The top four methods are essentially tied (AUC 0.851–0.862). Logistic
regression is clearly worst but still above 0.80. Elastic net performs
marginally best, consistent with moderate-dimensional clinical data where
some regularization helps. XGBoost was dropped due to R 4.2 incompatibility.
Total runtime: 27.3 minutes.

### CATE modeling (stage 5) — the headline result

The CATE comparison was run twice on the derivation set (n=670) with
cv_k=100 and 5 methods. Pass 1 used 35 clinical covariates only.
Pass 2 used 35 clinical + 8 log-transformed biomarkers (43 total).
Both passes used the same seed, same CV splits, same methods — only
the covariate vector differed.

**AUTOC comparison (×100; higher = more detected heterogeneity):**

| Method | Clinical only | + Biomarker | Change |
|---|---|---|---|
| Causal Forest | 2.77 | 2.80 | +0.03 (flat) |
| Local Linear RF R-learner | 0.87 | 2.09 | +1.22 (large) |
| Bayesian Causal Forest | 1.72 | 2.03 | +0.31 (meaningful) |
| Tian Elastic Net | 0.79 | 1.08 | +0.29 (meaningful) |
| Probit BART S-learner | -0.92 | -1.44 | inverted (method failure) |

**Kendall's tau (rank correlation; closer to 1 = better ranking):**

| Method | Clinical only | + Biomarker |
|---|---|---|
| Causal Forest | 0.96 | 0.91 |
| Local Linear RF R-learner | 0.46 | 0.94 |
| Bayesian Causal Forest | 0.93 | 0.95 |
| Tian Elastic Net | 0.55 | 0.60 |
| Probit BART S-learner | -0.60 | -0.59 |

**Interpretation:**

Heterogeneous treatment effects exist in CLOVERS — four of five CATE
methods detect positive AUTOC in both passes. This is consistent with
Kiernan et al.'s finding that molecular subphenotypes respond differently
to fluid resuscitation strategy.

Biomarkers carry real HTE-relevant information. The local linear RF
R-learner more than doubles its AUTOC when biomarkers are added
(0.87 → 2.09), and its Kendall's tau jumps from 0.46 to 0.94 —
meaning it goes from barely ranking patients correctly to near-perfect
ranking. Bayesian causal forest and Tian elastic net also show meaningful
improvement.

However, the causal forest — the best-performing method overall — barely
improves with biomarkers (2.77 → 2.80). It already finds most of the
heterogeneity from clinical variables alone. This suggests that a
sufficiently flexible nonparametric method can extract HTE-relevant
signal from clinical covariates that overlaps substantially with what
the biomarkers provide.

Probit BART produces negative AUTOC in both passes, meaning its predicted
treatment effect rankings are inverted relative to the truth. This method
is not functioning on this dataset and should be excluded from
interpretation.

Runtime: Pass 1 took 6.9 minutes, Pass 2 took 9.8 minutes.

### Variable importance (from confirmatory model)

When cf.CATE was trained on the full derivation set and applied to
validation with the +biomarker covariate set, the top 10 variables
by importance were:

| Rank | Variable | Importance |
|---|---|---|
| 1 | Temperature | 0.064 |
| 2 | sRAGE (biomarker) | 0.057 |
| 3 | Angiopoietin-2 (biomarker) | 0.056 |
| 4 | S/F ratio | 0.055 |
| 5 | sTNFR-1 (biomarker) | 0.052 |
| 6 | Sodium | 0.042 |
| 7 | sTREM-1 (biomarker) | 0.039 |
| 8 | Lactate (log) | 0.037 |
| 9 | Age | 0.037 |
| 10 | Bicarbonate | 0.035 |

Four of the top 10 variables are biomarkers (sRAGE, Ang-2, sTNFR-1,
sTREM-1). Notably, three of these four (Ang-2, sTNFR-1, and by
extension sRAGE as an endothelial/inflammatory marker) overlap with the
biomarkers Kiernan et al. used for their LCA subphenotyping, providing
convergent evidence that endothelial injury and inflammation markers drive
treatment effect heterogeneity in CLOVERS.

In the clinical-only model, the top variables were S/F ratio, temperature,
BUN, glucose (log), sodium, lactate (log), age, hemoglobin, bicarbonate,
and platelets (sqrt).

### Confirmatory analysis (stage 6)

The best CATE method (cf.CATE) was trained on the full derivation set
(n=670) and used to predict individualized treatment effects on the
held-out validation set (n=670). Validation patients were classified
into benefit, indeterminate, and harm tertiles based on their predicted
CATE scores.

**Clinical only model:**

| Subgroup | N | Mortality (liberal) | Mortality (restrictive) | Difference | p-value |
|---|---|---|---|---|---|
| Benefit (bottom tertile) | 223 | 7.5% | 6.9% | +0.6% | 0.50 |
| Indeterminate | 224 | 9.7% | 10.8% | -1.1% | 0.96 |
| Harm (top tertile) | 223 | 19.3% | 20.2% | -0.9% | 0.50 |

**Clinical + biomarker model:**

| Subgroup | N | Mortality (liberal) | Mortality (restrictive) | Difference | p-value |
|---|---|---|---|---|---|
| Benefit (bottom tertile) | 223 | 5.6% | 4.3% | +1.3% | 0.55 |
| Indeterminate | 224 | 11.1% | 9.4% | +1.8% | 0.83 |
| Harm (top tertile) | 223 | 20.0% | 23.9% | -3.9% | 0.71 |

No subgroup reached statistical significance in either model (all p > 0.50).
However, the models do correctly stratify patients by overall risk: the
harm tertile has approximately 20% mortality compared to 5-7% in the
benefit tertile, confirming that the CATE model captures meaningful
prognostic variation even if the subgroup-specific treatment effects
are not individually significant.

The null confirmatory result is consistent with underpowering rather than
absence of HTE. With approximately 93 events in the 670-patient validation
set split into thirds, each subgroup contains roughly 31 events and
approximately 110 patients per treatment arm — substantially less power
than Kiernan et al.'s interaction test, which used 1,289 patients in two
groups (not three). Future work with larger validation samples or
alternative confirmation strategies (e.g., benefit vs everyone else rather
than tertiles) may recover the signal.

### Confirmatory analysis

No subgroup reached significance (all p > 0.50). The model 
stratifies by risk (harm group ~20% mortality vs benefit group ~5-7%).
There are only ~93 events in 670 validation patients split into thirds so this may be attributed to underpowering. 
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



### Run order

# From the project root:

# Stage 0: setup
Rscript preprocessing/00_setup_renv.R

# Stage 1: data assembly
cd preprocessing
Rscript 01b_explore_all_datasets.R
Rscript 02_build_flat_file.R
Rscript 02b_build_variable_source_table.R
Rscript 02c_missing.R
Rscript 00_sofa_charlson_audit.R

# Stage 2: split
Rscript 03_split_derivation.R
Rscript 03b_comparasion.R

# Stage 3: imputation
cd ../scripts
Rscript 04_Impute.R

# Stage 4: risk modeling
Rscript 05_risk_modeling.R

# Stage 5: CATE modeling (SLOW: several hours)
Rscript 07_cate_modeling.R

# Stage 6: confirmatory analysis
Rscript 10_confirmatory_analysis.R

# Stage 7: Publication Materials
Rscript 11_publication_tables.R

# Supplementary
cd ../preprocessing
Rscript ../flowcharts/

See mermaid diagrams rendered related to the archecture, methodology, and participants. 

---

## Repo structure

clover_hte_biomarkers/
├── config/
│   └── config.R                    <- all paths; edit PROJECT_ROOT once
├── data/                           <- raw data (gitignored)
│   ├── Data/
│   │   ├── Curated datasets/       <- yw.csv, xvars_egdt.csv, xvars_other.csv
│   │   └── data/csv/               <- DATASET.csv, DERIVED.csv
│   └── Share.4.29.26.csv           <- biomarker file (V1+V2+V3)
├── cate-repo/                      <- Victor's shared CATE toolkit (unchanged)
├── preprocessing/                  <- stages 0-2 scripts
├── scripts/                        <- stages 3-6 and post-CATE output scripts
├── outputs/                        <- each script writes to its own subfolder
│   ├── 02_build_flat_file/
│   ├── 03_split_derivation_validation/
│   ├── 04_impute/
│   ├── 05_risk_modeling/
│   ├── 07_cate_modeling/
│   ├── 09_build_cate_outputs/
│   └── 10_confirmatory_analysis/
|.  └── 11_table1/
    
├── flowcharts/                     <- Mermaid participant flow diagrams
├── documents/                      <- eTables builder
├── relevant_literature/            <- reference PDFs (these are .gitignored in public repos)
└── README.md


## Flowcharts
https://mermaid.ai/open-source/intro/index.html?utm_source=mermaid_js&utm_medium=landing_pop_up&utm_campaign=docs_v1

Flowcharts were rendered using Mermaid. See the above documentation for more information.

## Contact

Contact Jaspreet Singh, Victor Talisa, or Faraaz Shah with questions on this analysis. 



# R Studio

Installation of new R packages 
Rscript -e 'install.packages("randomForest", repos="https://cloud.r-project.org")' 2>&1 | tail -20

Use of packages, including Random Forest requires installation of fortran (as random forest is written in fortran based code similiarly gboost has a Fortran/C build)
Rscript -e 'install.packages(c("glmnet","xgboost"), repos="https://cloud.r-project.org")'
Rscript -e 'library(stochtree); cat("stochtree OK\n")'
Rscript -e 'install.packages("grf", repos="https://cloud.r-project.org")'