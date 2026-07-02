# ============================================================================
# CLOVERS biomarker-CATE flat file -- complete variable inventory
# ============================================================================
# This script contains NO executable logic. It is a documented manifest of
# every variable that needs to go into the flat file before any modeling
# (no imputation, no risk modeling, no CATE -- per Victor's stage 1).
#
# Source files referenced below:
#   yw.csv          -- Curated datasets/yw.csv
#   xvars_egdt.csv  -- Curated datasets/xvars_egdt.csv
#   xvars_other.csv -- Curated datasets/xvars_other.csv
#   Share.4.29.26.csv (or Clover_4_29_26_Editable.csv, identical) -- biomarkers
# ============================================================================


# ---------------------------------------------------------------------------
# SUBJECT ID
# ---------------------------------------------------------------------------
# id                  -- patient identifier, format "C-####"
#                        present in ALL FOUR source files; this is the merge
#                        key. NOTE: the original clovers_forHTE_flat.R script
#                        dropped this column from its final output -- it must
#                        be kept this time so biomarkers can be merged in and
#                        so patients can be traced later.


# ---------------------------------------------------------------------------
# OUTCOME(S)                                              [source: yw.csv]
# ---------------------------------------------------------------------------
# inhosp90            -- primary outcome: in-hospital death (or failure to
#                        survive to discharge) within 90 days. Binary 0/1.
# sofa_diff           -- secondary outcome: change in SOFA score, baseline
#                        to Day 3. Continuous. ~310/1563 missing (no Day-3
#                        assessment recorded -- death/discharge/missed visit
#                        before reassessment).
#
# also present in yw.csv but NOT currently used in the flat file -- confirm
# with Victor whether any of these should be added as additional outcomes:
#   osfd, vfd, rrtfd, vasofd, icufd, hfd   (various organ/ventilator/ICU/
#                                            hospital "free days" outcomes)


# ---------------------------------------------------------------------------
# TREATMENT ASSIGNMENT                                    [source: yw.csv]
# ---------------------------------------------------------------------------
# w                   -- randomized treatment arm, binary 0/1.
#                        (confirm with Victor/trial docs which value = which
#                        arm -- crystalloid-liberal vs early-vasopressor)


# ---------------------------------------------------------------------------
# CONTINUOUS CLINICAL COVARIATES ("norm.vars")
# ---------------------------------------------------------------------------
# These 24 variables are the continuous/near-continuous covariates from the
# original clovers_forHTE_flat.R script. Several are log- or sqrt-transformed
# versions of a raw lab value (raw version listed in parentheses).
#
#   age                 -- age in years                      [xvars_egdt.csv]
#   temp                -- temperature                        [xvars_egdt.csv]
#   rr                  -- respiratory rate                   [xvars_egdt.csv]
#   hr                  -- heart rate                         [xvars_egdt.csv]
#   map                 -- mean arterial pressure              [xvars_egdt.csv]
#   sbp                 -- systolic blood pressure             [xvars_egdt.csv]
#   sofa                -- baseline SOFA score (with GCS)      [xvars_egdt.csv]
#                          (verified: exact match to DERIVED.csv$d_sofa_gcs)
#   albumin             -- serum albumin                       [xvars_egdt.csv]
#   ln_bili             -- log(bilirubin); raw = "bili"        [xvars_egdt.csv]
#                          (one raw value of exactly 0 must be set to NA
#                          before taking the log)
#   cr                  -- creatinine                          [xvars_egdt.csv]
#   bun                 -- blood urea nitrogen                 [xvars_egdt.csv]
#   ln_g                -- log(glucose); raw = "g"             [xvars_egdt.csv]
#   ln_lac              -- log(lactate); raw = "lac"            [xvars_egdt.csv]
#   sqrt_plt            -- sqrt(platelets); raw = "plt"         [xvars_egdt.csv]
#   ln_wbc              -- log(white blood cell count); raw = "wbc"
#                          (raw values of exactly 0 must be set to NA first)
#                                                               [xvars_egdt.csv]
#   hgb                 -- hemoglobin                          [xvars_egdt.csv]
#   na                  -- sodium                               [xvars_other.csv]
#   bicarb              -- bicarbonate                          [xvars_other.csv]
#   prefluid            -- pre-randomization IV fluid volume    [xvars_other.csv]
#   ln_bmi              -- log(BMI); raw = "bmi"                [xvars_other.csv]
#   gcs                 -- Glasgow Coma Scale                  [xvars_egdt.csv]
#   charlson            -- Charlson comorbidity score (age-adjusted)
#                                                               [xvars_egdt.csv]
#                          KNOWN BUG: patients with age exactly 80 get 0
#                          age-points instead of the expected 4. Affects 18
#                          patients. Fix before use (add back 4 points where
#                          age == 80), or flag to Victor to regenerate.
#                          Also note: tumor/liver/diabetes/kidney severity
#                          grades from the raw charl_* fields are NOT counted
#                          in this score, even though the CRF collects them.
#   o2sat               -- oxygen saturation                   [xvars_egdt.csv]
#   s2f                 -- SpO2/FiO2 ratio                     [xvars_other.csv]


# ---------------------------------------------------------------------------
# BINARY CLINICAL COVARIATES ("bin.vars")
# ---------------------------------------------------------------------------
# These 11 variables are the binary covariates from the original script.
#
#   site_lung           -- infection source: lung               [xvars_egdt.csv]
#   site_abdom          -- infection source: abdomen             [xvars_egdt.csv]
#   site_urine          -- infection source: urinary             [xvars_egdt.csv]
#   mv                  -- mechanical ventilation at baseline    [xvars_egdt.csv]
#   vaso                -- vasopressor use at baseline           [xvars_egdt.csv]
#   ards                -- ARDS present                          [xvars_other.csv]
#   dial                -- dialysis-dependent                    [xvars_other.csv]
#   chf                 -- congestive heart failure               [xvars_other.csv]
#   copd                -- COPD                                   [xvars_other.csv]
#   liver               -- liver disease                          [xvars_other.csv]
#   kidney              -- kidney disease                         [xvars_other.csv]


# ---------------------------------------------------------------------------
# BIOMARKERS (V1 / baseline draw only)      [source: Share.4.29.26.csv]
# ---------------------------------------------------------------------------
# The biomarker file has up to 3 visits per patient (V1/V2/V3 = baseline,
# 24hr, 72hr per the CLOVERS Lab Manual). For a BASELINE covariate set,
# filter to visit == "V1" before merging.
#
# 1,340 of the 1,563 clinical-cohort patients have V1 biomarkers (the other
# 223 have no biomarker draw at all -- they will be dropped from any
# biomarker-inclusive model).
#
# All 8 are right-skewed; log-transform each before modeling, same
# convention as the clinical labs above (bili/g/lac/wbc/bmi).
#
#   il1_pg_ml           -- IL-1, pg/mL
#   ang1_pg_ml          -- Angiopoietin-1, pg/mL
#   ang2_pg_ml          -- Angiopoietin-2, pg/mL
#   tnfr_uplex          -- TNFR-1, pg/mL (U-Plex platform)
#   tnfr_calibrated_rplex -- TNFR-1, pg/mL (calibrated/R-Plex platform)
#                          NOTE: two different measurements of the same
#                          analyte on two platforms. Pick ONE for modeling --
#                          tnfr_calibrated_rplex was used in the merged
#                          dataset built so far, but confirm this choice
#                          with Victor before finalizing.
#   il6_pg_ml           -- IL-6, pg/mL
#   strem1_pg_ml        -- sTREM-1, pg/mL
#   kim1_pg_ml          -- KIM-1, pg/mL
#   srage_pg_ml         -- sRAGE, pg/mL
#
# LLOD/ULOD handling: already applied in this file (half-LLOD/half-ULOD
# imputation per Biomarker_LLOD_ULOD_2026_04_29.docx) -- do not re-impute
# detection-limit values; treat the values in this file as final.
#
# Missingness at V1 (within the 1,340-patient biomarker cohort) is low,
# all <=1%:
#   il1: 10/1340, ang1: 7/1340, ang2: 7/1340, tnfr_uplex: 7/1340,
#   tnfr_calibrated_rplex: 7/1340, il6: 7/1340, strem1: 7/1340,
#   kim1: 0/1340, srage: 11/1340


# ---------------------------------------------------------------------------
# VARIABLES SEEN IN SOURCE FILES BUT DELIBERATELY EXCLUDED FROM THE
# ORIGINAL FLAT FILE -- confirm with Victor whether any should be reinstated
# ---------------------------------------------------------------------------
#   pao2                -- excluded in original script: "too much missing"
#   covid               -- excluded in original script: "too much missing"
#   female              -- present in xvars_egdt.csv, NOT in original norm/bin
#                          lists. Likely an oversight -- ask Victor.
#   site_blood, site_skin, site_other -- these do NOT exist as separate
#                          columns in xvars_egdt.csv; only site_lung,
#                          site_abdom, site_urine are present. "Other/
#                          unknown source" is implicitly the reference
#                          category when all three flags are 0.
#   visit, sample        -- biomarker file housekeeping columns, not
#                          covariates -- exclude from the model


# ---------------------------------------------------------------------------
# TOTAL COUNT
# ---------------------------------------------------------------------------
#   1  id
#   2  outcomes        (inhosp90, sofa_diff)
#   1  treatment       (w)
#  24  continuous clinical covariates
#  11  binary clinical covariates
#   8  biomarkers (V1, pick one TNFR-1 column -> 8, or keep both -> 9)
# ---
#  47  total columns (clinical-only model uses 1+2+1+24+11 = 39 columns,
#                     i.e. 35 covariates; adding biomarkers brings it to
#                     43 covariates)