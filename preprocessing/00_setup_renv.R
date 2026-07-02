# ============================================================================
# 00_setup_renv.R
#
# PURPOSE: Run this ONCE, yourself, after installing all required packages.
# It creates renv.lock -- R's equivalent of a lockfile (like uv.lock for
# Python) -- recording the exact version of every package this project
# uses. Anyone else (Victor, a future you) can then run renv::restore() to
# get the identical package versions, instead of "whatever's on CRAN today".
#
# WORKFLOW:
#   1. You run this script once, after install.packages()-ing everything
#      the project needs (see the list below).
#   2. It creates an renv.lock file in the project root. Commit that file
#      to GitHub along with your scripts.
#   3. Whoever clones the repo (Victor, etc.) runs renv::restore() ONE TIME
#      -- not this script -- and gets the exact same package versions you
#      used. They do NOT need to re-run this setup script.
#
# This only needs to be re-run if you add or upgrade a package later.
# ============================================================================

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# ---- make sure every package this project actually uses is installed ----
# (renv::snapshot() only records what's installed AND used in the project's
# .R files -- so install everything first, then snapshot.)

required_packages <- c(
  "readxl",        # reading .xlsx files (01b_explore_all_datasets.R)
  "grf",           # causal_forest (cf.CATE, in cate-repo/scoring.methods_v4.R)
  "stochtree",     # BART-based methods (pbart.slearn.CATE, in cate-repo/)
  "glmnet",        # elastic net / LASSO methods (cate-repo/)
  "randomForest",  # missforest_train.R and rf-based CATE methods
  "caret",         # used inside several scoring methods (cate-repo/)
  "recipes",       # used inside several scoring methods (cate-repo/)
  "miselect",      # penalized variable selection under MI (cate-repo/)
  "ic.infer",       # GATES monotonicity test (gates.inference.R)
  "scales",        # plotting helper (gates.inference.R)
  "sandwich",      # robust standard errors (gates.inference.R)
  "bcf",           # Bayesian causal forest (bcf.CATE, if used)
  "nnet"           # only needed if you end up using a nnet.CATE.* method
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:\n")
  print(missing_packages)
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
} else {
  cat("All required packages are already installed.\n")
}

# ---- initialize renv for this project (only does anything the first time) ----
source(file.path(dirname(getwd()), "config", "config.R"))

if (!file.exists(file.path(PROJECT_ROOT, "renv.lock"))) {
  cat("Initializing renv for this project...\n")
  renv::init(project = PROJECT_ROOT, bare = TRUE)
} else {
  cat("renv already initialized for this project.\n")
}

# ---- snapshot: write the exact installed versions to renv.lock ----
renv::snapshot(project = PROJECT_ROOT, prompt = FALSE)

cat("\n============================================================\n")
cat("Done. renv.lock has been written to:\n  ", file.path(PROJECT_ROOT, "renv.lock"), "\n")
cat("Commit this file (and the renv/ folder it created) to your GitHub repo.\n")
cat("Anyone else can then run renv::restore() ONCE to get these exact\n")
cat("package versions -- they do not need to run this script.\n")
cat("============================================================\n")