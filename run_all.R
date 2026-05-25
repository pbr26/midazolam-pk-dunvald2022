# =============================================================================
# run_all.R
# Project: Reproduction of Dunvald et al. (2022)
# Purpose: Run the entire analysis pipeline from start to finish
#
# Usage:
#   Rscript run_all.R
#   OR in RStudio: source("run_all.R")
# =============================================================================

library(here)

cat("============================================================\n")
cat("  Dunvald 2022 Clinical PK Reproduction — Full Pipeline\n")
cat("============================================================\n\n")

# Helper to run each script with timing and error catching
run_script <- function(path, label) {
  cat(sprintf("Running: %s\n", label))
  cat(sprintf("  File: %s\n", path))
  start <- Sys.time()
  tryCatch({
    source(here(path))
    elapsed <- round(as.numeric(Sys.time() - start, units = "secs"), 1)
    cat(sprintf("  ✓ Done in %.1f s\n\n", elapsed))
  }, error = function(e) {
    cat(sprintf("  ✗ FAILED: %s\n\n", e$message))
    stop(sprintf("Pipeline halted at: %s", label))
  })
}

run_script("R/01_load_data.R",           "Step 1: Load and validate data")
run_script("R/02_descriptive_stats.R",   "Step 2: Descriptive statistics (Table 1)")
run_script("R/03_pk_plots.R",            "Step 3: PK plots (Figure 2 + Figure 3)")
run_script("R/04_nca_ncappc.R",          "Step 4: NCA with ncappc (Table 2)")
run_script("R/04b_nca_pknca.R",          "Step 4b: NCA cross-validation with PKNCA")
run_script("R/05_gmr_testing.R",         "Step 5: GMR, 95% CI, hypothesis tests")
run_script("R/06_sample_size.R",         "Step 6: Sample size calculation")

cat("============================================================\n")
cat("  All scripts complete.\n")
cat("  Now render the report:\n")
cat("    quarto render report/dunvald2022_reproduction.qmd\n")
cat("============================================================\n")
