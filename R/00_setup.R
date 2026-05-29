# =============================================================================
# 00_setup.R
# Project: Reproduction of Dunvald et al. (2022) Clinical PK Tutorial
# Purpose: Install all required packages and initialise renv environment
# Author:  Pramod BR
# Date:    2026-05-25
#
# Reference:
#   Dunvald A-CD, et al. Tutorial: Statistical analysis and reporting
#   of clinical pharmacokinetic studies.
#   Clin Transl Sci. 2022;15:1856-1866. doi:10.1111/cts.13305
# =============================================================================

# ── WHY renv? ─────────────────────────────────────────────────────────────────
# renv creates a project-local package library and a lockfile (renv.lock).
# Anyone who clones this repo runs renv::restore() and gets the exact same
# package versions you used. This is what "reproducible" means in practice.
# Without it, a collaborator on a newer R version may get different results.
# ─────────────────────────────────────────────────────────────────────────────

# Step 1: Bootstrap renv (only needed once per project)

# This installs renv itself and creates the project library scaffold
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::init()  

# Creates renv/ folder and renv.lock. Say YES when prompted.

# ── WHY each package? ─────────────────────────────────────────────────────────

# ncappc:
#   The exact NCA package used in Dunvald 2022 (Appendix S1, Section II).
#   Uses the linear-up/log-down trapezoidal method to calculate AUClast,
#   AUCinf, Cmax, t½, CL/F, Tmax. Output is NONMEM-compatible.
#   Reference: Acharya C et al. Comput Methods Programs Biomed. 2016;127:83-93

# PKNCA:
#   Modern alternative NCA package. We use it to CROSS-VALIDATE results from
#   ncappc (same numbers = confidence in both). Output is SDTM-ready for FDA
#   submissions. Used by many pharma companies today.

# ggplot2:
#   Reproduce Figure 2 (mean ± SD concentration-time curve with log-scale
#   inset) and Figure 3 (spaghetti plot of individual CL/F). The grammar of
#   graphics makes publication-quality plots straightforward.

# dplyr + tidyr:
#   Data wrangling: reshaping wide → long (tidy) format, filtering BLQ values,
#   computing summary statistics per group and subject.

# tableone:
#   Generates Table 1 (demographic + clinical characteristics). Standard in
#   clinical research — produces median (IQR) for continuous variables and
#   n (%) for categorical, matching the paper's reporting style.

# pwr:
#   Sample size and power calculations (Appendix S1, Section VII).
#   Reproduces the "5 subjects required" calculation from the paper.

# scales:
#   Needed for log-scale axis formatting in ggplot2 (Figure 2 inset).

# knitr + kableExtra:
#   Renders clean Tables 1 and 2 inside the Quarto report with proper
#   formatting (bold headers, striped rows, footnotes).

# here:
#   Makes file paths relative to the project root, so the code works on any
#   machine regardless of where the repo is cloned.
#   here("data", "raw", "mock_data.csv") always resolves correctly.

# quarto (check only — Quarto is external software, not an R package):
#   The report is written in Quarto (.qmd). You need Quarto CLI installed
#   separately: https://quarto.org/docs/get-started/

# ─────────────────────────────────────────────────────────────────────────────

# Step 2: Install all project packages into the renv library

packages <- c(
  # Core PK analysis
  "ncappc",      # Published Paper's NCA package
  "PKNCA",       # Cross-validation NCA

  # Data manipulation
  "dplyr",      
  "tidyr",
  "readr",      # Fast CSV reading

  # Visualisation
  "ggplot2",
  "scales",      # Log axis labels
  "patchwork",   # Combine Fig 2 main + inset into one plot

  # Reporting
  "tableone",    # Table 1 demographics
  "knitr",
  "kableExtra",  # Styled tables in report

  # Sample size
  "pwr",

  # Utility
  "here"         # Portable file paths
)

install.packages(packages)

# Step 3: Snapshot the installed versions into renv.lock
# This is the file you commit to GitHub. Collaborators run renv::restore()
renv::snapshot()

# ── VERIFY everything loaded correctly ───────────────────────────────────────
cat("\n=== Package load check ===\n")
for (pkg in packages) {
  status <- tryCatch({
    library(pkg, character.only = TRUE)
    paste("OK  -", packageVersion(pkg))
  }, error = function(e) {
    paste("FAIL -", e$message)
  })
  cat(sprintf("%-12s : %s\n", pkg, status))
}

# ── SESSION INFO (save for reproducibility record) ────────────────────────────

# This records R version, OS, all loaded packages and their versions.
# Commit this output alongside your renv.lock.

sink("session_info.txt")
cat("Project: Dunvald 2022 Clinical PK Reproduction\n")
cat("Date:    ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
sessionInfo()
sink()

cat("\n✓ Setup complete. session_info.txt written.\n")
cat("✓ renv.lock created — commit this file to GitHub.\n")
cat("\nNext step: open R/01_load_data.R\n")
