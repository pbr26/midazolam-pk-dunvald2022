# =============================================================================
# 04b_nca_pknca.R
# Project : Reproduction of Dunvald et al. (2022)
# Purpose : Cross-validate NCA results using PKNCA package
#           Compare output to ncappc (04_nca_ncappc.R)
# Author  : [Your Name]
# Date    : 2026-05-26
#
# WHAT THIS SCRIPT DOES:
#   The paper uses ncappc. We ALSO run PKNCA — a modern NCA package used
#   in industry and regulatory submissions — and compare results.
#   If both packages give the same AUClast, Cmax, t½ → confidence is high.
#   This cross-validation is YOUR original contribution beyond the paper.
#
# WHY PKNCA?
#   PKNCA output is formatted for SDTM (Study Data Tabulation Model) —
#   the standard required for FDA/EMA regulatory submissions.
#   Learning PKNCA today means your skills are industry-ready.
# =============================================================================

library(PKNCA)   # Modern NCA package — SDTM-compatible
library(dplyr)
library(tidyr)
library(readr)
library(here)

# MY UNDERSTANDING: ___________________________________________________________

# ── Load data ─────────────────────────────────────────────────────────────────
long_data <- read_csv(here("data","processed","mock_data_long.csv"),
                      show_col_types=FALSE)

# Load ncappc results for comparison
nca_ncappc <- read_csv(here("output","tables","nca_parameters_individual.csv"),
                       show_col_types=FALSE)

# ── Prepare PKNCA input ───────────────────────────────────────────────────────
# PKNCA needs two separate objects:
#   PKNCAconc  — concentration data
#   PKNCAdose  — dose data
# Then combined into a PKNCAdata object

conc_data <- long_data %>%
  filter(BLQ == 0 | TIME == 0) %>%
  mutate(SID = paste0(ID, "_", ARM)) %>%
  select(SID, TIME, CONC, ARM, ID)

dose_data <- conc_data %>%
  distinct(SID, ARM, ID) %>%
  mutate(DOSE = 2, TIME = 0)   # Single 2mg dose at t=0

# PKNCAconc: tell PKNCA which columns are concentration, time, and subject
conc_obj <- PKNCAconc(
  data     = conc_data,
  formula  = CONC ~ TIME | SID,
  # formula reads: CONC as a function of TIME, grouped by SID
  # | SID means: treat each SID as a separate PK profile
)

# PKNCAdose: tell PKNCA about dosing
dose_obj <- PKNCAdose(
  data    = dose_data,
  formula = DOSE ~ TIME | SID
)

# Combine into PKNCAdata
pk_data <- PKNCAdata(conc_obj, dose_obj)

# ── Specify which parameters to calculate ────────────────────────────────────
# By default PKNCA calculates everything — we restrict to what we need
intervals_df <- data.frame(
  start      = 0,
  end        = 12,         # Our sampling window
  auclast    = TRUE,       # AUC to last measurable concentration
  cmax       = TRUE,       # Maximum concentration
  tmax       = TRUE,       # Time of Cmax
  half.life  = TRUE,       # Terminal half-life
  cl.obs     = TRUE        # Apparent clearance (CL/F for oral)
)

# Apply intervals to all subjects
pk_data$intervals <- merge(
  pk_data$intervals[, c("SID")],
  intervals_df
)

# ── Run PKNCA ─────────────────────────────────────────────────────────────────
cat("Running PKNCA...\n")
pknca_result <- pk.nca(pk_data)
cat("✓ PKNCA complete\n")

# ── Extract results ───────────────────────────────────────────────────────────
pknca_table <- as.data.frame(pknca_result$result) %>%
  filter(PPTESTCD %in% c("auclast","cmax","tmax","half.life","cl.obs")) %>%
  select(SID, Parameter = PPTESTCD, Value = PPORRES) %>%
  pivot_wider(names_from = Parameter, values_from = Value) %>%
  separate(SID, into = c("ID","ARM"), sep = "_", extra = "merge") %>%
  mutate(
    ID  = as.integer(ID),
    cl.obs = cl.obs * 1000   # Convert to L/h
  ) %>%
  rename(
    AUClast_pknca = auclast,
    Cmax_pknca    = cmax,
    Tmax_pknca    = tmax,
    t_half_pknca  = half.life,
    CLF_pknca     = cl.obs
  ) %>%
  arrange(ARM, ID)

# ── Compare ncappc vs PKNCA ───────────────────────────────────────────────────
comparison <- nca_ncappc %>%
  rename(
    AUClast_ncappc = AUClast,
    Cmax_ncappc    = Cmax,
    t_half_ncappc  = t_half,
    CLF_ncappc     = CLF
  ) %>%
  left_join(pknca_table, by = c("ID","ARM")) %>%
  mutate(
    AUC_diff_pct  = round(abs(AUClast_ncappc - AUClast_pknca) / AUClast_ncappc * 100, 1),
    Cmax_diff_pct = round(abs(Cmax_ncappc - Cmax_pknca) / Cmax_ncappc * 100, 1),
    t_half_diff   = round(abs(t_half_ncappc - t_half_pknca), 3)
  )

cat("\n=== ncappc vs PKNCA comparison (first 6 rows) ===\n")
print(head(comparison %>% select(ID, ARM,
                                  AUClast_ncappc, AUClast_pknca, AUC_diff_pct,
                                  Cmax_ncappc, Cmax_pknca, Cmax_diff_pct), 6))

cat("\n── AGREEMENT SUMMARY ───────────────────────────────────────────────────\n")
cat("Mean AUClast difference:", round(mean(comparison$AUC_diff_pct, na.rm=TRUE),2), "%\n")
cat("Mean Cmax difference:   ", round(mean(comparison$Cmax_diff_pct, na.rm=TRUE),2), "%\n")
cat("Mean t½ difference:     ", round(mean(comparison$t_half_diff, na.rm=TRUE),4), "h\n")
cat("Acceptable threshold    : <5% difference for NCA cross-validation\n")
cat("────────────────────────────────────────────────────────────────────────\n")

# Save
write_csv(comparison, here("output","tables","nca_crossvalidation.csv"))
cat("\n✓ Saved: output/tables/nca_crossvalidation.csv\n")
cat("\n=== 04b_nca_pknca.R COMPLETE ===\n")
cat("Next: commit to GitHub then open R/05_gmr_testing.R\n")
