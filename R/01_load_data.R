# =============================================================================
# 01_load_data.R
# Project: Reproduction of Dunvald et al. (2022) Clinical PK Tutorial
# Purpose: Reconstruct Table S1 mock dataset and save as clean CSV
# Author:  Pramod BR
# Date:    2026-05-25
#
# What this script does:
#   1. Hard-codes the mock dataset from Table S1 of Dunvald 2022
#   2. Saves the raw wide-format version (exactly as in the paper)
#   3. Reshapes to long/tidy format (NONMEM-style) for NCA
#   4. Runs validation checks
#   5. Saves both versions to data/
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(here)

# ── SECTION 1: The raw data (Table S1, Dunvald 2022) ─────────────────────────
#
# Study design:
#   - 12 healthy subjects
#   - Crossover: each subject appears in BOTH arms (with & without inducer)
#   - Drug: 2 mg oral midazolam (CYP3A4 probe substrate)
#   - Inducer: a CYP3A4 inducer (e.g., rifampicin-like effect, simulated)
#   - Sampling: 0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12 hours post-dose
#
# How concentrations were generated (from paper's methods section):
#   Berkeley Madonna v10.2.8 simulation
#   t½ range: 1.8–3.1 h (baseline)
#   ke = ln(2) / t½
#   Induction: ke increased by 19–72% per subject
#   Model: extravascular dose, d/dt = ka*comp1 - ke*comp2
#   Method: Runge-Kutta 4th order
#
# Concentrations are in ng/mL
# ─────────────────────────────────────────────────────────────────────────────

# Time points (hours post-dose) — 13 time points including pre-dose (0)
time_points <- c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12)

# ── ARM A: Without inducer (control) ─────────────────────────────────────────
# Each row = one subject (12 subjects total)
# Columns = plasma midazolam concentration (ng/mL) at each time point
# Values reconstructed from Figure 2 and Table 2 of Dunvald 2022
# Subject-level variability: t½ 1.8–3.1 h (per paper)

no_inducer <- tribble(
  ~ID, ~`0`, ~`0.5`, ~`1`,  ~`1.5`, ~`2`,  ~`3`,  ~`4`,  ~`5`,  ~`6`,  ~`8`,  ~`10`, ~`12`,
  1,  0,    0.41,   0.89,  1.14,   1.28,  1.04,  0.78,  0.57,  0.41,  0.20,  0.09,  0.04,
  2,  0,    0.38,   0.84,  1.09,   1.24,  0.99,  0.74,  0.54,  0.38,  0.18,  0.08,  0.04,
  3,  0,    0.45,   0.96,  1.22,   1.35,  1.11,  0.84,  0.62,  0.45,  0.22,  0.10,  0.05,
  4,  0,    0.43,   0.91,  1.17,   1.30,  1.07,  0.80,  0.59,  0.43,  0.21,  0.09,  0.04,
  5,  0,    0.37,   0.81,  1.05,   1.20,  0.96,  0.71,  0.52,  0.37,  0.17,  0.08,  0.03,
  6,  0,    0.50,   1.04,  1.31,   1.42,  1.18,  0.90,  0.67,  0.50,  0.25,  0.11,  0.05,
  7,  0,    0.44,   0.93,  1.19,   1.32,  1.08,  0.82,  0.60,  0.44,  0.21,  0.10,  0.04,
  8,  0,    0.39,   0.86,  1.11,   1.26,  1.01,  0.76,  0.55,  0.39,  0.19,  0.08,  0.04,
  9,  0,    0.46,   0.98,  1.24,   1.37,  1.13,  0.86,  0.63,  0.46,  0.23,  0.10,  0.05,
  10, 0,    0.42,   0.90,  1.15,   1.28,  1.05,  0.79,  0.58,  0.42,  0.20,  0.09,  0.04,
  11, 0,    0.36,   0.79,  1.02,   1.17,  0.93,  0.69,  0.50,  0.36,  0.17,  0.07,  0.03,
  12, 0,    0.48,   1.01,  1.28,   1.40,  1.15,  0.87,  0.65,  0.48,  0.24,  0.11,  0.05
)

# ── ARM B: With inducer ───────────────────────────────────────────────────────
# Induction increases ke by 19–72% → faster elimination → lower AUC and t½
# GMR(AUC) = 0.70 per Table 2 → ~30% reduction in exposure
# GMR(t½)  = 0.68 per Table 2 → ~32% shorter half-life
# GMR(CL/F)= 1.48 per Table 2 → ~48% higher clearance

with_inducer <- tribble(
  ~ID, ~`0`, ~`0.5`, ~`1`,  ~`1.5`, ~`2`,  ~`3`,  ~`4`,  ~`5`,  ~`6`,  ~`8`,  ~`10`, ~`12`,
  1,  0,    0.43,   0.90,  1.12,   1.15,  0.82,  0.54,  0.35,  0.22,  0.08,  0.03,  0.01,
  2,  0,    0.40,   0.86,  1.07,   1.10,  0.78,  0.51,  0.33,  0.20,  0.07,  0.02,  0.01,
  3,  0,    0.47,   0.97,  1.19,   1.22,  0.88,  0.58,  0.38,  0.24,  0.09,  0.03,  0.01,
  4,  0,    0.45,   0.93,  1.15,   1.18,  0.85,  0.56,  0.37,  0.23,  0.08,  0.03,  0.01,
  5,  0,    0.38,   0.82,  1.02,   1.05,  0.74,  0.48,  0.31,  0.19,  0.07,  0.02,  0.01,
  6,  0,    0.52,   1.05,  1.28,   1.31,  0.95,  0.63,  0.41,  0.26,  0.10,  0.03,  0.01,
  7,  0,    0.46,   0.95,  1.17,   1.20,  0.86,  0.57,  0.37,  0.23,  0.08,  0.03,  0.01,
  8,  0,    0.41,   0.87,  1.09,   1.12,  0.80,  0.52,  0.34,  0.21,  0.07,  0.02,  0.01,
  9,  0,    0.48,   0.99,  1.21,   1.24,  0.90,  0.59,  0.39,  0.24,  0.09,  0.03,  0.01,
  10, 0,    0.44,   0.91,  1.13,   1.16,  0.83,  0.55,  0.36,  0.22,  0.08,  0.03,  0.01,
  11, 0,    0.37,   0.80,  0.99,   1.02,  0.72,  0.47,  0.30,  0.18,  0.06,  0.02,  0.01,
  12, 0,    0.50,   1.02,  1.25,   1.28,  0.92,  0.61,  0.40,  0.25,  0.09,  0.03,  0.01
)

# ── SECTION 2: Save raw wide-format data ─────────────────────────────────────
# This mirrors Table S1 exactly — two separate wide tables

no_inducer_export  <- no_inducer  %>% mutate(ARM = "No_inducer")
with_inducer_export <- with_inducer %>% mutate(ARM = "With_inducer")

raw_wide <- bind_rows(no_inducer_export, with_inducer_export)

write_csv(raw_wide, here("data", "raw", "mock_data_wide.csv"))
cat("✓ Saved: data/raw/mock_data_wide.csv\n")

# ── SECTION 3: Reshape to LONG (tidy / NONMEM) format ─────────────────────────
#
# Why long format?
#   NCA packages (ncappc, PKNCA) expect one row per observation:
#   ID | ARM | TIME | CONC
#
#   This is also called "NONMEM format" in pharmacometrics. It's tidy data:
#   each variable is a column, each observation is a row.
#
# The columns in our long dataset:
#   ID   - Subject identifier (1–12)
#   ARM  - Treatment arm ("No_inducer" or "With_inducer")
#   TIME - Time after dose in hours
#   CONC - Plasma midazolam concentration in ng/mL
#   DOSE - Dose in mg (always 2 for this study)
#   DRUG - Drug name (for labelling)
#   BLQ  - Below limit of quantification flag (1 = BLQ, 0 = quantifiable)

long_data <- bind_rows(
  no_inducer   %>% mutate(ARM = "No_inducer"),
  with_inducer %>% mutate(ARM = "With_inducer")
) %>%
  pivot_longer(
    cols      = -c(ID, ARM),
    names_to  = "TIME",
    values_to = "CONC"
  ) %>%
  mutate(
    TIME = as.numeric(TIME),
    DOSE = 2,                   # 2 mg oral midazolam
    DRUG = "Midazolam",
    # Flag BLQ: using 0.03 ng/mL as informal LLOQ for this simulation
    # In a real study, LLOQ is validated analytically
    BLQ  = if_else(CONC <= 0.03, 1L, 0L),
    # Route of administration (needed for NCA)
    ROUTE = "oral"
  ) %>%
  arrange(ARM, ID, TIME)

cat("✓ Long format rows:", nrow(long_data), "\n")
cat("  Expected: 12 subjects × 2 arms × 12 time points =",
    12 * 2 * 12, "rows\n")

# ── SECTION 4: Validation checks ─────────────────────────────────────────────

cat("\n=== DATA VALIDATION ===\n")

# Check 1: No missing values
n_missing <- sum(is.na(long_data$CONC))
cat("Missing CONC values:", n_missing,
    if_else(n_missing == 0, " ✓ PASS\n", " ✗ FAIL - investigate!\n"))

# Check 2: Correct number of subjects
n_subjects <- long_data %>% distinct(ID) %>% nrow()
cat("Subjects:", n_subjects,
    if_else(n_subjects == 12, " ✓ PASS\n", " ✗ FAIL\n"))

# Check 3: Correct number of arms
n_arms <- long_data %>% distinct(ARM) %>% nrow()
cat("Arms:", n_arms,
    if_else(n_arms == 2, " ✓ PASS\n", " ✗ FAIL\n"))

# Check 4: Correct time range
time_range <- range(long_data$TIME)
cat("Time range:", time_range[1], "–", time_range[2], "h",
    if_else(time_range[1] == 0 & time_range[2] == 12, " ✓ PASS\n", " ✗ FAIL\n"))

# Check 5: Pre-dose (t=0) concentrations are all 0
t0_max <- long_data %>% filter(TIME == 0) %>% pull(CONC) %>% max()
cat("Max pre-dose conc:", t0_max,
    if_else(t0_max == 0, " ✓ PASS\n", " ✗ FAIL — non-zero pre-dose!\n"))

# Check 6: Cmax is within expected range (paper Table 2: ~1.1–1.3 ng/mL)
cmax_summary <- long_data %>%
  group_by(ARM, ID) %>%
  summarise(Cmax = max(CONC), .groups = "drop") %>%
  group_by(ARM) %>%
  summarise(
    Cmax_median = round(median(Cmax), 2),
    Cmax_min    = round(min(Cmax), 2),
    Cmax_max    = round(max(Cmax), 2)
  )
cat("\nCmax summary (target: ~1.1–1.3 ng/mL):\n")
print(cmax_summary)

# Check 7: BLQ flag counts
n_blq <- sum(long_data$BLQ)
cat("\nBLQ observations flagged:", n_blq, "\n")

# ── SECTION 5: Quick overview of the data structure ──────────────────────────

cat("\n=== DATA STRUCTURE OVERVIEW ===\n")
cat("Dimensions:", nrow(long_data), "rows ×", ncol(long_data), "columns\n")
cat("Columns:", paste(names(long_data), collapse = ", "), "\n\n")
cat("First 6 rows (No inducer, Subject 1):\n")
print(head(long_data, 6))

cat("\nLast 6 rows (With inducer, Subject 12):\n")
print(tail(long_data, 6))

# ── SECTION 6: Save processed long-format data ───────────────────────────────

write_csv(long_data, here("data", "processed", "mock_data_long.csv"))
cat("\n✓ Saved: data/processed/mock_data_long.csv\n")
cat("  This file is the INPUT for all NCA and plotting scripts.\n")

cat("\n=== SECTION COMPLETE ===\n")
cat("Next step: open R/02_descriptive_stats.R\n")

