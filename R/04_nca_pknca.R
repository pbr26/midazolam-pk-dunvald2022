# =============================================================================
# 04_nca_pknca.R
# Project : Reproduction of Dunvald et al. (2022)
# Purpose : Noncompartmental Analysis using PKNCA
#           Reproduces Table 2 — AUClast, Cmax, t½, CL/F, Tmax
#           Computes GMR with 95% CI (primary DDI result)
# Author  : Pramod BR
# Date    : 2026
#
# WHY PKNCA INSTEAD OF ncappc?
#   ncappc v1.0.0 changed its API completely — old argument names no longer
#   work. PKNCA is the modern standard used in pharma regulatory submissions.
#   It produces identical NCA results and its output is SDTM-ready.
#
# WHAT THIS SCRIPT PRODUCES:
#   output/tables/nca_parameters_individual.csv  — 24 rows (12 subj x 2 arms)
#   output/tables/gmr_results.csv                — GMR + 95% CI
#   output/tables/table2_nca_summary.csv         — final Table 2
#
# TARGET VALUES FROM TABLE 2 (Dunvald 2022):
#   AUClast GMR : 0.70 (0.65-0.75)
#   Cmax    GMR : 0.88 (0.86-0.91)
#   t½      GMR : 0.68 (0.63-0.73)
#   CL/F    GMR : 1.48 (1.38-1.59)
# =============================================================================


# ── STEP 1: Load libraries ────────────────────────────────────────────────────
# Run each line one at a time with Cmd+Enter
# No output = success for each library() call

library(PKNCA)   # NCA calculations — modern, SDTM-ready
library(dplyr)   # Data manipulation
library(tidyr)   # pivot_wider, separate
library(readr)   # read_csv, write_csv
library(here)    # Portable file paths

cat("✓ All libraries loaded\n")

# MY UNDERSTANDING: ___________________________________________________________


# ── STEP 2: Load the dataset ──────────────────────────────────────────────────
# This file was created by 01_load_data.R
# If you get a file not found error here — run source("R/01_load_data.R") first

long_data <- read_csv(
  here("data", "processed", "mock_data_long.csv"),
  show_col_types = FALSE
)

cat("✓ Dataset loaded:", nrow(long_data), "rows\n")
cat("  Expected: 288 rows (12 subjects x 2 arms x 12 timepoints)\n")


# ── STEP 3: Prepare concentration data for PKNCA ─────────────────────────────
#
# WHAT THIS DOES:
#   Creates a unique subject ID per arm called SID
#   Subject 1 No_inducer  → "1_No_inducer"
#   Subject 1 With_inducer → "1_With_inducer"
#   This gives PKNCA 24 separate PK profiles to analyse
#
#   Also removes BLQ observations (except at t=0 which is correctly 0)
#   BLQ = Below Limit of Quantification — not a real measured value

conc_df <- long_data %>%
  filter(BLQ == 0 | TIME == 0) %>%
  mutate(SID = paste0(ID, "_", ARM)) %>%
  select(SID, TIME, CONC, DOSE) %>%
  arrange(SID, TIME)

cat("\n✓ Concentration data prepared\n")
cat("  Rows:", nrow(conc_df), "\n")
cat("  Unique profiles:", n_distinct(conc_df$SID),
    "(expected 24 = 12 subjects x 2 arms)\n")
cat("\n  First 4 rows:\n")
print(head(conc_df, 4))


# ── STEP 4: Prepare dose data ─────────────────────────────────────────────────
# PKNCA needs a separate dose dataframe
# Dose = 2 mg oral midazolam at TIME = 0 for every profile

dose_df <- conc_df %>%
  distinct(SID) %>%
  mutate(TIME = 0, DOSE = 2)

cat("\n✓ Dose data prepared\n")
cat("  Rows:", nrow(dose_df), "(one per profile = 24)\n")


# ── STEP 5: Build PKNCA objects ───────────────────────────────────────────────
#
# PKNCA works in 3 objects:
#   PKNCAconc  — tells PKNCA: "CONC is concentration, TIME is time, SID groups"
#   PKNCAdose  — tells PKNCA: "DOSE is dose, given at TIME=0, per SID"
#   PKNCAdata  — combines both + tells PKNCA what parameters to calculate
#
# The formula CONC ~ TIME | SID means:
#   CONC as a function of TIME, grouped by SID

conc_obj <- PKNCAconc(conc_df, CONC ~ TIME | SID)
dose_obj <- PKNCAdose(dose_df, DOSE ~ TIME | SID)

cat("\n✓ PKNCA concentration object built\n")
cat("✓ PKNCA dose object built\n")


# ── STEP 6: Set calculation intervals ─────────────────────────────────────────
#
# This tells PKNCA:
#   - Calculate from time 0 to time 12 (our full sampling window)
#   - Calculate these specific parameters: TRUE means calculate it
#
# auclast   = AUC from t=0 to last quantifiable concentration
# cmax      = Maximum observed concentration
# tmax      = Time of maximum concentration
# half.life = Terminal elimination half-life (t½)
# cl.obs    = Apparent clearance CL/F (for oral/extravascular dosing)

my_intervals <- data.frame(
  start     = 0,
  end       = 12,
  auclast   = TRUE,
  cmax      = TRUE,
  tmax      = TRUE,
  half.life = TRUE,
  cl.obs    = TRUE
)

# Build the final PKNCAdata object with intervals
pk_data <- PKNCAdata(
  data.conc = conc_obj,
  data.dose = dose_obj,
  intervals = my_intervals
)

cat("✓ PKNCA data object built with intervals\n")


# ── STEP 7: Run the NCA ───────────────────────────────────────────────────────
#
# This is the main calculation step.
# pk.nca() runs NCA for all 24 profiles automatically.
# Takes 20-30 seconds — you will see some messages, that is normal.
#
# WHAT HAPPENS INSIDE:
#   For each profile:
#   1. Reads Cmax = highest observed CONC value
#   2. Reads Tmax = TIME at which Cmax occurred
#   3. Calculates AUClast using linear-up/log-down trapezoidal method
#   4. Fits regression line to terminal log-linear phase → gets t½
#   5. Calculates CL/F = Dose / AUCinf

cat("\n── Running NCA (20-30 seconds) ──────────────────────────\n")
nca_result <- pk.nca(pk_data)
cat("✓ NCA complete\n")
cat("────────────────────────────────────────────────────────\n")


# ── STEP 8: Extract results ───────────────────────────────────────────────────
#
# PKNCA stores results in a long format — one row per parameter per profile
# We filter to the 5 parameters we need, then pivot to wide format
# (one row per profile, one column per parameter)

nca_raw <- as.data.frame(nca_result$result)

cat("\nParameters calculated by PKNCA:\n")
cat(paste(unique(nca_raw$PPTESTCD), collapse=", "), "\n")

nca_clean <- nca_raw %>%
  # Keep only the 5 parameters that match Table 2
  filter(PPTESTCD %in% c("auclast","cmax","tmax","half.life","cl.obs")) %>%
  select(SID, PPTESTCD, PPORRES) %>%

  # Reshape: one row per SID, one column per parameter
  pivot_wider(
    names_from  = PPTESTCD,
    values_from = PPORRES
  ) %>%

  # Split SID "1_No_inducer" back into ID=1, ARM=No_inducer
  separate(SID,
           into  = c("ID", "ARM"),
           sep   = "_",
           extra = "merge") %>%

  mutate(
    ID = as.integer(ID),

    # Convert CL/F units to L/h
    # PKNCA returns in mg/(ng/mL) = 1000 L — multiply by 1000
    cl.obs = cl.obs * 1000
  ) %>%

  # Rename columns to match our project convention
  rename(
    AUClast = auclast,
    Cmax    = cmax,
    Tmax    = tmax,
    t_half  = half.life,
    CLF     = cl.obs
  ) %>%
  arrange(ARM, ID)

cat("\n✓ Results extracted\n")
cat("  Rows:", nrow(nca_clean), "(expected 24)\n")
cat("\n  First 6 rows:\n")
print(head(nca_clean, 6))


# ── STEP 9: Compute Table 2 summary statistics ────────────────────────────────
#
# Table 2 reports: median (IQR) per arm for each parameter
# IQR = interquartile range = 25th to 75th percentile
#
# WHY MEDIAN not MEAN?
#   PK parameters with small samples (n=12) can be skewed by one outlier
#   Median is resistant to outliers — it is the preferred measure

fmt <- function(x) {
  x <- na.omit(x)
  paste0(
    round(median(x), 2),
    " (", round(quantile(x, 0.25), 2),
    "-",  round(quantile(x, 0.75), 2), ")"
  )
}

cat("\n=== TABLE 2 — median (IQR) per arm ===\n")
summary_table <- nca_clean %>%
  group_by(ARM) %>%
  summarise(
    AUClast_nghmL = fmt(AUClast),
    Cmax_ngmL     = fmt(Cmax),
    t_half_h      = fmt(t_half),
    CLF_Lh        = fmt(CLF),
    .groups       = "drop"
  )
print(summary_table)

cat("\n── COMPARE TO PAPER (Table 2, Dunvald 2022) ──────────────\n")
cat("Paper No inducer  AUClast: 5.99 (5.57-6.97)\n")
cat("Paper With inducer AUClast: 4.16 (3.88-5.08)\n")
cat("──────────────────────────────────────────────────────────\n")


# ── STEP 10: Compute GMR with 95% CI ─────────────────────────────────────────
#
# GMR = Geometric Mean Ratio = with inducer / without inducer
#
# HOW IT IS CALCULATED:
#   1. Log-transform both arms' values
#   2. Run a paired t-test (same subjects in both arms = paired)
#   3. The mean difference on log scale = log(GMR)
#   4. Back-transform: exp(log(GMR)) = GMR
#   5. exp(95% CI on log scale) = 95% CI on GMR scale
#
# INTERPRETATION:
#   GMR = 0.70 → 30% reduction in exposure (induction)
#   GMR = 1.00 → no change
#   GMR = 1.50 → 50% increase (inhibition)
#
# BIOEQUIVALENCE BOUNDARY: 0.80 to 1.25
#   If ENTIRE 95% CI is within 0.80-1.25 → no clinically relevant DDI
#   If CI crosses outside → clinically significant DDI

no_ind   <- nca_clean %>% filter(ARM == "No_inducer")   %>% arrange(ID)
with_ind <- nca_clean %>% filter(ARM == "With_inducer") %>% arrange(ID)

gmr_calc <- function(param, x1, x2) {
  x1 <- na.omit(x1)
  x2 <- na.omit(x2)
  tt <- t.test(log(x2), log(x1), paired = TRUE)
  tibble(
    Parameter  = param,
    GMR        = round(exp(tt$estimate), 2),
    CI_lower   = round(exp(tt$conf.int[1]), 2),
    CI_upper   = round(exp(tt$conf.int[2]), 2),
    GMR_CI     = paste0(
                   round(exp(tt$estimate), 2),
                   " (", round(exp(tt$conf.int[1]), 2),
                   "-",  round(exp(tt$conf.int[2]), 2), ")"
                 ),
    p_value    = round(tt$p.value, 4),
    Outside_BE = ifelse(
                   exp(tt$conf.int[1]) < 0.80 |
                   exp(tt$conf.int[2]) > 1.25,
                   "Yes - significant DDI",
                   "No - within boundary"
                 )
  )
}

gmr_results <- bind_rows(
  gmr_calc("AUClast", no_ind$AUClast, with_ind$AUClast),
  gmr_calc("Cmax",    no_ind$Cmax,    with_ind$Cmax),
  gmr_calc("t_half",  no_ind$t_half,  with_ind$t_half),
  gmr_calc("CL_F",    no_ind$CLF,     with_ind$CLF)
)

cat("\n=== GMR RESULTS ===\n")
print(gmr_results %>% select(Parameter, GMR_CI, p_value, Outside_BE))

auc <- gmr_results %>% filter(Parameter == "AUClast")

cat("\n── PRIMARY RESULT ──────────────────────────────────────────\n")
cat("AUClast GMR:", auc$GMR,
    "(95% CI:", auc$CI_lower, "-", auc$CI_upper, ")\n")
cat("Exposure reduction:", round((1 - auc$GMR) * 100, 0), "%\n")
cat("Significant DDI?  ", auc$Outside_BE, "\n")
cat("────────────────────────────────────────────────────────────\n")

cat("\n── COMPARE TO PAPER ────────────────────────────────────────\n")
cat("Paper AUClast GMR : 0.70 (0.65-0.75)\n")
cat("Your  AUClast GMR :", auc$GMR_CI, "\n")
cat("Paper Cmax    GMR : 0.88 (0.86-0.91)\n")
cat("Your  Cmax    GMR :",
    gmr_results %>% filter(Parameter=="Cmax") %>% pull(GMR_CI), "\n")
cat("Paper t½      GMR : 0.68 (0.63-0.73)\n")
cat("Your  t½      GMR :",
    gmr_results %>% filter(Parameter=="t_half") %>% pull(GMR_CI), "\n")
cat("────────────────────────────────────────────────────────────\n")


# ── STEP 11: Tmax — Wilcoxon test ─────────────────────────────────────────────
#
# WHY WILCOXON for Tmax?
#   Tmax is categorical — it can only be one of the sampling timepoints
#   0.5, 1, 1.5, 2, 3, 4... It is not normally distributed.
#   Wilcoxon rank-sum test does not assume normality.
#   The paper (p.1860) specifically uses Wilcoxon for Tmax.

wt <- wilcox.test(
  no_ind$Tmax,
  with_ind$Tmax,
  paired = TRUE,
  exact  = FALSE
)

cat("\n=== TMAX ANALYSIS ===\n")
cat("No inducer   Tmax: median",
    median(no_ind$Tmax, na.rm=TRUE), "h",
    "range", min(no_ind$Tmax, na.rm=TRUE),
    "-", max(no_ind$Tmax, na.rm=TRUE), "h\n")
cat("With inducer Tmax: median",
    median(with_ind$Tmax, na.rm=TRUE), "h",
    "range", min(with_ind$Tmax, na.rm=TRUE),
    "-", max(with_ind$Tmax, na.rm=TRUE), "h\n")
cat("Wilcoxon p-value:", round(wt$p.value, 4), "\n")
cat("Significant (p<0.05)?",
    ifelse(wt$p.value < 0.05, "Yes", "No"), "\n")


# ── STEP 12: Build final Table 2 ──────────────────────────────────────────────

table2_final <- tibble(
  Parameter = c(
    "AUClast (ng.h/mL)",
    "Cmax (ng/mL)",
    "t1/2 (h)",
    "CL/F (L/h)",
    "Tmax (h)"
  ),
  `Without inducer (median IQR)` = c(
    fmt(no_ind$AUClast),
    fmt(no_ind$Cmax),
    fmt(no_ind$t_half),
    fmt(no_ind$CLF),
    paste0(median(no_ind$Tmax, na.rm=TRUE),
           " (", min(no_ind$Tmax, na.rm=TRUE),
           "-",  max(no_ind$Tmax, na.rm=TRUE), ")")
  ),
  `With inducer (median IQR)` = c(
    fmt(with_ind$AUClast),
    fmt(with_ind$Cmax),
    fmt(with_ind$t_half),
    fmt(with_ind$CLF),
    paste0(median(with_ind$Tmax, na.rm=TRUE),
           " (", min(with_ind$Tmax, na.rm=TRUE),
           "-",  max(with_ind$Tmax, na.rm=TRUE), ")*")
  ),
  `GMR (95% CI)` = c(
    gmr_results %>% filter(Parameter=="AUClast") %>% pull(GMR_CI),
    gmr_results %>% filter(Parameter=="Cmax")    %>% pull(GMR_CI),
    gmr_results %>% filter(Parameter=="t_half")  %>% pull(GMR_CI),
    gmr_results %>% filter(Parameter=="CL_F")    %>% pull(GMR_CI),
    "NA"
  )
)

cat("\n=== FINAL TABLE 2 ===\n")
print(table2_final)


# ── STEP 13: Save all outputs ─────────────────────────────────────────────────

write_csv(nca_clean,
          here("output", "tables", "nca_parameters_individual.csv"))
write_csv(gmr_results,
          here("output", "tables", "gmr_results.csv"))
write_csv(table2_final,
          here("output", "tables", "table2_nca_summary.csv"))

cat("\n✓ Saved: output/tables/nca_parameters_individual.csv\n")
cat("✓ Saved: output/tables/gmr_results.csv\n")
cat("✓ Saved: output/tables/table2_nca_summary.csv\n")
cat("\n=== PHASE 3 COMPLETE ===\n")
cat("Check output/tables/ for your three new files\n")
cat("Then run: quarto render report/dunvald2022_reproduction.qmd\n")
