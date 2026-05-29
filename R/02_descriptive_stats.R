# =============================================================================
# 02_descriptive_stats.R
# Project : Reproduction of Dunvald et al. (2022)
# Purpose : Compute Table 1 — descriptive / demographic statistics
# Author  : Pramod BR
# Date    : 2026-05-26
#
# WHAT THIS SCRIPT DOES (plain English):
#   Before reporting any PK results, we describe the study population —
#   their age, weight, sex — so a reader knows who the results apply to.
#   This is always Table 1 in a clinical paper.
#
#   In Dunvald 2022 subjects are SIMULATED, so we create realistic
#   demographics matching a typical healthy-volunteer CYP3A DDI study.
#   In a real study you would load this from a demographics CSV file.
#
# PAPER REFERENCE: Section "Considerations for descriptive statistics"
#   and Table 1 template (p.1861, Dunvald 2022)
# =============================================================================


# ── STEP 0: Load libraries ────────────────────────────────────────────────────
# library() activates an installed package for this R session.
# Think of install.packages() as downloading an app,
# and library() as opening it. You install once, but open every session.

library(dplyr)      # Data manipulation: filter(), mutate(), group_by(), summarise()
library(tidyr)      # Data reshaping (used in later scripts too)
library(readr)      # Fast CSV reading/writing: read_csv(), write_csv()
library(tableone)   # Generates clinical Table 1 automatically
library(knitr)      # kable() — renders tables in the Quarto report
library(kableExtra) # Adds styling to kable tables (stripes, bold headers)
library(here)       # Portable file paths — here() always resolves from project root

# MY UNDERSTANDING: (write your own explanation here after reading)
# ___________________________________________________________________________


# ── STEP 1: Load the dataset created in 01_load_data.R ───────────────────────

# here("data", "processed", "mock_data_long.csv") builds the path:
#   Windows: C:\Projects\midazolam-pk-dunvald2022\data\processed\mock_data_long.csv
#   Mac/Linux: ~/Projects/midazolam-pk-dunvald2022/data/processed/mock_data_long.csv
# It works on ANY machine without editing — that is the whole point of here()

long_data <- read_csv(
  here("data", "processed", "mock_data_long.csv"),
  show_col_types = FALSE   # Suppresses the column-type printout in Console
)

cat("Dataset loaded:", nrow(long_data), "rows x", ncol(long_data), "cols\n")
# Expected: 288 rows (12 subjects x 2 arms x 12 timepoints), 8 columns


# ── STEP 2: Simulate demographics ────────────────────────────────────────────
#
# CLINICAL PK CONCEPT — Why do demographics matter in PK?
# ─────────────────────────────────────────────────────────
# Demographics are NOT just background — they directly affect PK parameters:
#
# BODY WEIGHT → affects Volume of Distribution (Vd = how widely drug spreads)
#   Heavier person → drug distributes into more tissue → lower Cmax
#
# AGE → affects hepatic blood flow and enzyme expression
#   Older patients often have lower CYP3A4 activity → higher AUC, longer t½
#
# SEX → women have ~20-40% higher CYP3A4 activity than men on average
#   This causes meaningfully different midazolam clearance between sexes
#
# By reporting demographics, we tell the reader:
#   "Our results apply to people like THESE 12 subjects."
# If all 12 were elderly women, results may not apply to young men.
#
# Dunvald 2022 Box 2 explicitly urges including both sexes in PK studies.
# ─────────────────────────────────────────────────────────────────────────

set.seed(42)
# set.seed() makes random number generation REPRODUCIBLE.
# Without it: every time you run this script, you get slightly different ages/weights.
# With set.seed(42): the script produces IDENTICAL numbers on every machine, every time.
# This is essential for reproducible research.
# The number 42 is arbitrary — any integer works.

n_subjects <- 12

demo <- tibble(
  # tibble() creates a modern data frame
  # It prints more cleanly than base data.frame() and handles strings better

  ID = 1:n_subjects,
  # 1:12 creates the sequence 1, 2, 3, ..., 12

  # AGE: rnorm(n, mean, sd) draws n values from a normal distribution
  # round(..., 0) rounds to whole years
  # Healthy volunteer DDI studies typically enrol 18-45 year olds
  Age_years = round(rnorm(n_subjects, mean = 28, sd = 5), 0),

  # SEX: balanced — 6 male, 6 female
  # rep("Male", 6) creates: "Male" "Male" "Male" "Male" "Male" "Male"
  # c() combines two rep() vectors into one 12-element vector
  Sex = c(rep("Male", 6), rep("Female", 6)),

  # WEIGHT: different distributions for male vs female
  # ifelse(condition, value_if_TRUE, value_if_FALSE) works element-by-element
  Weight_kg = round(
    ifelse(
      Sex == "Male",
      rnorm(n_subjects, mean = 75, sd = 8),   # Men: mean 75 kg
      rnorm(n_subjects, mean = 63, sd = 7)    # Women: mean 63 kg
    ), 1),

  # HEIGHT: needed to compute BMI
  Height_cm = round(
    ifelse(Sex == "Male",
           rnorm(n_subjects, mean = 178, sd = 6),
           rnorm(n_subjects, mean = 165, sd = 5)), 0),

  # eGFR: kidney function (mL/min/1.73m²)
  # Normal is >90. Healthy volunteers always have normal eGFR.
  # CYP3A4 is a liver enzyme — why report eGFR?
  # Because impaired kidneys affect drug protein binding and active secretion,
  # which can alter free drug concentration even for hepatically-cleared drugs.
  eGFR_ml_min = round(rnorm(n_subjects, mean = 105, sd = 12), 0)

) %>%
  # mutate() adds new columns or transforms existing ones
  # The pipe %>% means: "take the result above and pass it into the next function"
  mutate(
    # BMI computed AFTER height and weight columns exist
    # Formula: weight(kg) / height(m)^2
    # Height_cm / 100 converts centimetres to metres
    BMI = round(Weight_kg / (Height_cm / 100)^2, 1),

    # Clamp age to realistic healthy volunteer range
    # pmax(18, x) replaces any value below 18 with 18 (lower bound)
    # pmin(45, x) replaces any value above 45 with 45 (upper bound)
    # pmax/pmin are "parallel" max/min — they work element-by-element
    Age_years   = pmax(18, pmin(45, Age_years)),
    eGFR_ml_min = pmax(80, pmin(140, eGFR_ml_min))
  )

cat("\n=== Demographics preview ===\n")
print(demo)


# ── STEP 3: Compute summary statistics manually ───────────────────────────────
#
# CLINICAL PK CONCEPT — median vs mean, and WHY it matters for small samples
# ───────────────────────────────────────────────────────────────────────────
# With 12 subjects, ONE outlier (e.g. one very heavy person) pulls the mean
# away from where most subjects actually are. The median is resistant to this.
#
# IQR = interquartile range = distance from 25th to 75th percentile
# It tells you where the MIDDLE 50% of your data falls.
#
# Example: Age median 27, IQR 24-32 means:
#   - Half the subjects were between 24 and 32 years old
#   - One subject aged 45 does not distort this at all
#
# Dunvald 2022 p.1861: "median is preferred in situations with small
# sample sizes, non-normal distribution, or if there are outliers."
#
# Reporting format: "27 (24-32)" — always median (Q1-Q3)
# NEVER use SEM (standard error of the mean) in descriptive stats —
# the paper explicitly warns against this (p.1861)
# ───────────────────────────────────────────────────────────────────────────

# Helper function 1: formats continuous variable as "median (Q1-Q3)"
# We define this ONCE and reuse it — avoids copy-paste errors
fmt_median_iqr <- function(x) {
  # function(x) defines a custom function that takes one input: x
  q1  <- round(quantile(x, 0.25, na.rm = TRUE), 1)  # 25th percentile
  med <- round(median(x,           na.rm = TRUE), 1)  # 50th percentile
  q3  <- round(quantile(x, 0.75, na.rm = TRUE), 1)  # 75th percentile
  paste0(med, " (", q1, "-", q3, ")")
  # paste0() joins strings with NO separator between them
  # Result: "27 (24.2-31.8)"
}

# Helper function 2: formats categorical variable as "n (%)"
fmt_n_pct <- function(x, level) {
  # x = the full vector (e.g. demo$Sex)
  # level = which category to count (e.g. "Male")
  n_level <- sum(x == level, na.rm = TRUE)   # How many match this level
  pct     <- round(n_level / length(x) * 100, 0)
  paste0(n_level, " (", pct, "%)")
  # Result: "6 (50%)"
}

# Build Table 1 — one row per characteristic
# This is a crossover self-controlled study, so there is only ONE column
# (the same 12 subjects appear in both arms — their demographics don't change)
table1 <- tibble(
  Characteristic = c(
    "Number of subjects, n",
    "Age, years",
    "Sex — Male, n (%)",
    "Sex — Female, n (%)",
    "Body weight, kg",
    "BMI, kg/m²",
    "Height, cm",
    "eGFR, mL/min/1.73m²"
  ),
  Value = c(
    as.character(n_subjects),             # "12"
    fmt_median_iqr(demo$Age_years),       # "28 (24-32)"
    fmt_n_pct(demo$Sex, "Male"),          # "6 (50%)"
    fmt_n_pct(demo$Sex, "Female"),        # "6 (50%)"
    fmt_median_iqr(demo$Weight_kg),       # e.g. "72 (65-78)"
    fmt_median_iqr(demo$BMI),             # e.g. "23 (21-25)"
    fmt_median_iqr(demo$Height_cm),       # e.g. "174 (168-179)"
    fmt_median_iqr(demo$eGFR_ml_min)      # e.g. "105 (96-115)"
  ),
  Note = c(
    "",
    "Median (IQR)",
    "", "",
    "Median (IQR)",
    "Median (IQR)",
    "Median (IQR)",
    "Median (IQR)"
  )
)

cat("\n=== TABLE 1: Demographic and clinical characteristics ===\n")
print(table1)


# ── STEP 4: Cross-check using the tableone package ────────────────────────────
# tableone automates Table 1 generation — it is used in virtually all
# clinical trials and epidemiology papers.
# We use it here as a CROSS-CHECK: if our manual table matches tableone,
# we know our helper functions are correct.

cat_vars <- "Sex"  # Variables to treat as categorical (get n% not mean±SD)

table1_pkg <- CreateTableOne(
  vars       = c("Age_years", "Sex", "Weight_kg", "BMI", "eGFR_ml_min"),
  data       = demo,
  factorVars = cat_vars
)

cat("\n=== tableone cross-check ===\n")
# nonnormal specifies which variables to report as median (IQR)
print(table1_pkg,
      nonnormal      = c("Age_years", "Weight_kg", "BMI", "eGFR_ml_min"),
      showAllLevels  = TRUE,
      quote          = FALSE,
      noSpaces       = TRUE)


# ── STEP 5: Compute mean ± SD concentration per timepoint per arm ─────────────
#
# CLINICAL PK CONCEPT — Why mean concentrations?
# ────────────────────────────────────────────────
# Figure 2 in the paper plots mean ± SD at each timepoint.
# This answers: "What does a TYPICAL patient's concentration look like?"
# The SD bars answer: "How much do patients vary around that typical value?"
#
# We compute this HERE (not in the plotting script) because:
# 1. Separation of concerns — data prep vs visualisation
# 2. The summary data can be reused in the Quarto report directly
# 3. Easier to check: print and inspect before plotting
# ────────────────────────────────────────────────

mean_conc <- long_data %>%

  filter(BLQ == 0) %>%
  # BLQ == 1 means concentration was Below the Limit of Quantification
  # These are not real measured values — the assay couldn't detect the drug
  # Including them in the mean would pull concentrations DOWN artificially
  # So we exclude them before computing means
  # (In a real study, BLQ handling follows FDA/EMA guidance — often set to
  # LLOQ/2 for NCA, but excluded from descriptive means)

  group_by(ARM, TIME) %>%
  # group_by() splits the data into groups for the next summarise()
  # Here: one group per ARM-TIME combination
  # With 2 arms and 12 timepoints, this creates 24 groups

  summarise(
    mean_conc = mean(CONC),         # Arithmetic mean concentration (ng/mL)
    sd_conc   = sd(CONC),           # Standard deviation
    n         = n(),                # Number of observations in this group
    se_conc   = sd_conc / sqrt(n),  # Standard error (for reference)
    .groups   = "drop"
    # .groups = "drop" removes the grouping structure after summarise()
    # Without this, the result stays grouped and behaves unexpectedly in plots
  )

cat("\n=== Mean concentrations by ARM and TIME (first 8 rows) ===\n")
print(head(mean_conc, 8))
cat("  Total rows:", nrow(mean_conc), "(expected: 24 = 2 arms x 12 timepoints)\n")


# ── STEP 6: Save all outputs ──────────────────────────────────────────────────

write_csv(demo,       here("data", "processed", "demographics.csv"))
write_csv(table1,     here("output", "tables",  "table1_demographics.csv"))
write_csv(mean_conc,  here("data", "processed", "mean_concentrations.csv"))

cat("\n✓ Saved: data/processed/demographics.csv\n")
cat("✓ Saved: output/tables/table1_demographics.csv\n")
cat("✓ Saved: data/processed/mean_concentrations.csv\n")
cat("\n=== 02_descriptive_stats.R COMPLETE ===\n")
cat("Next: open R/03_pk_plots.R\n")
