
# Reproducible Clinical PK Analysis in R — Educational Project Based on Dunvald et al. (2022)

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org/)
[![Quarto](https://img.shields.io/badge/Quarto-%3E%3D1.3-orange)](https://quarto.org/)
[![License: CC BY-NC-ND](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)
[![Reproducible](https://img.shields.io/badge/renv-locked-green)](renv.lock)

> [!IMPORTANT]
> ## Disclaimer
>
> This repository is an independent educational and portfolio-based pharmacokinetic analysis project created for learning, reproducibility, and scientific skill development purposes.
>
> The work is inspired by and references:
>
> **Dunvald A-CD, Iversen DB, Svendsen ALO, et al. Tutorial: Statistical analysis and reporting of clinical pharmacokinetic studies. Clinical and Translational Science. 2022;15:1856–1866.**
>
> All rights to the original publication, study design, and source material remain with the respective authors and publishers.
>
> This repository does not claim ownership of the original publication or associated intellectual property. Any recreated analyses, scripts, visualizations, or workflows were independently developed for educational demonstration and reproducible research practice only.
>
> This project is not intended for clinical decision-making, regulatory submission, or medical guidance.

## Overview

This repository reproduces the complete statistical analysis from:

> Dunvald A-CD, Iversen DB, Svendsen ALO, et al.  
> **Tutorial: Statistical analysis and reporting of clinical pharmacokinetic studies.**  
> *Clinical and Translational Science.* 2022;15:1856–1866.  
> DOI: [10.1111/cts.13305](https://doi.org/10.1111/cts.13305)

The paper provides a mock dataset simulating a clinical pharmacokinetic drug–drug interaction (DDI) study: 12 healthy subjects administered 2 mg oral **midazolam** (CYP3A4 probe substrate) with and without a CYP3A inducer, in a crossover design.

## 🚧 Project Status

**Ongoing Live Learning Project**

This repository is being actively developed as part of a continuous learning journey in:

- Clinical Pharmacokinetics (PK)
- Noncompartmental Analysis (NCA)
- Reproducible Research
- R-based PK workflows
- Pharmacometric data analysis

New analyses, validation steps, visualizations, and workflow improvements may be added over time as the project evolves.

## Learning Objectives

This project was developed to practice and demonstrate:

- Clinical pharmacokinetic data analysis workflows
- Noncompartmental analysis (NCA)
- Bioequivalence and GMR calculations
- Publication-quality PK visualization in R
- Reproducible research practices using `renv` and Quarto
- Transparent scientific reporting and workflow documentation

## What this project reproduces

| Output | Paper reference | Script |
|--------|----------------|--------|
| Table 1 — Demographics | Section 4 | `R/02_descriptive_stats.R` |
| Figure 2 — Mean ± SD concentration–time curve | Section 4 | `R/03_pk_plots.R` |
| Table 2 — NCA parameters (AUClast, Cmax, t½, CL/F) | Section 4 | `R/04_nca_ncappc.R` |
| GMR + 95% CI (bioequivalence assessment) | Section 3.1 | `R/05_gmr_testing.R` |
| Figure 3 — Spaghetti plot of individual CL/F | Section 4 | `R/03_pk_plots.R` |
| Sample size calculation | Appendix S1 VII | `R/06_sample_size.R` |

**Bonus:** NCA cross-validated with `PKNCA` package alongside the paper's `ncappc` method.

## Study design

```
12 healthy subjects
        │
   Randomisation
   ┌────┴────┐
   │         │
Midazolam  Midazolam        ← Crossover (washout period between)
+ inducer  alone
   │         │
  NCA       NCA
   └────┬────┘
        │
   GMR + 95% CI
   Bioequivalence check (0.80–1.25)
```

## How to reproduce

### Prerequisites
- R ≥ 4.2 — [Download](https://www.r-project.org/)
- RStudio ≥ 2023 — [Download](https://posit.co/download/rstudio-desktop/)
- Quarto CLI ≥ 1.3 — [Download](https://quarto.org/docs/get-started/)

### Steps

```r
# 1. Clone the repository
# git clone https://github.com/pbr26/midazolam-pk-dunvald2022.git

# 2. Open the .Rproj file in RStudio (double-click it)

# 3. Restore the exact package environment
renv::restore()

# 4. Run scripts in order
source("R/00_setup.R")        # Package install (first time only)
source("R/01_load_data.R")    # Reconstruct dataset
source("R/02_descriptive_stats.R")
source("R/03_pk_plots.R")
source("R/04_nca_ncappc.R")
source("R/05_gmr_testing.R")
source("R/06_sample_size.R")

# 5. Render the full report
quarto::quarto_render("report/dunvald2022_reproduction.qmd")
```

The rendered HTML report opens automatically and is also saved to `docs/index.html`.

## Project structure

```
midazolam-pk-dunvald2022/
├── R/
│   ├── 00_setup.R              # Package installation, renv init
│   ├── 01_load_data.R          # Reconstruct Table S1, reshape to long format
│   ├── 02_descriptive_stats.R  # Table 1: demographics
│   ├── 03_pk_plots.R           # Figure 2 (mean ± SD) + Figure 3 (spaghetti)
│   ├── 04_nca_ncappc.R         # NCA with ncappc (paper method)
│   ├── 04b_nca_pknca.R         # NCA with PKNCA (cross-validation)
│   ├── 05_gmr_testing.R        # GMR, 95% CI, paired t-test, Wilcoxon
│   └── 06_sample_size.R        # Power and sample size calculation
├── data/
│   ├── raw/
│   │   └── mock_data_wide.csv  # Table S1 exactly as in paper
│   └── processed/
│       └── mock_data_long.csv  # NONMEM-style tidy format (analysis input)
├── report/
│   └── dunvald2022_reproduction.qmd   # Quarto report source
├── output/
│   ├── figures/                # Figure 2, Figure 3 (PNG + PDF)
│   └── tables/                 # Table 1, Table 2 (CSV)
├── docs/
│   └── index.html              # Rendered report (GitHub Pages)
├── renv.lock                   # Exact package versions (reproducibility)
├── session_info.txt            # R version + OS at time of analysis
└── README.md
```

## Key results (reproduced from Table 2)

| Parameter | Without inducer | With inducer | GMR (95% CI) |
|-----------|----------------|-------------|--------------|
| AUClast (ng·h/mL) | 5.99 (5.57–6.97) | 4.16 (3.88–5.08) | **0.70 (0.65–0.75)** |
| Cmax (ng/mL) | 1.28 (1.24–1.34) | 1.12 (1.10–1.20) | 0.88 (0.86–0.91) |
| t½ (h) | 2.16 (2.00–2.59) | 1.47 (1.37–1.81) | 0.68 (0.63–0.73) |
| CL/F (L/h) | 3236.67 | 4780.99 | 1.48 (1.38–1.59) |

GMR for AUClast = **0.70** — outside the 0.80–1.25 bioequivalence boundary, confirming a statistically and clinically significant CYP3A induction effect.

## Live report

[View the rendered report →](https://pbr26.github.io/midazolam-pk-dunvald2022/)

## Skills demonstrated

- Noncompartmental analysis (NCA) — `ncappc` and `PKNCA`
- Drug–drug interaction (DDI) assessment — GMR with 95% CI
- Bioequivalence analysis — no-effect boundary interpretation
- Regulatory-style reporting (Table 1, Table 2, spaghetti plot)
- Reproducible research — `renv`, Quarto, GitHub Actions CI/CD

## License

The original paper is open access under CC BY-NC-ND 4.0.  
This reproduction code is shared for educational purposes.  
Please cite the original paper (Dunvald et al. 2022) if you build on this work.

## Citation

```bibtex
@article{dunvald2022,
  author  = {Dunvald, Ann-Cathrine Dalg{\aa}rd and others},
  title   = {Tutorial: Statistical analysis and reporting of clinical
             pharmacokinetic studies},
  journal = {Clinical and Translational Science},
  year    = {2022},
  volume  = {15},
  pages   = {1856--1866},
  doi     = {10.1111/cts.13305}
}gi
```
=======
# midazolam-pk-dunvald2022
Reproduction of Dunvald 2022 clinical PK tutorial in R

## Acknowledgment

This project is an independent educational effort and is not affiliated with, endorsed by, or sponsored by the original authors, journal, or publishers.
>>>>>>> 8c4fd1451b11b9497e456fa58aac0ce51bb56f54


