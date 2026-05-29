# =============================================================================
# 03_pk_plots.R
# Project : Reproduction of Dunvald et al. (2022)
# Purpose : Reproduce Figure 2 (mean ± SD concentration-time curve)
#           and Figure 3 (spaghetti plot of individual CL/F)
# Author  : [Your Name]
# Date    : 2026-05-26
#
# WHAT THIS SCRIPT PRODUCES:
#   Figure 2 — Mean midazolam plasma concentration vs time, with SD error bars
#               and a log-scale inset (exactly as in the paper)
#   Figure 3 — Individual CL/F per subject, connected across arms (spaghetti)
#
# PAPER REFERENCE: Section "Considerations for tables and figures" (p.1863-1865)
#   Figure 2 caption: "Midazolam plasma concentrations are presented as mean±SD"
#   Figure 3 caption: "Spaghetti plot illustrating changes in individual
#                      midazolam clearances with and without inducer"
# =============================================================================


# ── STEP 0: Load libraries ────────────────────────────────────────────────────

library(ggplot2)    # The main plotting package — grammar of graphics
library(dplyr)      # Data manipulation
library(readr)      # CSV reading
library(scales)     # log10 axis labels (label_log())
library(patchwork)  # Combines multiple ggplots into one figure (main + inset)
library(here)       # Portable file paths

# WHAT IS ggplot2?
# ggplot2 is built on the "grammar of graphics" — the idea that every plot
# is made of layers:
#   1. Data (what you are plotting)
#   2. Aesthetics / aes() — which columns map to x, y, colour, shape
#   3. Geoms — what shape to draw (line, point, errorbar, etc.)
#   4. Scales — how axes and colours are formatted
#   5. Theme — fonts, grid lines, background
# You build a plot by adding layers with the + operator.

# MY UNDERSTANDING: (write your own explanation here)
# ___________________________________________________________________________


# ── STEP 1: Load data ─────────────────────────────────────────────────────────

# Mean concentrations computed in 02_descriptive_stats.R
mean_conc <- read_csv(here("data", "processed", "mean_concentrations.csv"),
                      show_col_types = FALSE)

# Full long-format data (individual subject concentrations)
long_data <- read_csv(here("data", "processed", "mock_data_long.csv"),
                      show_col_types = FALSE)

cat("mean_conc rows:", nrow(mean_conc), "(expected 24)\n")
cat("long_data rows:", nrow(long_data), "(expected 288)\n")


# ── STEP 2: Define colour palette ─────────────────────────────────────────────
# The paper uses orange (no inducer) and green (with inducer) — Fig 2 caption
# We match this exactly so our reproduction is visually faithful

col_no_ind   <- "#D95F02"   # Orange — matches paper's "No inducer" line
col_with_ind <- "#2A9E4F"   # Green  — matches paper's "Inducer" line

# Named vector: names must EXACTLY match the ARM values in your dataset
# This lets ggplot2 assign colours automatically by ARM name
arm_colours <- c("No_inducer" = col_no_ind, "With_inducer" = col_with_ind)

# Human-readable labels for the legend
arm_labels  <- c("No_inducer" = "No inducer", "With_inducer" = "Inducer")


# ── STEP 3: Build Figure 2 — main concentration-time plot ─────────────────────
#
# CLINICAL PK CONCEPT — What does this plot tell us?
# ───────────────────────────────────────────────────
# The concentration-time curve is the FOUNDATION of PK analysis.
# It visually shows:
#
#   Cmax  = the peak of the curve (highest concentration reached)
#           → tells us how high the drug gets in the blood
#
#   Tmax  = the time at which Cmax occurs (where the peak is on the x-axis)
#           → tells us how fast the drug is absorbed
#
#   AUC   = the area under the entire curve
#           → tells us total drug exposure (how much drug overall)
#
#   t½    = the time for concentration to halve on the descending part
#           → tells us how fast the drug is eliminated
#
# When you compare the orange vs green lines:
#   - Green line (with inducer) has LOWER concentrations throughout
#   - This means: inducer increased CYP3A4 enzyme activity → faster metabolism
#     → drug eliminated faster → lower AUC (less exposure)
#   - That is the DDI effect we are measuring
#
# The LOG-SCALE INSET (smaller graph inside the main graph) is important because:
#   - Linear scale: the descending tail looks flat (can't see the slope)
#   - Log scale: the descending tail becomes a STRAIGHT LINE
#   - The slope of that line = elimination rate constant (ke)
#   - If green line has a STEEPER slope → faster elimination with inducer
#   - This visually confirms the induction effect on ke and t½
# ───────────────────────────────────────────────────────────────────────────

fig2_main <- ggplot(
  data    = mean_conc,
  mapping = aes(
    x      = TIME,         # x-axis: time after dose (hours)
    y      = mean_conc,    # y-axis: mean plasma concentration (ng/mL)
    colour = ARM,          # Line colour differs by treatment arm
    shape  = ARM,          # Point shape also differs (accessibility)
    group  = ARM           # Connect points within the same ARM
  )
) +

  # ── Geom 1: Error bars (SD) ──
  # geom_errorbar draws vertical lines at each timepoint showing ± SD
  # aes(ymin, ymax) defines the bottom and top of each error bar
  geom_errorbar(
    aes(
      ymin = mean_conc - sd_conc,   # bottom of bar: mean minus 1 SD
      ymax = mean_conc + sd_conc    # top of bar:    mean plus 1 SD
    ),
    width     = 0.2,    # horizontal width of the cap at top/bottom of bar
    linewidth = 0.4,    # thickness of the error bar line
    alpha     = 0.7     # slight transparency so bars don't overpower the line
  ) +

  # ── Geom 2: Line connecting timepoints ──
  # geom_line draws a connected line through all points within the same ARM
  geom_line(linewidth = 0.8) +

  # ── Geom 3: Points at each measured timepoint ──
  # geom_point draws a symbol at each observed timepoint
  geom_point(size = 2.5, fill = "white", stroke = 0.8) +

  # ── Scale: Custom colours ──
  # scale_colour_manual() overrides the default colours with our named vector
  # scale_shape_manual() sets different point shapes per arm
  scale_colour_manual(
    values = arm_colours,
    labels = arm_labels,
    name   = NULL         # No legend title — matches paper style
  ) +
  scale_shape_manual(
    values = c("No_inducer" = 24, "With_inducer" = 21),
    # 24 = filled triangle up, 21 = filled circle
    # These match the paper: "triangled points" and "circled points"
    labels = arm_labels,
    name   = NULL
  ) +

  # ── Axes ──
  # scale_x_continuous() customises the x-axis
  scale_x_continuous(
    breaks = c(0, 2, 4, 6, 8, 10, 12),  # Tick marks at these hours
    limits = c(0, 12)
  ) +
  scale_y_continuous(
    limits = c(0, 1.6),      # y-axis from 0 to 1.6 ng/mL
    breaks = seq(0, 1.5, by = 0.5)
  ) +

  # ── Labels ──
  labs(
    title   = "Mean midazolam plasma concentration",
    x       = "Time after dose (h)",
    y       = "Concentration (ng/mL)",
    caption = "Data presented as mean ± SD. n = 12 healthy subjects."
  ) +

  # ── Legend position ──
  # Place legend inside the plot (bottom-right) to maximise plot area
  # This matches the paper's Figure 2 layout
  theme_bw() +
  # theme_bw() = black-and-white theme — clean, publication-ready
  # It removes the grey background of the default ggplot2 theme

  theme(
    legend.position      = c(0.85, 0.75),  # (x, y) in plot coordinates (0-1)
    legend.background    = element_rect(fill = "white", colour = "grey80",
                                        linewidth = 0.3),
    legend.key.size      = unit(0.8, "lines"),
    legend.text          = element_text(size = 10),
    plot.title           = element_text(size = 12, face = "bold"),
    axis.title           = element_text(size = 11),
    axis.text            = element_text(size = 10),
    panel.grid.minor     = element_blank()  # Remove minor grid lines
  )


# ── STEP 4: Build the log-scale inset ─────────────────────────────────────────
#
# CLINICAL PK CONCEPT — Why log scale for the terminal phase?
# ─────────────────────────────────────────────────────────────
# Drug elimination follows FIRST-ORDER KINETICS:
#   C(t) = C0 * e^(-ke * t)
#
# On a LINEAR scale: this produces a curve (hard to analyse)
# On a LOG scale:    ln(C) = ln(C0) - ke*t  → a STRAIGHT LINE
#
# The slope of that line = -ke (elimination rate constant)
# ke = ln(2) / t½    →    t½ = ln(2) / ke = 0.693 / ke
#
# Steeper slope = larger ke = faster elimination = shorter t½
# When you look at the log-scale inset:
#   Green line (with inducer) has a STEEPER SLOPE
#   This VISUALLY confirms that the inducer shortened t½
#   (confirmed numerically in Table 2: t½ GMR = 0.68)
# ─────────────────────────────────────────────────────────────

fig2_inset <- ggplot(
  data    = mean_conc %>% filter(TIME > 0, mean_conc > 0),
  # Exclude TIME == 0 (pre-dose = 0 ng/mL; log(0) is undefined)
  # Exclude any zero concentrations (cause log scale errors)
  mapping = aes(x = TIME, y = mean_conc, colour = ARM, group = ARM)
) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = ARM), size = 1.8) +

  scale_y_log10(
    # log10 scale on y-axis
    # label_log() formats tick labels as 10^x notation
    labels = label_log(digits = 2),
    breaks = c(0.03, 0.1, 0.3, 1.0)
  ) +
  scale_colour_manual(values = arm_colours, guide = "none") +
  # guide = "none" removes the colour legend from the inset (avoids duplication)
  scale_shape_manual(values = c("No_inducer" = 24, "With_inducer" = 21),
                     guide = "none") +

  scale_x_continuous(breaks = c(0, 4, 8, 12), limits = c(0, 12)) +

  labs(x = "Time (h)", y = "Conc (ng/mL)") +

  theme_bw(base_size = 8) +   # base_size = 8 makes everything smaller for inset
  theme(
    panel.grid.minor  = element_blank(),
    axis.title        = element_text(size = 7),
    axis.text         = element_text(size = 6),
    plot.background   = element_rect(fill = "white", colour = "grey70",
                                     linewidth = 0.3)
  )


# ── STEP 5: Combine main plot + inset using patchwork ─────────────────────────
# patchwork's inset_element() places one plot inside another at specified coords
# xmin/xmax/ymin/ymax are in normalised plot coordinates (0 = left/bottom, 1 = right/top)

fig2_final <- fig2_main +
  inset_element(
    fig2_inset,
    left   = 0.42,   # inset starts 42% from the left
    bottom = 0.35,   # inset starts 35% from the bottom
    right  = 0.98,   # inset ends 98% from the left
    top    = 0.98    # inset ends 98% from the bottom
  )

# Print to the RStudio Plots panel
print(fig2_final)
cat("Figure 2 rendered in Plots panel.\n")


# ── STEP 6: Save Figure 2 ─────────────────────────────────────────────────────
# Save as both PNG (for web/LinkedIn) and PDF (for publication/report)

ggsave(
  filename = here("output", "figures", "fig2_conc_time.png"),
  plot     = fig2_final,
  width    = 7,      # inches — standard single-column journal figure
  height   = 5,
  dpi      = 300     # 300 dpi = publication quality
)

ggsave(
  filename = here("output", "figures", "fig2_conc_time.pdf"),
  plot     = fig2_final,
  width    = 7,
  height   = 5
)

cat("✓ Saved: output/figures/fig2_conc_time.png\n")
cat("✓ Saved: output/figures/fig2_conc_time.pdf\n")


# ── STEP 7: Build Figure 3 — Spaghetti plot ───────────────────────────────────
#
# CLINICAL PK CONCEPT — What is a spaghetti plot and why do we need it?
# ───────────────────────────────────────────────────────────────────────
# Figure 2 shows the AVERAGE. But averages hide individual variability.
# Example: if one subject had almost NO induction effect while another
# had a huge effect, the mean would show a moderate effect.
# The spaghetti plot shows EACH INDIVIDUAL LINE — one line per subject.
#
# In Figure 3 we plot CL/F (apparent clearance) per subject, not concentration.
# Why CL/F and not AUC?
#   CL/F = Dose / AUC  (for oral dosing with unknown bioavailability F)
#   Higher CL/F = drug cleared faster = more enzyme activity
#   Inducer INCREASES CYP3A4 activity → INCREASES CL/F
#
# So we expect: green dots (with inducer) should be HIGHER than orange dots
# for every subject. A line connecting each subject's two dots shows
# the DIRECTION and MAGNITUDE of the DDI for that person.
#
# If all lines go UP from left to right: consistent induction across subjects
# If one line goes DOWN: that subject may be a CYP3A5 expresser or outlier
# This is clinically meaningful — individual variability is a key PK finding
# ───────────────────────────────────────────────────────────────────────────

# We need CL/F per subject per arm.
# CL/F is computed in 04_nca_ncappc.R — for now we SIMULATE it from
# the known AUClast values to be able to build the plot today.
#
# From Table 2 of Dunvald 2022:
#   Without inducer: median CL/F = 3236.67 L/h
#   With inducer:    median CL/F = 4780.99 L/h
#   Individual values shown in Figure 3 range approximately 2000–6000 L/h

set.seed(42)  # Same seed = same random values every run

clf_data <- tibble(
  ID  = rep(1:12, 2),                           # Each subject appears twice
  ARM = c(rep("No_inducer", 12), rep("With_inducer", 12)),
  # Simulate individual CL/F values consistent with Table 2 summary stats
  CLF = c(
    # No inducer: mean ~3237, range ~2700-3500 from Table 2 IQR
    round(rnorm(12, mean = 3237, sd = 350), 0),
    # With inducer: mean ~4781, range ~3872-5131 from Table 2 IQR
    round(rnorm(12, mean = 4781, sd = 500), 0)
  )
)

# Clamp to physiologically plausible range for oral midazolam CL/F
clf_data <- clf_data %>%
  mutate(CLF = pmax(2000, pmin(6200, CLF)))

cat("\n=== CL/F summary ===\n")
clf_data %>%
  group_by(ARM) %>%
  summarise(
    median_CLF = round(median(CLF), 0),
    Q1 = round(quantile(CLF, 0.25), 0),
    Q3 = round(quantile(CLF, 0.75), 0)
  ) %>%
  print()


# ── Build the spaghetti plot ──────────────────────────────────────────────────

# Create a labeller for cleaner x-axis text
arm_x_labels <- c("No_inducer" = "No inducer", "With_inducer" = "Inducer")

fig3 <- ggplot(
  data    = clf_data,
  mapping = aes(
    x     = ARM,     # x-axis: the two treatment arms (categorical)
    y     = CLF,     # y-axis: individual CL/F values
    group = ID       # group = ID means: draw one line per subject
  )
) +

  # ── Geom 1: Lines connecting each subject's two values ──
  # Because group = ID above, geom_line() draws one line per subject
  # This is what makes it a "spaghetti" plot — the lines look like pasta
  geom_line(
    colour    = "grey50",  # Grey lines for individual subjects
    linewidth = 0.6,
    alpha     = 0.7        # Slightly transparent so overlapping lines are visible
  ) +

  # ── Geom 2: Points at each measurement ──
  # colour = ARM means: colour depends on which arm the observation is from
  geom_point(
    aes(colour = ARM),
    size  = 3,
    alpha = 0.9
  ) +

  scale_colour_manual(
    values = arm_colours,
    labels = arm_labels,
    name   = NULL
  ) +

  # ── Customise x-axis labels ──
  scale_x_discrete(labels = arm_x_labels) +

  # ── Y-axis ──
  scale_y_continuous(
    limits = c(1800, 6500),
    breaks = seq(2000, 6000, by = 1000),
    labels = scales::comma    # Formats 3000 as "3,000" not "3000"
  ) +

  labs(
    title   = "Individual midazolam apparent clearance (CL/F)",
    x       = NULL,                          # No x-axis title needed
    y       = expression("CL/F (L h"^{-1}~")"),
    # expression() allows mathematical notation: L h^-1
    caption = "Each line represents one subject (n = 12)."
  ) +

  theme_bw() +
  theme(
    legend.position  = "none",               # No legend needed (x-axis is self-explanatory)
    plot.title       = element_text(size = 12, face = "bold"),
    axis.title.y     = element_text(size = 11),
    axis.text        = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()     # Remove vertical grid lines
  )

print(fig3)
cat("Figure 3 rendered in Plots panel.\n")


# ── STEP 8: Save Figure 3 ────────────────────────────────────────────────────

ggsave(
  filename = here("output", "figures", "fig3_spaghetti_clf.png"),
  plot     = fig3,
  width    = 5,
  height   = 5,
  dpi      = 300
)

ggsave(
  filename = here("output", "figures", "fig3_spaghetti_clf.pdf"),
  plot     = fig3,
  width    = 5,
  height   = 5
)

cat("✓ Saved: output/figures/fig3_spaghetti_clf.png\n")
cat("✓ Saved: output/figures/fig3_spaghetti_clf.pdf\n")


# ── STEP 9: What to look for in your plots ────────────────────────────────────
cat("\n=== WHAT TO VERIFY IN YOUR PLOTS ===\n")
cat("Figure 2 — Concentration-time curve:\n")
cat("  ✓ Orange line (No inducer) should be HIGHER than green (Inducer)\n")
cat("  ✓ Both lines peak around t = 2h (Tmax from Table 2)\n")
cat("  ✓ Error bars show SD — some overlap is normal\n")
cat("  ✓ Inset (log scale): green line has STEEPER slope = shorter t½\n")
cat("  ✓ At t=0: both lines start at 0 (pre-dose = no drug)\n\n")
cat("Figure 3 — Spaghetti plot:\n")
cat("  ✓ EVERY line should go UP from left (No inducer) to right (Inducer)\n")
cat("  ✓ CL/F range: ~2000-5000 L/h without inducer, ~3500-6200 with\n")
cat("  ✓ Some subjects show a larger increase than others (inter-subject variability)\n")
cat("  ✓ No lines should go DOWN — that would suggest an inhibitor, not inducer\n")

cat("\n=== 03_pk_plots.R COMPLETE ===\n")
cat("Next: open R/04_nca_ncappc.R\n")
