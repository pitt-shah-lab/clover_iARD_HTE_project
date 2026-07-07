# ============================================================================
# 09_build_iard_plot.R
#
# Generates the individual absolute risk difference (iARD) figure for
# the CLOVERS confirmatory analysis, similar to Talisa et al. (CCX)
# Figure showing iARDs by subtype for ADRENAL/EGDT.
#
# Left panel:  Group-level ATE by tertile (box/diamond)
# Right panel: Individual patient iARDs by tertile (strip/swarm)
#
# One figure per model (clinical only, clinical + biomarker).
#
# Input:
#   outputs/10_confirmatory_analysis/confirmatory_preds_clinical.csv
#   outputs/10_confirmatory_analysis/confirmatory_preds_with_biomarker.csv
#
# Output:
#   outputs/09_build_cate_outputs/iard_plot_clinical.png
#   outputs/09_build_cate_outputs/iard_plot_with_biomarker.png
#   outputs/09_build_cate_outputs/iard_plot_combined.png
# ============================================================================

this_script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
source(file.path(dirname(dirname(this_script)), "config", "config.R"))
out_dir <- make_output_subdir("09_build_cate_outputs")

library(ggplot2)

# ---- Load predictions ----
preds_clin <- read.csv(file.path(OUTPUTS_DIR,
  "10_confirmatory_analysis", "confirmatory_preds_clinical.csv"))
preds_bio  <- read.csv(file.path(OUTPUTS_DIR,
  "10_confirmatory_analysis", "confirmatory_preds_with_biomarker.csv"))

cat("Clinical predictions:", nrow(preds_clin), "patients\n")
cat("Biomarker predictions:", nrow(preds_bio), "patients\n")

# ---- Assign tertile subgroups ----
assign_subgroup <- function(df) {
  q33 <- quantile(df$score, 1/3)
  q67 <- quantile(df$score, 2/3)
  df$subgroup <- ifelse(df$score < q33, "Benefit",
                 ifelse(df$score > q67, "Harm", "Indeterminate"))
  df$subgroup <- factor(df$subgroup, levels = c("Benefit", "Indeterminate", "Harm"))
  df
}

preds_clin <- assign_subgroup(preds_clin)
preds_bio  <- assign_subgroup(preds_bio)

# ---- Compute group-level ATEs ----
compute_group_ate <- function(df) {
  do.call(rbind, lapply(levels(df$subgroup), function(grp) {
    sub <- df[df$subgroup == grp, ]
    rate_w1 <- mean(sub$y[sub$w == 1])
    rate_w0 <- mean(sub$y[sub$w == 0])
    ate <- rate_w1 - rate_w0
    n <- nrow(sub)
    se <- sqrt(rate_w1*(1-rate_w1)/sum(sub$w==1) + rate_w0*(1-rate_w0)/sum(sub$w==0))
    data.frame(subgroup = grp, ate = ate, se = se, n = n,
               lo = ate - 1.96*se, hi = ate + 1.96*se)
  }))
}

ate_clin <- compute_group_ate(preds_clin)
ate_bio  <- compute_group_ate(preds_bio)

# ---- Color palette ----
cols <- c("Benefit" = "#2E86AB", "Indeterminate" = "#A23B72", "Harm" = "#E8A838")

# ---- Build the combined figure ----
# Convert iARD scores to percentage scale
preds_clin$iard_pct <- preds_clin$score * 100
preds_bio$iard_pct  <- preds_bio$score * 100
ate_clin$ate_pct    <- ate_clin$ate * 100
ate_clin$lo_pct     <- ate_clin$lo * 100
ate_clin$hi_pct     <- ate_clin$hi * 100
ate_bio$ate_pct     <- ate_bio$ate * 100
ate_bio$lo_pct      <- ate_bio$lo * 100
ate_bio$hi_pct      <- ate_bio$hi * 100

# ---- PANEL FUNCTION ----
make_iard_plot <- function(preds_df, ate_df, title, filename) {

  # Jitter x positions for strip plot
  set.seed(42)
  preds_df$x_jitter <- as.numeric(preds_df$subgroup) + runif(nrow(preds_df), -0.25, 0.25)

  p <- ggplot() +

    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +

    # Individual iARDs (strip plot)
    geom_point(data = preds_df,
               aes(x = x_jitter, y = iard_pct, color = subgroup),
               size = 1.2, alpha = 0.5, shape = 18) +

    # Group ATE diamond + error bars
    geom_errorbar(data = ate_df,
                  aes(x = as.numeric(factor(subgroup, levels = c("Benefit","Indeterminate","Harm"))),
                      ymin = lo_pct, ymax = hi_pct),
                  width = 0.15, linewidth = 0.8, color = "black") +
    geom_point(data = ate_df,
               aes(x = as.numeric(factor(subgroup, levels = c("Benefit","Indeterminate","Harm"))),
                   y = ate_pct, fill = subgroup),
               shape = 23, size = 4, color = "black", stroke = 0.8) +

    scale_color_manual(values = cols, guide = "none") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = c("Benefit", "Indeterminate", "Harm")) +

    labs(
      title = title,
      x = "CATE-predicted subgroup (tertile)",
      y = "Individual absolute risk difference (%)\n(liberal \u2212 restrictive)"
    ) +

    # Annotations
    annotate("text", x = 3.6, y = max(preds_df$iard_pct) * 0.8,
             label = "Increasing benefit\nfrom liberal",
             hjust = 0.5, size = 3, color = "grey40") +
    annotate("segment", x = 3.6, xend = 3.6,
             y = max(preds_df$iard_pct) * 0.5,
             yend = max(preds_df$iard_pct) * 0.7,
             arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
    annotate("text", x = 3.6, y = min(preds_df$iard_pct) * 0.8,
             label = "Increasing harm\nfrom liberal",
             hjust = 0.5, size = 3, color = "grey40") +
    annotate("segment", x = 3.6, xend = 3.6,
             y = min(preds_df$iard_pct) * 0.5,
             yend = min(preds_df$iard_pct) * 0.7,
             arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +

    coord_cartesian(xlim = c(0.5, 4.0)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_dir, filename), p, width = 7, height = 6, dpi = 300)
  cat("Wrote:", filename, "\n")
  p
}

p1 <- make_iard_plot(preds_clin, ate_clin,
                      "iARD by CATE subgroup — Clinical only (n = 670 validation)",
                      "iard_plot_clinical.png")

p2 <- make_iard_plot(preds_bio, ate_bio,
                      "iARD by CATE subgroup — Clinical + biomarker (n = 670 validation)",
                      "iard_plot_with_biomarker.png")

# ---- Combined side-by-side ----
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_combined <- p1 + p2 + plot_layout(ncol = 2) +
    plot_annotation(title = "Individual absolute risk differences by CATE-predicted subgroup",
                    subtitle = "Diamonds = group ATE with 95% CI. Each point = one patient's predicted iARD.",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  ggsave(file.path(out_dir, "iard_plot_combined.png"), p_combined,
         width = 14, height = 6, dpi = 300)
  cat("Wrote: iard_plot_combined.png\n")
} else {
  cat("Install 'patchwork' for the combined side-by-side plot.\n")
}

cat("\nDone. All plots in:", out_dir, "\n")