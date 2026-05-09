#!/usr/bin/env Rscript
# =============================================================================
# plot_band_pi_panel.R
# =============================================================================
# Two-panel PDF for one inversion:
#   Top    : pi11 / pi22 / pi12 windowed across the interval
#   Bottom : mean per-sample tP for Homo_1 / Het / Homo_2
# Inversion interval is shaded. Mean_tP panel is omitted if no file exists.
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

opt_list <- list(
  make_option("--inv",         type = "character"),
  make_option("--chr",         type = "character"),
  make_option("--start",       type = "integer"),
  make_option("--end",         type = "integer"),
  make_option("--window-tsv",  type = "character"),
  make_option("--mean-tp-tsv", type = "character", default = NULL),
  make_option("--out-pdf",     type = "character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(
  !is.null(opt$inv), !is.null(opt$chr),
  !is.null(opt$start), !is.null(opt$end),
  !is.null(opt$`window-tsv`), !is.null(opt$`out-pdf`)
)

w <- fread(opt$`window-tsv`, sep = "\t", header = TRUE)
basis_lab <- if ("site_basis" %in% colnames(w)) {
  paste0(" (site_basis = ", w$site_basis[1], ")")
} else ""

w_long <- melt(
  w,
  id.vars = c("chromo", "win_center"),
  measure.vars = c("mean_pi11", "mean_pi22", "mean_pi12"),
  variable.name = "metric", value.name = "pi"
)
w_long[, metric := factor(
  metric,
  levels = c("mean_pi11", "mean_pi22", "mean_pi12"),
  labels = c("pi11", "pi22", "pi12")
)]

band_cols <- c(pi11 = "#1F77B4", pi22 = "#D62728", pi12 = "#2CA02C")

p_top <- ggplot(w_long, aes(x = win_center / 1e6, y = pi,
                            colour = metric, group = metric)) +
  annotate("rect", xmin = opt$start / 1e6, xmax = opt$end / 1e6,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = band_cols, name = NULL) +
  labs(
    title = sprintf("%s — band-level pi (pre-polarization)%s",
                    opt$inv, basis_lab),
    subtitle = sprintf("%s : %d–%d", opt$chr, opt$start, opt$end),
    x = sprintf("Position on %s (Mb)", opt$chr),
    y = expression("pi"[11]~","~"pi"[22]~","~"pi"[12])
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.subtitle = element_text(colour = "grey40", size = 9))

# Bottom panel: mean per-sample tP, if file exists
have_tp <- !is.null(opt$`mean-tp-tsv`) && file.exists(opt$`mean-tp-tsv`)
if (have_tp) {
  tp <- fread(opt$`mean-tp-tsv`, sep = "\t", header = TRUE)
  if (nrow(tp) > 0L) {
    tp[, group := factor(group, levels = c("Homo_1", "Het", "Homo_2"))]
    grp_cols <- c(Homo_1 = "#1F77B4", Het = "#7F7F7F", Homo_2 = "#D62728")
    p_bot <- ggplot(tp, aes(x = WinCenter / 1e6, y = mean_tP,
                            colour = group, group = group,
                            ymin = mean_tP - se_tP, ymax = mean_tP + se_tP,
                            fill = group)) +
      annotate("rect", xmin = opt$start / 1e6, xmax = opt$end / 1e6,
               ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
      geom_ribbon(alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.6) +
      scale_colour_manual(values = grp_cols, name = NULL) +
      scale_fill_manual(values = grp_cols, guide = "none") +
      labs(
        x = sprintf("Position on %s (Mb)", opt$chr),
        y = "mean per-sample tP (pestPG, per-site)",
        caption = "Mean of within-individual diversity by karyotype group; not formal group pi."
      ) +
      theme_classic(base_size = 11) +
      theme(legend.position = "top",
            plot.caption = element_text(colour = "grey40", size = 8, hjust = 0))
    p <- p_top / p_bot + plot_layout(heights = c(1, 1))
  } else {
    p <- p_top
  }
} else {
  p <- p_top
}

dir.create(dirname(opt$`out-pdf`), recursive = TRUE, showWarnings = FALSE)
ggsave(opt$`out-pdf`, p, width = 9, height = if (have_tp) 7 else 4, units = "in")
message(sprintf("[plot_band_pi_panel] wrote %s", opt$`out-pdf`))
