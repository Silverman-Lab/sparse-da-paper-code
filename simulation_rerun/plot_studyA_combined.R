################################################################################
## plot_studyA_combined.R -- the published-style Study A figure:
##   3 rows (FDR, TPR, F0.5) x 4 columns (N, D, lambda, pi), one shared legend
##   on the right listing every method.
##
## Colours, line types and method ORDER are copied from the original
##   simulation/src/06072025/Results combined updated.R
## specifically from the LAST definition of each in that file (method_colors at
## line 72, line_types at line 89, neworder at line 122) -- it redefines
## method_colors four times and line_types twice, and only the final one is the
## version that produced the published figures. Do not copy from the earlier
## blocks; they are superseded drafts kept in place.
##
## The one rename: the original calls the knockoff method "CKF", this pipeline
## calls it "MRCKF". Everything else matches by name.
##
## Reads results/main/<tag>/aggregated_metrics.csv, written by aggregate_main.R.
## Usage: Rscript plot_studyA_combined.R [tag]
################################################################################

.bootstrap_config <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) {
    cfg <- file.path(p, "simulation_rerun", "config.R")
    if (file.exists(cfg)) return(cfg)
    p <- dirname(p)
  }
  stop("could not locate simulation_rerun/config.R above ", getwd())
}
source(.bootstrap_config())

suppressWarnings(suppressMessages({
  library(ggplot2); library(patchwork); library(cowplot)
}))

args <- commandArgs(trailingOnly = TRUE)
tag  <- if (length(args) >= 1) args[1] else MAIN_TAG
## Optional 2nd arg: a suffix such as "_width", selecting an aggregate produced
## under an alternative scale rule (see redecide_scale_rule.R). Output figure
## names carry the same suffix, so the two rules' figures coexist.
##
## Defaults to the SCENARIO suffix the tag carries ("_pos20" for the
## supplementary 20/80 sign split, "" for the main design), which is what
## aggregate_main.R named its outputs with -- so plotting a supplementary tree
## reads the right table and cannot overwrite the main figure.
sfx <- if (length(args) >= 2) args[2] else scenario_sfx_of(tag)
agg_file <- file.path(RES_MAIN, tag, paste0("aggregated_metrics", sfx, ".csv"))
if (!file.exists(agg_file)) stop("no aggregated_metrics.csv at ", agg_file,
                                 " -- run aggregate_main.R first")
agg <- read.csv(agg_file, stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## Plotted pi range
## ---------------------------------------------------------------------------
## pi = 0.75 IS simulated and stays in aggregated_metrics.csv and every table --
## this cap is presentation only. It is excluded from the figure because at three
## quarters of taxa differentially abundant the scale model is not identifiable,
## so the point measures the breakdown rather than the method: the hybrid rule
## selects Laplace in 50/50 replicates (margin -3.16) and TPR falls to 0.104
## while the bootstrap arm still holds 0.418. Plotted, it compresses the y-range
## of the whole pi column around a setting no real study occupies.
PI_PLOT_MAX <- 0.65
.drop <- agg$sweep == "overPi" & agg$lever_value > PI_PLOT_MAX
if (any(.drop)) {
  message("pi > ", PI_PLOT_MAX, " simulated but omitted from the figure: ",
          paste(sort(unique(agg$lever_value[.drop])), collapse = ", "))
  agg <- agg[!.drop, ]
}

## ---------------------------------------------------------------------------
## Conventions -- colours, line types, legend order
## ---------------------------------------------------------------------------
## Moved to src/plot_conventions.R so plot_runtime.R draws the same methods the
## same way by construction rather than by a second copy of the values.
source(file.path(HERE, "src", "plot_conventions.R"))

## agg keeps EVERY method; each figure below selects its own subset. The main
## figure shows the HYBRID SSRV rule only -- the Laplace-only and bootstrap-only
## arms would otherwise draw three near-identical red lines on top of it -- and
## the decomposition figure shows exactly those three arms and nothing else.
missing <- setdiff(neworder, unique(as.character(agg$method)))
if (length(missing)) message("methods absent from this result set: ",
                             paste(missing, collapse = ", "))

## ---------------------------------------------------------------------------
## Decomposition figure conventions
## ---------------------------------------------------------------------------
## The hybrid rule selected Laplace in ~99% of replicates, so those two curves
## coincide almost everywhere. Drawing hybrid THICK and Laplace as a thin dash
## ON TOP makes the coincidence legible -- the dashes ride along the red line --
## instead of one silently hiding the other and looking like a missing series.
decomp_methods <- c("Sparse_SSRV", "Sparse_SSRV_L", "Sparse_SSRV_B")
decomp_labels  <- c(Sparse_SSRV   = "Hybrid",
                    Sparse_SSRV_L = "Laplace only",
                    Sparse_SSRV_B = "Bootstrap only")
decomp_colors  <- c("Hybrid"         = "red",
                    "Laplace only"   = "black",
                    "Bootstrap only" = "#0072B2")
decomp_lines   <- c("Hybrid"         = "solid",
                    "Laplace only"   = "dashed",
                    "Bootstrap only" = "solid")
decomp_widths  <- c("Hybrid"         = 1.4,
                    "Laplace only"   = 0.6,
                    "Bootstrap only" = 0.9)

## ---------------------------------------------------------------------------
## Panels
## ---------------------------------------------------------------------------
## Column order: D, N, lambda, pi. The names are the panel x-axis labels and are
## deliberately plain English rather than the symbol the code sweeps -- a reader
## should not have to know that D is the feature count or that pi is the
## proportion of true associations to read the figure.
##
## BOTH figures below are built from this vector, so the order and the labels
## cannot differ between the all-methods figure and the scale-model one.
##
## LOGX matches aggregate_main.R -- the three count-like levers are plotted on a
## log axis, pi linearly. LOGX and REVX are looked up BY NAME, so their order is
## cosmetic; they are kept in step with SWEEPS so a reordering reads as one edit.
## From config.R, so this figure and aggregate_main.R's per-sweep figures carry
## the same axis names. The order of SWEEP_LABEL is the panel order.
SWEEPS <- SWEEP_LABEL
LOGX   <- c(overD = TRUE, overN = TRUE, overLambda = TRUE, overPi = FALSE)
## Reversed x on the pi panel, matching the older manuscript figure, which used
## xlim(c(max(pi), min(pi))). pi decreasing left-to-right puts the hardest
## setting at the left edge, so every curve is read in the direction of
## increasing information. scale_x_reverse(), not xlim() -- xlim() DROPS rows
## outside its range before the line is drawn, which would silently break a
## series; reversing the scale keeps every point.
REVX   <- c(overD = FALSE, overN = FALSE, overLambda = FALSE, overPi = TRUE)

## Drop any sweep with fewer than two settings, matching aggregate_main.R's rule
## and for the same reason: expand_sweeps() attaches the baseline to ALL FOUR
## sweeps, so a sweep that has not been run still holds exactly one point. Drawn,
## that is a lone dot in a panel captioned like a curve -- indistinguishable at a
## glance from a real result that happens to be flat. Dropping the column says
## "not run" unambiguously.
##
## The figure therefore has as many columns as there are sweeps with data, and
## the panel grid below sizes itself from that rather than from a fixed 4.
.n_lev <- tapply(agg$lever_value, agg$sweep, function(v) length(unique(v)))
.have  <- names(.n_lev)[.n_lev >= 2]
if (length(setdiff(names(SWEEPS), .have)))
  message("sweeps with <2 settings, omitted from the figure: ",
          paste(setdiff(names(SWEEPS), .have), collapse = ", "))
SWEEPS <- SWEEPS[names(SWEEPS) %in% .have]
LOGX   <- LOGX[names(LOGX) %in% .have]
REVX   <- REVX[names(REVX) %in% .have]
if (!length(SWEEPS)) stop("no sweep has two or more settings yet -- nothing to plot")
ROWS   <- c(fdr = "FDR", tdr = "TPR", f0.5 = "F0.5 Score")

## Fixed 0-1 y-axis on every panel, so the two figures are directly comparable
## and a value can be read off without checking the axis first. coord_cartesian
## (not scale_y_continuous limits) because it CLIPS the view rather than dropping
## rows -- with limits=, any point outside the range is deleted before the line
## is drawn, which silently breaks a series instead of running it off the edge.
## Set to NULL for free per-panel scales, which is what makes small differences
## (D = 20 in the decomposition figure) legible.
YLIM <- c(0, 1)

## Size of the per-panel x-axis titles (the sweep names).
AXIS_X_SIZE <- 15

panel <- function(sw, metric, show_x, dat = agg, cols = method_colors,
                  ltys = line_types, lwds = NULL, lab = "Method") {
  d <- dat[dat$sweep == sw, ]
  if (!nrow(d)) return(NULL)
  ggplot(d, aes(x = lever_value, y = .data[[metric]],
                group = method, color = method, linetype = method)) +
    ## Nominal-level reference on the FDR panels only. FIRST in the layer order,
    ## so every method curve draws on top of it rather than the other way round.
    ##
    ## grey40, not the grey80 used for the scale-free methods: five methods
    ## (DESeq2, baySeq, edgeR, limma, MRCKF) are drawn in grey80, and a grey80
    ## rule sitting among them would read as a sixth flat method. Darker and
    ## thinner separates it; horizontal-and-dashed does the rest.
    { if (metric == "fdr")
        geom_hline(yintercept = PVAL, colour = "grey40",
                   linetype = "dashed", linewidth = 0.4) } +
    { if (is.null(lwds)) geom_line(linewidth = 0.7)
      else list(aes(linewidth = method), geom_line(),
                scale_linewidth_manual(name = lab, values = lwds, guide = "none")) } +
    geom_point(size = 1.1) +
    scale_color_manual(name = lab, values = cols, drop = FALSE) +
    scale_linetype_manual(name = lab, values = ltys, drop = FALSE) +
    { if (isTRUE(LOGX[[sw]])) scale_x_log10()
      else if (isTRUE(REVX[[sw]])) scale_x_reverse()
      else NULL } +
    { if (!is.null(YLIM)) coord_cartesian(ylim = YLIM) else NULL } +
    labs(x = if (show_x) SWEEPS[[sw]] else NULL, y = ROWS[[metric]]) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.title = element_text(size = 12),
          ## The x titles name the sweep and are the figure's captions in
          ## practice, so they carry more weight than the metric names on y.
          axis.title.x = element_text(size = AXIS_X_SIZE,
                                      margin = margin(t = 6)))
}

## One panel carries the legend; it is extracted and placed on the right so the
## twelve panels stay the same size. override.aes re-pairs colour with line type
## in the key -- without it the keys draw as solid colour swatches, and in the
## main figure the two colour-sharing pairs become indistinguishable.
build_figure <- function(dat, cols, ltys, order_levels, lwds = NULL,
                         lab = "Method", out_stem, legend_rel = 0.75) {
  dat <- dat[dat$method %in% order_levels, ]
  dat$method <- factor(as.character(dat$method), levels = order_levels)

  pan <- function(sw, m, show_x)
    panel(sw, m, show_x, dat = dat, cols = cols, ltys = ltys, lwds = lwds, lab = lab)

  ## The FIRST sweep, whichever it is -- naming one hardcoded sweep breaks the
  ## legend the moment SWEEPS is reordered or that sweep is dropped for having
  ## too few settings.
  legend_src <- pan(names(SWEEPS)[1], "fdr", TRUE) +
    theme(legend.position = "right",
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 11),
          legend.key.size = unit(1.4, "lines")) +
    guides(color = guide_legend(override.aes = list(
             linetype = unname(ltys[order_levels]))),
           linetype = "none")
  legend <- cowplot::get_legend(legend_src)

  plots <- list()
  for (m in names(ROWS))
    for (sw in names(SWEEPS))
      plots[[length(plots) + 1]] <- pan(sw, m, show_x = (m == "f0.5"))
  plots <- Filter(Negate(is.null), plots)

  ## ncol = length(SWEEPS), not 4: with a sweep omitted, a hardcoded 4 wraps the
  ## 3 x 3 panels into a ragged 4-wide block and silently puts FDR and TPR panels
  ## on the same row.
  ncol_sw <- length(SWEEPS)
  final_plot <- plot_grid(wrap_plots(plots, ncol = ncol_sw), legend,
                          ncol = 2, rel_widths = c(2 * ncol_sw, legend_rel)) +
    ## An explicit white background: ggsave writes a TRANSPARENT one otherwise,
    ## which most viewers render as black.
    theme(plot.background = element_rect(fill = "white", colour = NA))

  dir.create(RES_FIG, recursive = TRUE, showWarnings = FALSE)
  for (ext in c("pdf", "png")) {
    f <- file.path(RES_FIG, paste0(out_stem, ".", ext))
    ## Width follows the column count so panels keep their aspect ratio whether
    ## the figure has two sweeps or four (16in was sized for 4).
    ggsave(f, final_plot, width = 4 + 3 * ncol_sw, height = 9,
           dpi = if (ext == "png") 150 else 300, bg = "white")
    message("wrote ", f)
  }
  invisible(final_plot)
}

## ---- figure 1: all methods, hybrid SSRV only ------------------------------
build_figure(agg, method_colors, line_types, neworder,
             out_stem = paste0("studyA_combined_fdr_tpr_f0.5", sfx))

## ---- figure 2: SSRV only, the three scale models --------------------------
d2 <- agg[agg$method %in% decomp_methods, ]
d2$method <- decomp_labels[as.character(d2$method)]
build_figure(d2, decomp_colors, decomp_lines, unname(decomp_labels),
             lwds = decomp_widths, lab = "Scale model",
             out_stem = paste0("studyA_decomp_scale_model", sfx), legend_rel = 0.9)
