################################################################################
## plot_runtime.R -- computation time by method, one figure per sweep.
##
##   results/figures/runtime_over_{D,N,lambda}<sfx>.pdf / .png
##   results/main/<tag>/runtime_over_{D,N,lambda}<sfx>.csv   the numbers behind them
##
## Generalised from an overD-only version on 2026-09-05. The overD filenames are
## Renamed timing_* -> runtime_* on 2026-09-05; the superseded timing_* files
## were deleted rather than left beside them as stale near-duplicates.
##
## Usage: Rscript plot_runtime.R [tag] [suffix]
##
## Replaces simulation/src/06072025/running time comparison.R, which cannot be
## re-pointed at this pipeline:
##   * it hardcodes /storage/work/... paths;
##   * it reads "final execution time.csv"; the worker here writes
##     "execution time.csv" (run_methods_rerun.R);
##   * its grid is D in {20,50,100,200,350,500,700,1000} at lambda = 2500*D,
##     against this study's {20,50,100,300,500,700,1000} at lambda ~ 1666.67*D.
##     THE X-VALUES ARE THEREFORE NOT THE PUBLISHED ONES and the curves are not
##     comparable to the old figure point for point. Say so in the caption.
##   * it indexes methods POSITIONALLY (runtime$X[3:12], elapsed[3:13]), so a
##     method that fails or is added shifts every label silently. Keyed by name
##     here, for the same reason run_all_methods() keys truth by taxon name.
##
## Colours and line types come from src/plot_conventions.R -- the same object the
## main Study A figure uses, so a method is the same colour in both figures by
## construction rather than by a second copy of the palette.
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
source(file.path(HERE, "src", "plot_conventions.R"))

suppressWarnings(suppressMessages({ library(ggplot2) }))

args <- commandArgs(trailingOnly = TRUE)
tag  <- if (length(args) >= 1) args[1] else MAIN_TAG
## Defaults to the scenario suffix the TAG carries ("_pos20" for the
## supplementary 20/80 sign split, "" for the main design), matching
## aggregate_main.R -- results/figures/ is shared by every tag, so the two
## scenarios' runtime figures must not land on the same filename.
sfx  <- if (length(args) >= 2) args[2] else scenario_sfx_of(tag)
TAG_DIR <- file.path(RES_MAIN, tag)
if (!dir.exists(TAG_DIR)) stop("no results at ", TAG_DIR)

## ---------------------------------------------------------------------------
## execution time.csv row name -> the palette's method name
## ---------------------------------------------------------------------------
## Sparse SSRV is `sparse_resample` + `sparse`: the Dirichlet array is not
## optional, so charging only the fit would understate it against methods that
## take raw counts. Same treatment as the published figure's
## sum(runtime$elapsed[1:2]).
##
## `sparse_bw_*` rows -- the supplementary nine-variant grid -- are DELIBERATELY
## excluded. They are three fits of one method, so including them would compare
## a sensitivity analysis against other methods' single fits. Their cost is
## reported separately below.
SPARSE_PARTS <- c("sparse_resample", "sparse")
TIME_MAP <- c(
  aldex_clr = "ALDEx2_clr", aldex_def = "ALDEx2_def",
  deseq2    = "DESeq2",     ancombc   = "ANCOM-BC",
  ancom     = "ANCOM",      ancombc2  = "ANCOM-BC2",
  linda     = "LinDA",      edgeR     = "edgeR",
  limma     = "limma",      baySeq    = "baySeq",
  mrckf     = "MRCKF")

dirs <- list.dirs(TAG_DIR, recursive = FALSE)
dirs <- dirs[grepl("^results_", basename(dirs))]
if (!length(dirs)) stop("no replicate directories under ", TAG_DIR)

read_time <- function(d) {
  f <- file.path(d, "execution time.csv")
  if (!file.exists(f)) return(NULL)
  tt <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(tt) || !nrow(tt)) return(NULL)
  info <- tryCatch(read.csv(file.path(d, "replicate_info.csv"), stringsAsFactors = FALSE),
                   error = function(e) NULL)
  if (is.null(info)) return(NULL)

  el <- setNames(tt$elapsed, tt$X)
  rows <- data.frame(method = unname(TIME_MAP[names(el)]),
                     secs   = unname(el), stringsAsFactors = FALSE)
  rows <- rows[!is.na(rows$method), , drop = FALSE]

  sp <- el[intersect(SPARSE_PARTS, names(el))]
  if (length(sp) == length(SPARSE_PARTS))
    rows <- rbind(rows, data.frame(method = "Sparse_SSRV", secs = sum(sp)))

  ## What the supplementary grid costs, carried alongside rather than folded in.
  bwv <- el[grepl("^sparse_bw_", names(el))]
  attr(rows, "secs_variants") <- if (length(bwv)) sum(bwv) else 0

  cbind(rows, N = info$N, D = info$D, lambda = info$lambda, pi = info$pi,
        seed = info$seed, dir = basename(d),
        secs_variants = attr(rows, "secs_variants"))
}

TT <- do.call(rbind, Filter(Negate(is.null), lapply(dirs, read_time)))
if (is.null(TT) || !nrow(TT)) stop("no execution time.csv found under ", TAG_DIR)

## --- MiMIX: not in execution time.csv, it is a separate Julia fit -----------
## run_main.R records it per replicate in the run manifest.
logs <- list.files(TAG_DIR, "^run_log_.*\\.csv$", full.names = TRUE)
if (length(logs)) {
  lg <- do.call(rbind, lapply(logs, function(f)
    tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)))
  lg <- lg[!is.na(lg$secs_mimix), , drop = FALSE]
  if (nrow(lg)) {
    ## One row per dir: a replicate re-run appears in more than one manifest, and
    ## the last write is the one whose results are on disk.
    lg <- lg[!duplicated(lg$dir, fromLast = TRUE), , drop = FALSE]
    TT <- rbind(TT, data.frame(method = "MiMIX", secs = lg$secs_mimix,
                               N = lg$N, D = lg$D, lambda = lg$lambda, pi = lg$pi,
                               seed = lg$rep, dir = lg$dir, secs_variants = 0))
  }
} else message("no run_log_*.csv -- MiMIX will be absent from the figure")

## ---------------------------------------------------------------------------
## One figure per sweep
## ---------------------------------------------------------------------------
## Sweep membership comes from main_sweeps_of() for each sweep, exactly as the
## overD-only version did, so a timing figure cannot drift from the study's own
## definition of that sweep -- including the baseline, which belongs to all four.
##
## All three levers are count-like and span 1.5-3 decades, so all three x axes
## are log, matching the metric figures. pi is deliberately absent: it does not
## change the size of the problem, so a pi timing curve would be a flat line per
## method and say nothing the other three do not.
## Axis names and x scaling copied from the previous manuscript's
##   simulation/src/06072025/running time comparison.R
## which titled them "Number of Features" and "Sample Size" against
## "Average Running Time (s)", and used a plain LINEAR x (it called no scale_x_*).
##
## overLambda is the exception and is flagged rather than assumed: that sweep has
## no counterpart in the old script, and its levers span 1e4-1e7, so a linear x
## would stack five of the six settings against the left axis. It keeps log10.
TIMING_SWEEPS <- list(
  overD      = list(var = "D",      lab = "Number of Features",
                    stem = "runtime_over_D",      xlog = FALSE),
  overN      = list(var = "N",      lab = "Sample Size",
                    stem = "runtime_over_N",      xlog = FALSE),
  overLambda = list(var = "lambda", lab = "Sequencing Depth",
                    stem = "runtime_over_lambda", xlog = TRUE)
)

agg_timing <- function(sw, spec) {
  keep <- vapply(seq_len(nrow(TT)), function(i)
    sw %in% main_sweeps_of(TT$N[i], TT$D[i], TT$lambda[i], TT$pi[i]), logical(1))
  TD <- TT[keep, , drop = FALSE]
  if (!nrow(TD)) { message("no replicates on ", sw, " -- skipped"); return(invisible(NULL)) }

  ## Same rule as aggregate_main.R: a sweep with one setting is a lone dot in a
  ## panel captioned like a curve. Skip rather than draw it.
  if (length(unique(TD[[spec$var]])) < 2) {
    message(sw, ": only one setting present -- skipped"); return(invisible(NULL))
  }

  agg <- do.call(rbind, lapply(split(TD, list(TD[[spec$var]], TD$method), drop = TRUE),
    function(z) data.frame(x = z[[spec$var]][1], method = z$method[1], n_rep = nrow(z),
                           secs = mean(z$secs, na.rm = TRUE),
                           secs_sd = sd(z$secs, na.rm = TRUE),
                           stringsAsFactors = FALSE)))
  agg <- agg[order(agg$x, agg$method), ]
  names(agg)[names(agg) == "x"] <- spec$var

  miss <- setdiff(neworder, unique(agg$method))
  if (length(miss))
    message("methods with no timing on ", sw, ": ", paste(miss, collapse = ", "))

  write.csv(agg, file.path(TAG_DIR, paste0(spec$stem, sfx, ".csv")), row.names = FALSE)
  attr(agg, "TD") <- TD
  agg
}

plot_timing <- function(sw, agg, spec, ybreaks) {
  TD <- attr(agg, "TD")

  ## LINEAR y, with evenly spaced breaks (changed from log10 on 2026-09-05 at the
  ## user's request). The trade-off is worth stating in the caption: the measured
  ## spread is four decades -- limma 0.02 s against ANCOM 1421 s at D = 1000 --
  ## so on an even scale every method under ~100 s sits on the axis and only the
  ## two or three most expensive are separable. That is the correct picture for
  ## "what does a replicate actually cost", and the wrong one for comparing cost
  ## SCALING between cheap methods, which is what a log axis showed as a slope.
  agg$method <- factor(agg$method, levels = neworder)
  ## Geometry and theme reproduce the previous manuscript figure exactly:
  ## theme_minimal() (theme_bw drew the black panel border), size-1 lines,
  ## size-2 points, axis.text 13, axis.title 21, legend text/title 12 at
  ## key.size 1.5 lines. No plot.title: the old script set "(a)"/"(b)" there,
  ## and those panel labels are added in Overleaf instead.
  ##
  ## Y LABELS: bare numbers, as in the old figure. The unit lives in the axis
  ## title ("Average Running Time (s)"), so repeating " s" on every break both
  ## departs from the old look and says the same thing twice.
  p <- ggplot(agg, aes(x = .data[[spec$var]], y = secs, group = method,
                       color = method, linetype = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(name = "Method", values = method_colors, drop = FALSE) +
    scale_linetype_manual(name = "Method", values = line_types, drop = FALSE) +
    ## Mantissa-and-exponent tick labels on the log x axis. The lambda breaks are
    ## 10,000 to 10,000,000, and spelled out in full they collide at this width;
    ## "1e4" ... "1e7" fits. Written from the break value rather than with
    ## format(scientific=TRUE), which yields "1e+05" / "5e+05" -- the plus sign
    ## and zero padding are exactly the width being saved here.
    { if (isTRUE(spec$xlog))
        scale_x_log10(breaks = sort(unique(agg[[spec$var]])),
                      labels = function(v) {
                        e <- floor(log10(v)); m <- round(v / 10^e, 3)
                        paste0(ifelse(m == 1, "1", format(m, trim = TRUE)), "e", e)
                      }) else NULL } +
    scale_y_continuous(breaks = ybreaks, limits = range(ybreaks)) +
    guides(color = guide_legend(override.aes = list(
             linetype = unname(line_types[neworder]))),
           linetype = "none") +
    labs(x = spec$lab, y = "Average Running Time (s)",
         color = "Method", linetype = "Method") +
    theme_minimal() +
    theme(legend.position = "right",
          axis.text  = element_text(size = 13),
          axis.title = element_text(size = 21),
          legend.key.size = unit(1.5, "lines"),
          legend.title = element_text(size = 12),
          legend.text  = element_text(size = 12))

  dir.create(RES_FIG, recursive = TRUE, showWarnings = FALSE)
  for (ext in c("pdf", "png")) {
    f <- file.path(RES_FIG, paste0(spec$stem, sfx, ".", ext))
    ## The old script saved these standalone at 6 x 4 in. That canvas CLIPS here:
    ## 13 methods at legend.key.size 1.5 lines make the legend taller than 4 in,
    ## so the "Method" legend title is pushed off the top and the 21pt y-axis
    ## title is cut at the left edge. (The old figure had 12 entries, and its
    ## side-by-side version was saved at 18 x 6, where the legend fitted.)
    ##
    ## Enlarged to 7.5 x 5 in -- the same 1.5 aspect ratio, and NOT a change to
    ## any type size: axis.text 13, axis.title 21 and the legend sizes are all
    ## still the old script's absolute pt values, so the figure looks the same,
    ## just uncropped. Revert to width = 6, height = 4 to match byte for byte.
    ggsave(f, p, width = 7.5, height = 5, dpi = if (ext == "png") 150 else 300, bg = "white")
    message("wrote ", f)
  }

  ## What the supplementary nine-variant grid adds, for the caption.
  v <- TD[TD$method == "Sparse_SSRV", ]
  if (nrow(v) && any(v$secs_variants > 0)) {
    vs <- do.call(rbind, lapply(split(v, v[[spec$var]]), function(z) {
      out <- data.frame(z[[spec$var]][1], mean(z$secs), mean(z$secs_variants))
      names(out) <- c(spec$var, "sparse_single", "variants_extra"); out
    }))
    cat("\n## ", sw, ": extra cost of the supplementary 3-bandwidth grid (s/replicate)\n\n", sep = "")
    print(vs, row.names = FALSE, digits = 4)
  }
  invisible(p)
}

## ---------------------------------------------------------------------------
## One y scale for all three figures
## ---------------------------------------------------------------------------
## Aggregate every sweep FIRST, then derive a single set of breaks from the
## global maximum, so the three figures are read against the same ruler and a
## bar height means the same thing in each.
##
## The cost, which the caption should state: the global max is ANCOM at D = 1000
## (1421 s), while the lambda sweep tops out at 173 s and the N sweep at 410 s.
## On a shared axis the lambda panel therefore occupies the bottom ~12% and the N
## panel the bottom ~29%, and within them the cheap methods are indistinguishable.
## That is the price of cross-sweep comparability and it is the trade the user
## asked for; per-sweep free scales are one argument away (pass ybreaks per sweep).
.aggs <- list()
for (.sw in names(TIMING_SWEEPS)) {
  a <- agg_timing(.sw, TIMING_SWEEPS[[.sw]])
  if (!is.null(a)) .aggs[[.sw]] <- a
}
if (!length(.aggs)) stop("no sweep produced timings")

.ymax    <- max(vapply(.aggs, function(a) max(a$secs, na.rm = TRUE), numeric(1)))
.ybreaks <- pretty(c(0, .ymax), n = 8)
message("shared y range: 0 to ", max(.ybreaks), " s (global max ",
        round(.ymax, 1), " s)")

for (.sw in names(.aggs)) plot_timing(.sw, .aggs[[.sw]], TIMING_SWEEPS[[.sw]], .ybreaks)
