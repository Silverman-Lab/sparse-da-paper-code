################################################################################
## aggregate_main.R -- Study A figures (goals 1 and 2).
##
##   results/figures/main_<metric>_<sweep>.pdf   13-method comparison, ONE Sparse
##                                               SSRV curve and it is the HYBRID
##   results/figures/decomp_<metric>_<sweep>.pdf hybrid vs laplace-only vs
##                                               bootstrap-only, three lines
##   results/figures/scale_choice_<sweep>.pdf    goal 1: P(bootstrap), margin,
##                                               and the superseded rule
##   results/main/<tag>/aggregated.RData         everything, for re-plotting
##
## FDR convention: a method that flags nothing has FDR = 0/0. It is reported as 0
## (config's FDR_WHEN_NO_DISCOVERIES), matching the original aggregators. The
## per-replicate metrics_summary.csv keeps the raw NaN, so the
## conditional-on-discovery version is recoverable without re-running. n_sig is
## plotted alongside every FDR panel, because "discovered nothing" and
## "discovered only true positives" both give 0 and are not the same result.
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
## classify() and f.score(), for the MIMIX block below. Sourced in the same order
## as shared.R. Not optional and not redundant with io_compat.R: io_compat defines
## an f.score INSIDE mrckf_metrics(), so it never reaches the top level, and
## classify() exists only here. Until MIMIX actually produced output this branch
## never ran, so the missing source could not surface.
source(file.path(OLD_SIM, "src", "common", "tidy_metrics.R"))
source(file.path(OLD_SIM, "src", "common", "io_compat.R"))

suppressWarnings(suppressMessages({ library(ggplot2) }))

args <- commandArgs(trailingOnly = TRUE)
tag  <- if (length(args) >= 1) args[1] else MAIN_TAG
RES_TAG_DIR <- file.path(RES_MAIN, tag)
if (!dir.exists(RES_TAG_DIR)) stop("no results at ", RES_TAG_DIR)

## ---------------------------------------------------------------------------
## Collect
## ---------------------------------------------------------------------------
dirs <- list.dirs(RES_TAG_DIR, recursive = FALSE)
## METRICS_FILE lets an alternative per-replicate metrics table be aggregated
## without touching the originals -- used by redecide_scale_rule.R, which writes
## metrics_summary_<rule>.csv after re-selecting the SSRV arm under a different
## scale rule. The arms are already fitted and saved, so changing the rule needs
## a re-decision and a re-aggregation, never a refit.
METRICS_FILE <- if (length(args) >= 2) args[2] else "metrics_summary.csv"
## Every output -- tables AND figures -- carries this suffix, so aggregating an
## alternative rule cannot silently overwrite the default rule's figures.
##
## TWO independent sources, concatenated. The SCENARIO comes first and is read off
## the TAG, not off the environment: results/figures/ is shared by every tag, so
## the supplementary 20/80 run's figures must be distinguishable from the main
## 80/20 run's whether or not SIM_PROP_POS happens to be set in the shell this
## aggregation is launched from. The second part is the alternative scale rule
## (redecide_scale_rule.R), which names a different metrics file.
sfx <- paste0(scenario_sfx_of(tag),
              if (METRICS_FILE == "metrics_summary.csv") "" else
                sub("^metrics_summary", "", sub("\\.csv$", "", METRICS_FILE)))
if (nzchar(sfx)) message("aggregating ", METRICS_FILE, " from tag ", tag,
                         " -> outputs suffixed '", sfx, "'")
dirs <- dirs[file.exists(file.path(dirs, METRICS_FILE))]
if (length(dirs) == 0) stop("no completed replicates in ", RES_TAG_DIR)
message(length(dirs), " completed replicates")

## Parameters come from replicate_info.csv, not from the directory name: pi is
## variable-width in the name ("0_1" vs "0_02"), which is exactly the kind of
## parsing the original had to do from both ends.
read_rep <- function(d) {
  info <- read.csv(file.path(d, "replicate_info.csv"), stringsAsFactors = FALSE)
  s    <- read.csv(file.path(d, METRICS_FILE), stringsAsFactors = FALSE)

  ## --- MIMIX, folded into the SAME classification rule as everything else ----
  ## The original computed MIMIX's TP/FP/TN/FN inline with its own convention.
  ## Here sig is "the 95% credible interval excludes 0", est is the posterior
  ## mean scaled to log2, and classify() does the rest -- so MIMIX sits on the
  ## same axis as the other twelve instead of a parallel one.
  ## Isolated in its own tryCatch: the caller wraps read_rep() as a whole, so
  ## anything thrown here would discard the replicate's OTHER TWELVE methods and
  ## report it as an unreadable directory. A MIMIX problem must cost the MiMIX
  ## row and nothing else -- which the "methods absent from this result set"
  ## message downstream then reports.
  s <- tryCatch({
  mo <- file.path(MIMIX_SRC, "SparseDOSSA-analysis-results", tag, basename(d))
  mb <- tryCatch(read_mimix_beta_summary(mo), error = function(e) NULL)
  taxa_file <- file.path(d, "taxa.csv")
  if (!is.null(mb) && file.exists(taxa_file)) {
    taxa <- readLines(taxa_file)
    taxa <- gsub('^"|"$', "", taxa)
    if (nrow(mb) == length(taxa)) {
      truth <- read.csv(file.path(d, "filtered_metadata_df.csv"), stringsAsFactors = FALSE)
      tl <- setNames(rep(0, length(taxa)), taxa)
      hit <- as.character(truth$feature_spiked)
      hit <- hit[hit %in% names(tl)]
      tl[hit] <- truth$effect_size[match(hit, as.character(truth$feature_spiked))]

      keep <- taxa != "other"
      est  <- mb[["mean"]][keep] * (1 / log(2))     # MiMIX beta is natural-log CLR
      sig  <- (mb[["2.5%"]][keep] * mb[["97.5%"]][keep]) > 0
      true_log2 <- unname(tl[taxa[keep]]) / log(2)

      cl <- classify(sig, est, true_log2)
      ntp <- sum(cl$tp); nfp <- sum(cl$fp); ntn <- sum(cl$tn); nfn <- sum(cl$fn)
      ppv <- ntp / (ntp + nfp); tdr <- ntp / (ntp + nfn)
      shift <- mean(est[true_log2 == 0], na.rm = TRUE)
      s <- rbind(s, data.frame(
        rep = info$seed, arm = "sample", method = "MiMIX",
        N = info$N, D_kept = info$D_kept, lambda = info$lambda, pi = info$pi,
        fdr = nfp / (ntp + nfp), tdr = tdr, typ1 = nfp / (ntn + nfp), ppv = ppv,
        f0.5 = f.score(0.5, ppv, tdr), f1 = f.score(1, ppv, tdr),
        n_sig = ntp + nfp, n_spiked = sum(true_log2 != 0),
        shift_hat = shift, centered_bias = NA_real_, centered_rmse = NA_real_,
        raw_bias = NA_real_, slope = NA_real_,
        theta_perp_true = info$theta_perp_true, scale_recovered = NA_real_,
        stringsAsFactors = FALSE)[, names(s), drop = FALSE])
    }
  }
  s
  }, error = function(e) {
    warning("MiMIX row skipped for ", basename(d), ": ", conditionMessage(e),
            call. = FALSE)
    s
  })

  s$N <- info$N; s$D <- info$D; s$lambda <- info$lambda; s$pi <- info$pi
  s$D_kept <- info$D_kept; s$seed <- info$seed
  s$theta_perp_true <- info$theta_perp_true

  sm <- tryCatch(read.csv(file.path(d, "scale_model.csv"), stringsAsFactors = FALSE),
                 error = function(e) NULL)

  ## theta^perp, one file per bandwidth: theta_perp.csv (nrd0, the main path)
  ## plus theta_perp_bw_<tag>.csv from the supplementary grid. Each holds three
  ## rows -- hybrid / laplace / bootstrap -- so one replicate contributes up to
  ## nine, which is exactly the variant grid.
  ##
  ## The TRUTH is joined here rather than inside sim_sparse_auto(): the fit is not
  ## told theta_perp_true and must not be, or the coverage it reports would not
  ## be the coverage of a real analysis.
  ## --- SUPPLEMENTARY: the nine Sparse SSRV variants, scored ------------------
  ## 3 arms x 3 bandwidths. Scored HERE rather than stashed as methods in
  ## run_all_methods(), so metrics_summary.csv and `agg` stay exactly what they
  ## were and no supplementary variant can reach a main figure by accident.
  ##
  ## Same classify() and f.score() the twelve methods go through -- the point of
  ## the table is comparability, so a parallel scoring rule would defeat it.
  vr <- tryCatch({
    truth_df <- read.csv(file.path(d, "filtered_metadata_df.csv"), stringsAsFactors = FALSE)
    arms <- c(hybrid = "sparse_auto", laplace = "sparse_L", bootstrap = "sparse_B")
    bws  <- c("nrd0", "nrd", "SJ-ste")
    do.call(rbind, lapply(bws, function(bw) {
      sfx_bw <- if (bw == "nrd0") "" else paste0("_bw_", bw)
      ## Which arm the rule chose AT THIS BANDWIDTH -- the superscript the
      ## real-data Task 3 table carries. A hybrid row must duplicate the arm row
      ## it names; that invariant is checked in the verification script.
      smf <- file.path(d, paste0("scale_model", sfx_bw, ".csv"))
      chosen <- if (file.exists(smf))
        tolower(read.csv(smf, stringsAsFactors = FALSE)$chosen_hybrid[1]) else NA_character_
      do.call(rbind, lapply(names(arms), function(a) {
        f <- file.path(d, paste0(arms[[a]], "_filterFALSE", sfx_bw, ".csv"))
        if (!file.exists(f)) return(NULL)
        tb <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
        tb <- tb[rownames(tb) != "other", , drop = FALSE]
        if (!nrow(tb)) return(NULL)
        tl  <- setNames(rep(0, nrow(tb)), rownames(tb))
        hit <- as.character(truth_df$feature_spiked); hit <- hit[hit %in% names(tl)]
        tl[hit] <- truth_df$effect_size[match(hit, as.character(truth_df$feature_spiked))]
        true_log2 <- unname(tl) / log(2)
        cl  <- classify(as.logical(tb$sig.q), tb$mean, true_log2)
        ntp <- sum(cl$tp); nfp <- sum(cl$fp); ntn <- sum(cl$tn); nfn <- sum(cl$fn)
        ppv <- ntp / (ntp + nfp); tdr <- ntp / (ntp + nfn)
        data.frame(bandwidth = bw, arm = a, chosen_at_bw = chosen,
                   fdr = nfp / (ntp + nfp), tdr = tdr, ppv = ppv,
                   typ1 = nfp / (ntn + nfp),
                   f0.5 = f.score(0.5, ppv, tdr), f1 = f.score(1, ppv, tdr),
                   n_sig = ntp + nfp, tp = ntp, fp = nfp,
                   n_spiked = sum(true_log2 != 0),
                   stringsAsFactors = FALSE)
      }))
    }))
  }, error = function(e) {
    warning("variant scoring failed for ", basename(d), ": ", conditionMessage(e),
            call. = FALSE); NULL
  })
  if (!is.null(vr) && nrow(vr))
    vr <- cbind(vr, info[rep(1, nrow(vr)), c("N","D","lambda","pi","seed","D_kept")])

  perp_files <- list.files(d, pattern = "^theta_perp(_bw_.*)?\\.csv$", full.names = TRUE)
  tp <- do.call(rbind, lapply(perp_files, function(f)
    tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)))
  if (!is.null(tp) && nrow(tp)) {
    tp$theta_perp_true <- info$theta_perp_true
    ## Same definition as theta_perp_check.R uses against the external
    ## measurements, so simulated and measured coverage are one quantity.
    tp$covered <- tp$lower <= tp$theta_perp_true & tp$theta_perp_true <= tp$upper
    tp$bias    <- tp$mean - tp$theta_perp_true
    tp <- cbind(tp, info[rep(1, nrow(tp)), c("N","D","lambda","pi","seed")])
  } else tp <- NULL

  list(summary = s,
       scale = if (is.null(sm)) NULL else cbind(sm, info[, c("N","D","lambda","pi","seed")]),
       perp  = tp,
       variants = vr)
}

parts <- lapply(dirs, function(d) tryCatch(read_rep(d), error = function(e) {
  warning("could not read ", basename(d), ": ", conditionMessage(e)); NULL
}))
parts <- Filter(Negate(is.null), parts)

S  <- do.call(rbind, lapply(parts, `[[`, "summary"))
SC <- do.call(rbind, Filter(Negate(is.null), lapply(parts, `[[`, "scale")))
TP <- do.call(rbind, Filter(Negate(is.null), lapply(parts, `[[`, "perp")))
VR <- do.call(rbind, Filter(Negate(is.null), lapply(parts, `[[`, "variants")))

## ---------------------------------------------------------------------------
## FDR convention (see header)
## ---------------------------------------------------------------------------
S$fdr_raw <- S$fdr
S$fdr[S$n_sig == 0] <- FDR_WHEN_NO_DISCOVERIES
S$ppv[S$n_sig == 0]  <- NA_real_        # precision is genuinely undefined
S$f0.5[is.na(S$f0.5) & S$n_sig == 0] <- 0
S$f1[is.na(S$f1)     & S$n_sig == 0] <- 0

## ---------------------------------------------------------------------------
## Re-attach each parameter set to every sweep it belongs to
## ---------------------------------------------------------------------------
expand_sweeps <- function(df) {
  do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    sw <- main_sweeps_of(df$N[i], df$D[i], df$lambda[i], df$pi[i])
    if (length(sw) == 0) return(NULL)
    z <- df[rep(i, length(sw)), , drop = FALSE]
    z$sweep <- sw
    z$lever_value <- vapply(sw, function(s)
      switch(s, overN = df$N[i], overD = df$D[i],
             overLambda = df$lambda[i], overPi = df$pi[i]), numeric(1))
    z
  }))
}
S  <- expand_sweeps(S)
SC <- expand_sweeps(SC)

METRICS <- c(fdr = "FDR", tdr = "TPR", f0.5 = "F0.5", n_sig = "number selected")

agg <- do.call(rbind, lapply(split(S, list(S$sweep, S$lever_value, S$method), drop = TRUE),
  function(z) {
    out <- data.frame(sweep = z$sweep[1], lever_value = z$lever_value[1],
                      method = z$method[1], n_rep = nrow(z), stringsAsFactors = FALSE)
    for (m in names(METRICS)) {
      out[[m]] <- mean(z[[m]], na.rm = TRUE)
      out[[paste0(m, "_se")]] <- sd(z[[m]], na.rm = TRUE) / sqrt(sum(!is.na(z[[m]])))
    }
    out$pct_zero_disc <- 100 * mean(z$n_sig == 0)
    out$D_kept <- mean(z$D_kept)
    out
  }))

## ---------------------------------------------------------------------------
## theta^perp coverage, and the nine-variant supplementary grid
## ---------------------------------------------------------------------------
## One row per (setting, bandwidth, arm), averaged over replicates.
##
## Coverage is reported WITH bias, RMSE and mean width, never alone: an interval
## can reach 0.95 by being wide and an interval can miss by being biased, and
## those are opposite problems with opposite fixes. For the Laplace arm the
## expected finding is under-coverage from BIAS -- the mode of the LFC density
## need not sit at the true load shift, and no variance formula addresses that
## (see the kde_mode_var() header in ../sparse_SSRV.R).
perp_agg <- NULL
if (!is.null(TP) && nrow(TP)) {
  TP$bandwidth[is.na(TP$bandwidth)] <- "nrd0"
  perp_agg <- do.call(rbind, lapply(
    split(TP, list(TP$N, TP$D, TP$lambda, TP$pi, TP$bandwidth, TP$arm), drop = TRUE),
    function(z) data.frame(
      N = z$N[1], D = z$D[1], lambda = z$lambda[1], pi = z$pi[1],
      bandwidth = z$bandwidth[1], arm = z$arm[1], n_rep = nrow(z),
      coverage      = mean(z$covered, na.rm = TRUE),
      coverage_se   = sd(z$covered, na.rm = TRUE) / sqrt(sum(!is.na(z$covered))),
      bias          = mean(z$bias, na.rm = TRUE),
      rmse          = sqrt(mean(z$bias^2, na.rm = TRUE)),
      mean_width    = mean(z$width, na.rm = TRUE),
      mean_est      = mean(z$mean, na.rm = TRUE),
      theta_perp_true = mean(z$theta_perp_true, na.rm = TRUE),
      ## how often the hybrid RULE landed on each arm at this bandwidth; NA on
      ## the arm rows themselves, where `selected` is the same value repeated
      p_bootstrap   = mean(z$selected == "bootstrap", na.rm = TRUE),
      stringsAsFactors = FALSE)))
  perp_agg <- perp_agg[order(perp_agg$N, perp_agg$D, perp_agg$lambda, perp_agg$pi,
                             perp_agg$bandwidth, perp_agg$arm), ]
  write.csv(perp_agg, file.path(RES_TAG_DIR, paste0("theta_perp_coverage", sfx, ".csv")),
            row.names = FALSE)
  message(sprintf("theta^perp coverage: %d cells over %d replicate-rows -> %s",
                  nrow(perp_agg), nrow(TP),
                  paste0("theta_perp_coverage", sfx, ".csv")))
} else {
  message("no theta_perp.csv found in any replicate -- coverage table skipped")
}

## ---------------------------------------------------------------------------
## The nine-variant supplementary table, averaged over replicates
## ---------------------------------------------------------------------------
## One row per (setting, bandwidth, arm). `p_chosen_at_bw` is the superscript the
## real-data Task 3 table carries as a letter: the fraction of replicates where
## the width rule picked THIS arm AT THIS BANDWIDTH. It is what makes the hybrid
## row readable -- a hybrid FDR between the two arm FDRs means the rule was
## switching between them, not that it found a third behaviour.
variant_agg <- NULL
if (!is.null(VR) && nrow(VR)) {
  VR$fdr[VR$n_sig == 0] <- FDR_WHEN_NO_DISCOVERIES
  VR$ppv[VR$n_sig == 0] <- NA_real_
  VR$f0.5[is.na(VR$f0.5) & VR$n_sig == 0] <- 0
  VR$f1[is.na(VR$f1)     & VR$n_sig == 0] <- 0
  variant_agg <- do.call(rbind, lapply(
    split(VR, list(VR$N, VR$D, VR$lambda, VR$pi, VR$bandwidth, VR$arm), drop = TRUE),
    function(z) data.frame(
      N = z$N[1], D = z$D[1], lambda = z$lambda[1], pi = z$pi[1],
      bandwidth = z$bandwidth[1], arm = z$arm[1], n_rep = nrow(z),
      fdr = mean(z$fdr, na.rm = TRUE), tdr = mean(z$tdr, na.rm = TRUE),
      f0.5 = mean(z$f0.5, na.rm = TRUE), n_sig = mean(z$n_sig, na.rm = TRUE),
      tp = mean(z$tp, na.rm = TRUE), fp = mean(z$fp, na.rm = TRUE),
      fdr_se = sd(z$fdr, na.rm = TRUE) / sqrt(sum(!is.na(z$fdr))),
      p_chosen_at_bw = mean(z$chosen_at_bw == z$arm[1], na.rm = TRUE),
      stringsAsFactors = FALSE)))
  variant_agg <- variant_agg[order(variant_agg$N, variant_agg$D, variant_agg$lambda,
                                   variant_agg$pi, variant_agg$bandwidth,
                                   variant_agg$arm), ]
  write.csv(variant_agg, file.path(RES_TAG_DIR, paste0("ssrv_variants", sfx, ".csv")),
            row.names = FALSE)
  message(sprintf("SSRV variant grid: %d cells over %d replicate-rows -> %s",
                  nrow(variant_agg), nrow(VR), paste0("ssrv_variants", sfx, ".csv")))
} else {
  message("no Sparse SSRV variant tables found -- supplementary grid skipped")
}

dir.create(RES_FIG, recursive = TRUE, showWarnings = FALSE)
save(S, SC, TP, VR, agg, perp_agg, variant_agg,
     file = file.path(RES_TAG_DIR, paste0("aggregated", sfx, ".RData")))
write.csv(agg, file.path(RES_TAG_DIR, paste0("aggregated_metrics", sfx, ".csv")),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## Which sweeps are far enough along to plot
## ---------------------------------------------------------------------------
## Study A is run one sweep at a time, so this executes against a partial result
## set. expand_sweeps() correctly attaches the baseline to ALL FOUR sweeps, which
## means that after the first sweep the other three each hold exactly one point.
## Plotting those produces three figures per metric showing a lone baseline dot,
## indistinguishable at a glance from a real curve that happens to be flat. The
## tables above keep every setting either way.
n_lev <- tapply(agg$lever_value, agg$sweep, function(v) length(unique(v)))
plot_sweeps <- names(n_lev)[n_lev >= 2]
if (length(plot_sweeps) < length(n_lev))
  message("sweeps with <2 settings so far, not plotted: ",
          paste(setdiff(names(n_lev), plot_sweeps), collapse = ", "))

## ---------------------------------------------------------------------------
## Figures
## ---------------------------------------------------------------------------
LOGX <- c(overN = TRUE, overD = TRUE, overLambda = TRUE, overPi = FALSE)

## MAIN figure: one Sparse SSRV curve, and it is the HYBRID rule. The L / B arms
## are deliberately absent -- hybrid is the proposed method, so it is what the
## comparison against competitors must show. The decomposition is its own figure.
MAIN_METHODS <- c("ALDEx2_def", "ALDEx2_clr", "DESeq2", "Sparse_SSRV", "ANCOM",
                  "ANCOM-BC", "ANCOM-BC2", "MiMIX", "LinDA", "edgeR", "limma",
                  "baySeq", "MRCKF")
DECOMP_METHODS <- c("Sparse_SSRV", "Sparse_SSRV_L", "Sparse_SSRV_B")
DECOMP_LAB <- c(Sparse_SSRV = "hybrid (selected)",
                Sparse_SSRV_L = "Laplace only",
                Sparse_SSRV_B = "bootstrap only")

panel <- function(d, m, title, lab) {
  ggplot(d, aes(lever_value, .data[[m]], colour = method, group = method)) +
    ## Nominal 0.05 level. FDR ONLY -- 0.05 is the rate the BH step is targeting,
    ## so it is the reference the FDR panel is read against; it means nothing on
    ## the TPR / F0.5 / n_sig panels. Added FIRST so it renders BEHIND the curves
    ## rather than cutting across them.
    { if (identical(m, "fdr"))
        geom_hline(yintercept = 0.05, colour = "grey55",
                   linetype = "dashed", linewidth = 0.4) else NULL } +
    geom_line() + geom_point(size = 1.2) +
    geom_errorbar(aes(ymin = .data[[m]] - .data[[paste0(m, "_se")]],
                      ymax = .data[[m]] + .data[[paste0(m, "_se")]]),
                  width = 0, alpha = 0.35) +
    { if (isTRUE(LOGX[[d$sweep[1]]])) scale_x_log10() else NULL } +
    labs(x = sweep_label(d$sweep[1]), y = METRICS[[m]], colour = NULL,
         title = paste0(title, " -- ", METRICS[[m]], " over ", d$sweep[1]),
         subtitle = lab) +
    theme_bw()
}

## Phase 0c runs SSRV only, so most of MAIN_METHODS is absent. Drop unused factor
## levels rather than carrying empty ones, which would otherwise produce a legend
## of 13 entries for a 1-method plot.
present_main <- intersect(MAIN_METHODS, unique(agg$method))
if (length(present_main) < length(MAIN_METHODS)) {
  message("methods absent from this result set (expected for an SSRV-only run): ",
          paste(setdiff(MAIN_METHODS, present_main), collapse = ", "))
}

for (sw in plot_sweeps) {
  for (m in names(METRICS)) {
    d <- agg[agg$sweep == sw & agg$method %in% present_main, ]
    if (nrow(d) == 0) next
    d$method <- factor(d$method, levels = present_main)
    ggsave(file.path(RES_FIG, sprintf("main_%s_%s%s.pdf", m, sw, sfx)),
           panel(d, m, "all methods",
                 "Sparse SSRV = hybrid rule; FDR = 0 when nothing is flagged (see n_sig panel)"),
           width = 7.5, height = 5)

    e <- agg[agg$sweep == sw & agg$method %in% DECOMP_METHODS, ]
    if (nrow(e) == 0) next
    e$method <- factor(DECOMP_LAB[e$method], levels = unname(DECOMP_LAB))
    ggsave(file.path(RES_FIG, sprintf("decomp_%s_%s%s.pdf", m, sw, sfx)),
           panel(e, m, "Sparse SSRV scale model",
                 "all three from ONE fit -- the hybrid table equals whichever arm was selected"),
           width = 7, height = 4.5)
  }
}

## ---------------------------------------------------------------------------
## Goal 1: which scale model, per setting
## ---------------------------------------------------------------------------
if (!is.null(SC) && nrow(SC) > 0) {
  wilson <- function(k, n, conf = 0.95) {
    if (n == 0) return(c(NA_real_, NA_real_))
    z <- qnorm(1 - (1 - conf) / 2); ph <- k / n; den <- 1 + z^2 / n
    ctr <- (ph + z^2 / (2 * n)) / den
    hw <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / den
    c(max(0, ctr - hw), min(1, ctr + hw))
  }
  sc <- do.call(rbind, lapply(split(SC, list(SC$sweep, SC$lever_value), drop = TRUE),
    function(z) {
      k <- sum(z$chosen == "Bootstrap"); n <- nrow(z); ci <- wilson(k, n)
      data.frame(sweep = z$sweep[1], lever_value = z$lever_value[1], n_rep = n,
                 p_bootstrap = k / n, lo = ci[1], hi = ci[2],
                 p_bootstrap_old = mean(z$chosen_old == "Bootstrap"),
                 margin = mean(z$margin, na.rm = TRUE),
                 laplace_fail = mean(z$n_na_draws),
                 resample_reject = if ("n_resample_reject" %in% names(z))
                   mean(z$n_resample_reject) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  write.csv(sc, file.path(RES_TAG_DIR, paste0("scale_choice_summary", sfx, ".csv")),
            row.names = FALSE)

  for (sw in intersect(plot_sweeps, unique(sc$sweep))) {
    d <- sc[sc$sweep == sw, ]
    p1 <- ggplot(d, aes(lever_value, p_bootstrap)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
      geom_line() + geom_point() +
      geom_line(aes(y = p_bootstrap_old), linetype = 2, colour = "firebrick") +
      geom_point(aes(y = p_bootstrap_old), colour = "firebrick", shape = 2) +
      { if (isTRUE(LOGX[[sw]])) scale_x_log10() else NULL } +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = sweep_label(sw), y = "P(bootstrap selected)",
           title = paste("scale model choice --", sw),
           subtitle = "solid = current rule (Wilson 95% CI); dashed red = superseded mean-sd rule, same fits") +
      theme_bw()
    p2 <- ggplot(d, aes(lever_value, margin)) +
      geom_hline(yintercept = 0, linetype = 2) + geom_line() + geom_point() +
      { if (isTRUE(LOGX[[sw]])) scale_x_log10() else NULL } +
      coord_cartesian(ylim = c(-1, 1)) +
      labs(x = sweep_label(sw), y = "mean margin (statB - statL)",
           subtitle = "-1 = Laplace wins for every taxon; 0 = the two models are indistinguishable") +
      theme_bw()
    pdf(file.path(RES_FIG, sprintf("scale_choice_%s%s.pdf", sw, sfx)), width = 7, height = 4.5)
    print(p1); print(p2); dev.off()
  }
}

message("figures -> ", RES_FIG)
message("tables  -> ", RES_TAG_DIR)
