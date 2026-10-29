# Backward-compatible readers for the two files foo.R now deletes after the MIMIX fit.
#
# Older result directories (anything generated before the disk-reduction change) still
# contain the full Y.csv and beta.tsv. Newer ones carry taxa.csv and beta_summary.csv
# instead. Both readers prefer the compact form and fall back to the original, so the
# aggregators work unchanged against old and new output alike.

# Stand-in for `Y <- read.csv("Y.csv")`.
#
# The aggregators use the result for exactly two things: an existence guard, and
# colnames() to match spiked feature names. So a zero-row frame carrying the taxa
# names as columns is a faithful drop-in. Returns NULL when neither file is present,
# which the callers' `if (!is.null(Y))` guard already handles.
read_Y_compat <- function(dir = ".") {
  taxa_path <- file.path(dir, "taxa.csv")
  if (file.exists(taxa_path)) {
    nm <- readLines(taxa_path, warn = FALSE)
    nm <- nm[nzchar(nm)]
    df <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = length(nm)))
    # read.csv applies make.names by default; match that so colnames() is identical
    # to what the old `read.csv("Y.csv")` path produced.
    names(df) <- make.names(nm, unique = TRUE)
    return(df)
  }
  y_path <- file.path(dir, "Y.csv")
  if (file.exists(y_path)) return(read.csv(y_path))
  NULL
}

# MRCKF metrics for one replicate, from its per-taxon MRCKF.out file.
#
# MRCKF does not follow the generic TP/FP rule the other methods use, for two reasons
# that are properties of the method, not choices:
#   * its statistic is a knockoff W (an importance score with no direction), so there is
#     no sign to check -- exactly as ANCOM's detected_0.7 is handled;
#   * recall is over the spiked taxa that survived the filter, and the true-negative
#     count includes the aggregated "other" row.
# These are sparse_MRCKF.R:121-134's definitions verbatim, so the numbers stay comparable
# to the CKF curves produced before CKF was merged into the main sweep.
mrckf_metrics <- function(mrckf.out, n_spiked_kept) {
  sel <- mrckf.out$selected
  spiked <- mrckf.out$effect_size != 0
  tp <- sum(sel & spiked); fp <- sum(sel & !spiked)
  tn <- sum(!sel & !spiked); num_selected <- tp + fp
  tdr <- tp / n_spiked_kept
  ppv <- tp / (tp + fp)
  f.score <- function(f, precision, recall) {
    (1 + f^2) * (precision * recall) / (f^2 * precision + recall)
  }
  c(num_selected = num_selected,
    FDP.S        = fp / num_selected,
    # FDP.Sp duplicated FDP.S in the original too: MRCKF() returns
    # predictorNames[result$S] for BOTH S and Sp, so the knockoff+ result was discarded.
    FDP.Sp       = fp / num_selected,
    mrckf_tdr    = tdr,
    mrckf_ppv    = ppv,
    mrckf_f1     = f.score(1,   ppv, tdr),
    mrckf_f0.5   = f.score(0.5, ppv, tdr),
    mrckf_typ1   = fp / (tn + fp),
    mrckf_fdr    = fp / num_selected)
}

# Average MRCKF across replicates for each sweep value, replacing the readRDS of the
# {D,N,lambda,pi}_out.RData that the separate sparse_MRCKF.R sweep used to produce.
#
# Reads the same per-taxon MRCKF.out file the worker writes for every other method and
# derives the metrics here, so the definitions live in ONE place rather than being
# duplicated across the four Results scripts -- which is how the differing TP/FP rules
# between ALDEx2 and ANCOM-BC crept into the existing aggregators.
#
# res_dir   the res_foo/<tag> directory holding the replicate folders
# vals      sweep values, in plot order
# dir_for   function(val, rep) -> the replicate directory NAME for that value
# n, d      used to build the "MRCKF.out N<n>D<d>.csv" filename
# reps      replicate indices to average over
#
# Sweep values with no readable replicate are dropped rather than returned as NA rows, so
# a partial run yields a shorter curve instead of a broken one.
collect_mrckf <- function(res_dir, vals, dir_for, n, d, reps = 1:100) {
  rows <- lapply(vals, function(v) {
    per <- lapply(reps, function(r) {
      dir_r <- file.path(res_dir, dir_for(v, r))
      # d varies along the D sweep, so glob rather than hardcode the suffix
      f <- Sys.glob(file.path(dir_r, "MRCKF.out N*.csv"))
      if (!length(f)) return(NULL)
      mo <- tryCatch(read.csv(f[1], row.names = 1), error = function(e) NULL)
      fm <- tryCatch(read.csv(file.path(dir_r, "filtered_metadata_df.csv")),
                     error = function(e) NULL)
      if (is.null(mo) || is.null(fm) || nrow(mo) == 0) return(NULL)
      m <- mrckf_metrics(mo, nrow(fm))
      # Timing is not derivable from the per-taxon file -- it now lives with every other
      # method's in execution time.csv, which is where `running time comparison.R` needs
      # to find it.
      et <- tryCatch(read.csv(file.path(dir_r, "execution time.csv"), row.names = 1),
                     error = function(e) NULL)
      m["mrckf.time"] <- if (!is.null(et) && "mrckf.time" %in% rownames(et) &&
                             "elapsed" %in% colnames(et))
                           et["mrckf.time", "elapsed"] else NA_real_
      m
    })
    per <- Filter(Negate(is.null), per)
    if (!length(per)) return(NULL)
    colMeans(do.call(rbind, per), na.rm = TRUE)
  })
  ok <- !vapply(rows, is.null, logical(1))
  if (!any(ok)) {
    warning("collect_mrckf: no MRCKF.out files found under ", res_dir, call. = FALSE)
    return(NULL)
  }
  cbind(vals = vals[ok], as.data.frame(do.call(rbind, rows[ok])))
}

# Stand-in for `readr::read_tsv(".../beta.tsv")` followed by a per-column quantile.
#
# Returns a data.frame with one row per taxon and columns "2.5%", "50%", "97.5%",
# "mean" -- already summarized. Callers should index the quantile columns directly
# rather than re-applying quantile(). Returns NULL if neither file exists.
read_mimix_beta_summary <- function(dir) {
  summary_path <- file.path(dir, "beta_summary.csv")
  if (file.exists(summary_path)) {
    return(read.csv(summary_path, row.names = 1, check.names = FALSE))
  }
  beta_path <- file.path(dir, "beta.tsv")
  if (file.exists(beta_path)) {
    beta <- as.matrix(read.delim(beta_path, check.names = FALSE))
    s <- t(apply(beta, 2, function(v) {
      c(quantile(v, c(0.025, 0.5, 0.975), na.rm = TRUE), mean = mean(v, na.rm = TRUE))
    }))
    return(as.data.frame(s, check.names = FALSE))
  }
  NULL
}
