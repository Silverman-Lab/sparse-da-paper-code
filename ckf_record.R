################################################################################
## ckf_record.R -- CKF (the MRCKF knockoff filter) run INSIDE each real-data DA
## script, instead of in a separate pass whose results are merged later.
##
## Why this exists. CKF used to live in six standalone scripts under
## ~/Research/sparse_DA_Roar, each of which re-parsed its study from scratch,
## rebuilt the phyloseq object, re-applied the taxon filter and re-derived the
## truth set. Six independent copies of the data prep is six chances to drift
## from the analysis script, and they did drift:
##
##   * prasse_mrckf.R scored CKF against a tax_scale from Mar 29 while the
##     analysis script had rewritten it -- 6 gold taxa vs 12.
##   * selex_mrckf.R computed FP with the same expression as TP (a missing "!"),
##     so every selex CKF cell in the summary table was the TP count printed
##     twice: 0.74/(0.74).
##   * vandeputte_mrckf.R and zemb_mrckf.R saved only into ~/Research/datarepo,
##     so the project tree looked like it had no CKF result at all.
##   * stammler_mrckf.R's mean>5 filter was applied here but commented out in
##     the analysis script, so CKF was scored on a different taxon set from
##     every other method.
##
## None of those are possible once CKF is handed the SAME `Y_raw` and `X` that
## sparse_auto() is given, a few lines later in the same script. That is this
## file's entire interface:
##
##   ckf <- run_ckf(Y_raw, X, truth = zemb_truth, dataset = "Zemb2020")
##   res_tpfp <- add_ckf_row(res_tpfp, ckf)          # CKF row, already merged
##   record_ckf(ckf, dir = DIR_ZEMB, file = "zemb_mrckf_out.RData")
##
## `add_ckf_row()` is the line that removes the merge from `real data results.R`:
## res_tpfp$aldex_gold comes off disk with CKF already in it, scored against the
## same truth as every other method in the table.
##
## Backward compatibility. record_ckf() still writes the legacy per-study file
## (zemb_mrckf_out.RData, jin_mrckf_out.RData, ...) in the legacy shape
## -- list(tpfp, mean, median, ratio_mean) -- so the existing
## `ckf_zemb$mean` / `colMeans(ckf_stammler[[5]])` lines keep working unchanged
## while the reporting script is being simplified.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("MRCKF", quietly = TRUE))
    stop("ckf_record.R needs the MRCKF package")
  if (!requireNamespace("compositions", quietly = TRUE))
    stop("ckf_record.R needs the compositions package")
})

CKF_RECORD_FILE <- "ckf_record.rds"
CKF_RECORD_CSV  <- "ckf_record.csv"

################################################################################
## ckf_W: the design matrix MRCKF wants, built from the analysis matrix.
##
## MRCKF() takes W as SAMPLES x TAXA with a pseudocount, while every DA script
## carries counts as `Y_raw`, TAXA x SAMPLES. The six Roar scripts each did this
## transpose by hand off a phyloseq object, and each got there differently
## depending on that study's `taxa_are_rows` -- t(otu_table())+1 in one,
## otu_table()+1 in the next. Deriving it from Y_raw instead makes the
## orientation a property of this file rather than of six separate scripts.
##
## THE LAST COLUMN OF W IS THE LOG-RATIO REFERENCE AND CAN NEVER BE SELECTED.
## MRCKF() sets p = ncol(W) - 1 and builds X[,t] = log(Z[,t]/W[,pstar]), so the
## final taxon has no design column at all. On a study whose truth set is a
## single taxon -- Zemb's Seq1, Jin's COTU-CM-2002 -- that taxon landing last
## makes CKF's true-positive count structurally zero, which would read as a
## method failure rather than as an artifact of column order. run_ckf() checks
## for exactly this and refuses to score silently.
################################################################################
ckf_W <- function(Y_raw, reference = NULL) {
  W <- t(as.matrix(Y_raw)) + 1                 # samples x taxa, pseudocount
  storage.mode(W) <- "double"
  if (is.null(colnames(W)))
    stop("Y_raw has no rownames -- CKF results could not be scored against any ",
         "truth set, because the selected taxa would have no names")

  ## Move a named taxon into the reference slot. Supplying one is the way to
  ## keep a truth taxon out of the position where it cannot be selected.
  if (!is.null(reference)) {
    j <- match(reference, colnames(W))
    if (is.na(j)) stop("reference taxon '", reference, "' is not in Y_raw")
    W <- W[, c(setdiff(seq_len(ncol(W)), j), j), drop = FALSE]
  }
  W
}

## X as MRCKF wants it: a one-column matrix of 0/1, in the same sample order as
## the columns of Y_raw.
##
## The DA scripts spell the contrast four different ways -- numeric 0/1 in selex
## and Vandeputte, logical in Zemb, a character Treatment in Stammler, a
## character TYPE in Prasse. Guessing which level means "treatment" is how a
## table ends up with TP and FP silently swapped, so a two-level character or
## factor has to name its treatment level explicitly.
ckf_X <- function(X, treatment = NULL) {
  if (is.logical(X)) return(as.matrix(as.numeric(X)))
  if (is.numeric(X)) {
    if (!all(X %in% c(0, 1)))
      stop("numeric X must be coded 0/1; got values ",
           paste(utils::head(sort(unique(X)), 5), collapse = ", "))
    return(as.matrix(as.numeric(X)))
  }
  lv <- sort(unique(as.character(X)))
  if (length(lv) != 2)
    stop("X has ", length(lv), " levels (", paste(lv, collapse = ", "),
         ") -- CKF handles a two-group contrast only")
  if (is.null(treatment))
    stop("X is a ", class(X)[1], " with levels '", lv[1], "' and '", lv[2],
         "' -- pass treatment = one of them so the contrast direction is not guessed")
  if (!treatment %in% lv)
    stop("treatment = '", treatment, "' is not a level of X (", paste(lv, collapse = ", "), ")")
  as.matrix(as.numeric(as.character(X) == treatment))
}

################################################################################
## ckf_once: one knockoff filter fit, returning BOTH selection sets.
##
## This reproduces MRCKF()'s single-response path, but only the deterministic
## part of it -- the acomp() closure and the log-ratio transform, six lines with
## no RNG -- and then calls the package's own runSingleKF() for everything
## stochastic. The knockoff draw and both knockoff.filter() calls therefore
## happen inside MRCKF exactly as they would if MRCKF() were called directly, in
## the same order, consuming the same random numbers. `S` is identical to what
## MRCKF() returns.
##
## Two defects in the installed package are reproduced rather than fixed, so
## these numbers stay comparable with every CKF figure already computed:
##
##   1. `ref` is bound to the RAW last column of W while the correctly closed
##      reference is computed into `res` and never used. So
##      X[i,t] = alr(i,t) - log(rowSum_i): every design column of row i carries
##      the same -log(n_i) offset, which the later column-wise scaling cannot
##      remove. MRCKF as installed is NOT scale-invariant. Noted identically in
##      simulation/src/common/mrckf_fit.R.
##
##   2. MRCKF() returns list(S = names[result$S], Sp = names[result$S]) -- `Sp`
##      is assigned from `result$S`, not `result$Sp`. In the single-response
##      setting every study here uses, MRCKF()'s knockoff+ set is therefore a
##      copy of its knockoff set, which is why the `tp2`/`fp2` vectors in all six
##      Roar scripts were duplicates of `tp`/`fp` and were never reported.
##      Calling runSingleKF() directly returns the real offset = 1 set, which is
##      what `offset = 1` below selects.
################################################################################
ckf_once <- function(W, x, q = 0.05, penalty = "lasso") {
  pstar <- ncol(W); p <- pstar - 1L
  ref <- W[, pstar]                                  # raw, per the package
  Z   <- as.matrix(compositions::acomp(W))
  Xd  <- matrix(0, nrow = nrow(W), ncol = p)
  for (t in seq_len(p)) Xd[, t] <- log(Z[, t] / ref)

  ## runSingleKF() prints a banner per call; 100 replicates of it buries the run
  ## log. The result is the return value, not the printing.
  utils::capture.output(
    res <- MRCKF::runSingleKF(Xd, as.vector(x), q, penalty),
    type = "output")

  nm <- colnames(W)
  list(S  = nm[res$S],      # offset = 0, the set MRCKF() reports
       Sp = nm[res$Sp])     # offset = 1, knockoff+ -- MRCKF() loses this one
}

################################################################################
## run_ckf: R replicates of the filter, scored against a truth set.
##
## WHY REPLICATES AT ALL. The knockoff filter is randomised twice over on every
## call -- create.second_order() draws a fresh knockoff matrix, and
## stat.glmnet_coefdiff() draws fresh cross-validation folds -- so one fit is one
## draw from a selection distribution, not an answer. With n in the 10-95 range
## and p in the hundreds, that distribution is wide, which is why these runs were
## always averaged over 100 fits.
##
## WHAT DRIVES THE WIDTH. knockoff.threshold() at offset = 0 takes the smallest
## t with
##
##      #{W <= -t} / max(1, #{W >= t})  <=  q
##
## and the textbook reading of that -- "at least 1/q taxa must clear the
## threshold, so raise q" -- turns out NOT to be what binds here. Measured on
## Zemb (D = 198, N = 9, 15 replicates, seed 1234):
##
##   * a replicate has only ~4.5 positive W statistics on average, and just
##     1 of 15 had none at all;
##   * 6 of 15 returned the empty set anyway, because some negative W was as
##     extreme as every positive one, putting the ratio at 1/k with k <= 5;
##   * and q = 0.20 selected EXACTLY the same set as q = 0.05 in all 15.
##
## With only a handful of positive statistics the ratio can essentially only be
## 0/k or >= 1/k, so every q between 0 and ~0.2 lands in the same gap and behaves
## identically. Raising q is not a usable lever at this n -- it would have to go
## past ~1/(#positive W), around 0.25, and by then the guarantee is being applied
## to five candidates and means very little. Whether a given replicate selects
## anything is a property of the knockoff DRAW, not of q.
##
## So the empty/non-empty split is real and is not tunable away per-replicate.
## What does work is pooling ACROSS draws:
##
##   derand_thresh=
##                derandomised knockoffs: keep the taxa selected in at least this
##                fraction of the replicates and report that ONE set, instead of
##                averaging R noisy TP/FP pairs. A taxon that clears the bar in
##                some draws accumulates frequency even when other draws return
##                nothing, which is exactly the failure mode above. The
##                replicates stop being a precision device and become the
##                estimator. Costs nothing extra -- same fits. `derand_tpfp`.
##   se_target=   stop once the Monte Carlo SE of both mean TP and mean FP is
##                below a tolerance (and at least R_min fits have run). The SE
##                falls as 1/sqrt(R), so a study that has settled by R = 30 pays
##                for 30 rather than 100. `R_used` records what it took.
##   q=           worth CHECKING per study rather than assuming, since the above
##                was measured on Zemb only. `empty_rate` at two values of q
##                answers it in one extra run; if it does not move, q is not the
##                lever for that study either.
##
## `offset = 1` (knockoff+) is orthogonal and is the conservative direction, not
## the stabilising one: on the same Zemb data it took the empty rate from 45% to
## 95%. Use it when the concern is false positives, not variance.
##
## Returns the legacy shape -- tpfp / mean / median / ratio_mean -- so anything
## already reading these files keeps working, plus the fields above.
################################################################################
run_ckf <- function(Y_raw, X, truth,
                    R = 100, q = 0.05, penalty = "lasso", offset = 0,
                    treatment = NULL, reference = NULL, seed = 1234,
                    derand_thresh = 0.5, se_target = NULL, R_min = 20,
                    dataset = NA_character_, pair = NA_character_,
                    verbose = TRUE) {

  stopifnot(offset %in% c(0, 1), R >= 1, q > 0, q < 1)
  W <- ckf_W(Y_raw, reference = reference)
  x <- ckf_X(X, treatment = treatment)
  if (nrow(W) != length(x))
    stop("Y_raw has ", ncol(Y_raw), " samples but X has ", length(x),
         " -- these must be the same samples, in the same order")

  truth <- as.character(truth)
  taxa  <- colnames(W)
  ref_taxon <- taxa[ncol(W)]

  ## The reference taxon has no design column, so it cannot be selected. If a
  ## true positive is sitting there, CKF's TP is zero by construction and the
  ## resulting table row would be a statement about column order, not method
  ## performance.
  ref_is_truth <- ref_taxon %in% truth
  if (ref_is_truth)
    warning("CKF reference taxon '", ref_taxon, "' is in the truth set for ",
            dataset, " -- MRCKF drops the last column of W, so this taxon can ",
            "NEVER be selected and TP is structurally 0. Re-run with ",
            "reference = <some other taxon> to move it out of the reference slot.",
            call. = FALSE, immediate. = TRUE)

  n_truth_present <- sum(truth %in% taxa)
  if (n_truth_present == 0)
    warning("none of the ", length(truth), " truth taxa for ", dataset,
            " survive into W -- every TP will be 0. Check that `truth` uses the ",
            "same ids as rownames(Y_raw).", call. = FALSE, immediate. = TRUE)

  lab <- paste0(dataset, if (is.na(pair)) "" else paste0(" / ", pair))
  if (verbose) {
    cat(sprintf("\n== CKF (%s) : %d taxa x %d samples | R=%d q=%.3g offset=%d | truth %d/%d present\n",
                lab, ncol(W), nrow(W), R, q, offset, n_truth_present, length(truth)))
    ## Printed on EVERY run, not only when it collides with the truth set: which
    ## taxon MRCKF dropped is a fact about the fit that has to be checkable
    ## against the final result, and it is not recoverable from the output
    ## afterwards. Also stored as `reference_taxon` in ckf_record.rds / .csv --
    ## see ckf_reference_check() for the one-call view across every study.
    cat(sprintf("   reference taxon (dropped, never selectable): %s%s\n",
                ref_taxon, if (ref_is_truth) "   *** IN THE TRUTH SET ***" else ""))
  }

  set.seed(seed)
  selected <- vector("list", R)
  tp <- fp <- rep(NA_integer_, R)
  t0 <- Sys.time()
  R_used <- R

  for (i in seq_len(R)) {
    fit <- ckf_once(W, x, q = q, penalty = penalty)
    S   <- if (offset == 1) fit$Sp else fit$S
    selected[[i]] <- S
    tp[i] <- sum(S %in% truth)
    fp[i] <- sum(!S %in% truth)

    ## Monte Carlo SE of the two means. sd()/sqrt(i) is the right stopping
    ## quantity here even though the underlying distribution is bimodal -- it is
    ## the precision of the reported mean, which is what R buys.
    if (!is.null(se_target) && i >= max(R_min, 2)) {
      se_tp <- stats::sd(tp[1:i]) / sqrt(i)
      se_fp <- stats::sd(fp[1:i]) / sqrt(i)
      if (is.finite(se_tp) && is.finite(se_fp) &&
          se_tp <= se_target && se_fp <= se_target) {
        R_used <- i
        if (verbose)
          cat(sprintf("   stopped at R=%d: SE(TP)=%.3f SE(FP)=%.3f both <= %.3f\n",
                      i, se_tp, se_fp, se_target))
        break
      }
    }
    if (verbose && (i %% 10 == 0))
      cat(sprintf("   %3d/%d  TP %.2f  FP %.2f  empty %d\n", i, R,
                  mean(tp[1:i]), mean(fp[1:i]),
                  sum(lengths(selected[1:i]) == 0)))
  }

  selected <- selected[seq_len(R_used)]
  tp <- tp[seq_len(R_used)]; fp <- fp[seq_len(R_used)]
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  tpfp <- cbind(TP = tp, FP = fp)

  ## Per-taxon selection frequency over the replicates, and the derandomised
  ## set. Every taxon ever selected gets an entry; taxa never selected are left
  ## out rather than stored as a wall of zeros.
  hits <- table(unlist(selected))
  freq <- if (length(hits)) sort(setNames(as.numeric(hits) / R_used, names(hits)),
                                 decreasing = TRUE) else numeric(0)
  derand <- names(freq)[freq >= derand_thresh]
  derand_tpfp <- c(TP = sum(derand %in% truth), FP = sum(!derand %in% truth))

  n_empty <- sum(lengths(selected) == 0)

  out <- list(
    ## ---- legacy shape: what the six Roar scripts wrote, unchanged -----------
    tpfp       = tpfp,
    mean       = colMeans(tpfp),
    median     = c(stats::median(tpfp[, 1]), stats::median(tpfp[, 2])),
    ratio_mean = mean(apply(tpfp, 1, function(z)
                       if (sum(z) == 0) 0 else if (z[2] == 0) 1 else z[1] / z[2])),
    ## ---- what the counts alone could not say --------------------------------
    selected     = selected,          # keeps re-scoring free; see score_ckf()
    freq         = freq,
    derand       = derand,
    derand_tpfp  = derand_tpfp,
    n_empty      = n_empty,
    empty_rate   = n_empty / R_used,
    se           = c(TP = stats::sd(tp) / sqrt(R_used),
                     FP = stats::sd(fp) / sqrt(R_used)),
    sd           = c(TP = stats::sd(tp), FP = stats::sd(fp)),
    ## ---- provenance ---------------------------------------------------------
    dataset = dataset, pair = as.character(pair),
    R = R, R_used = R_used, q = q, penalty = penalty, offset = offset,
    seed = seed, derand_thresh = derand_thresh,
    D = ncol(W), N = nrow(W),
    n_truth = length(truth), n_truth_present = n_truth_present,
    reference_taxon = ref_taxon, reference_is_truth = ref_is_truth,
    elapsed_sec = elapsed,
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

  if (verbose) {
    cat(sprintf("   mean  TP %.2f (sd %.2f, se %.3f)   FP %.2f (sd %.2f, se %.3f)\n",
                out$mean[1], out$sd[1], out$se[1], out$mean[2], out$sd[2], out$se[2]))
    cat(sprintf("   empty selections: %d/%d (%.0f%%)%s\n", n_empty, R_used,
                100 * out$empty_rate,
                if (out$empty_rate > 0.25)
                  "  <- the spread is the empty/large split, not noise; try a larger q"
                else ""))
    cat(sprintf("   derandomised (freq >= %.2f): %d taxa, TP %d FP %d\n",
                derand_thresh, length(derand), derand_tpfp[1], derand_tpfp[2]))
    cat(sprintf("   %.1f s for %d fits (%.1f s/fit)\n", elapsed, R_used, elapsed / R_used))
  }
  out
}

################################################################################
## score_ckf: re-score stored selections against a different truth set.
##
## Prasse needs this. Its gold standard is not one set but a list over the scale
## sensitivity parameter gamma, and the old prasse_mrckf.R handled that by
## collecting the 100 selection sets once and then looping over tax_scale. That
## is the right shape -- the filter does not depend on the truth, so refitting
## per gamma would be 100 wasted fits per level -- and it is what this does.
##
## Returns an object in the same shape as run_ckf(), so add_ckf_row() and
## record_ckf() take it unchanged.
################################################################################
score_ckf <- function(ckf, truth, pair = ckf$pair) {
  truth <- as.character(truth)
  tp <- vapply(ckf$selected, function(S) sum(S %in% truth),  integer(1))
  fp <- vapply(ckf$selected, function(S) sum(!S %in% truth), integer(1))
  tpfp <- cbind(TP = tp, FP = fp)

  out <- ckf
  out$tpfp        <- tpfp
  out$mean        <- colMeans(tpfp)
  out$median      <- c(stats::median(tpfp[, 1]), stats::median(tpfp[, 2]))
  out$ratio_mean  <- mean(apply(tpfp, 1, function(z)
                          if (sum(z) == 0) 0 else if (z[2] == 0) 1 else z[1] / z[2]))
  out$derand_tpfp <- c(TP = sum(ckf$derand %in% truth), FP = sum(!ckf$derand %in% truth))
  out$se          <- c(TP = stats::sd(tp) / sqrt(ckf$R_used),
                       FP = stats::sd(fp) / sqrt(ckf$R_used))
  out$sd          <- c(TP = stats::sd(tp), FP = stats::sd(fp))
  out$n_truth     <- length(truth)
  ## W is not retained on the object, so how many of THIS truth set were testable
  ## cannot be recomputed here. NA rather than a number derived from something
  ## else, so it is never read as "none were present".
  out$n_truth_present <- NA_integer_
  out$pair        <- as.character(pair)
  out$reference_is_truth <- ckf$reference_taxon %in% truth
  out
}

################################################################################
## ckf_row / add_ckf_row -- the two lines that remove the merge downstream.
##
## `real data results.R` currently builds every study's final table by hand:
##
##   res_zemb$aldex_gold_full <- rbind(res_zemb$aldex_gold, CKF = ckf_zemb$mean, ...)
##
## which means the reporting script has to know each study's CKF filename, its
## shape (a list for Zemb, a list-of-matrices indexed by pair for Stammler, a
## list over gamma for Prasse) and which element to take. add_ckf_row() puts the
## row into res_tpfp$aldex_gold at analysis time instead, so the file that comes
## off disk already carries it and the reporting script just reads a row named
## "CKF" like any other method.
##
## `stat` chooses what the row reports:
##   "mean"    the average over replicates -- what the existing tables show
##   "derand"  the derandomised set's TP/FP -- one selection, not an average
##   "median"  the median over replicates
################################################################################
ckf_row <- function(ckf, stat = c("mean", "derand", "median"), label = "CKF") {
  stat <- match.arg(stat)
  v <- switch(stat,
              mean   = ckf$mean,
              median = ckf$median,
              derand = ckf$derand_tpfp)
  m <- matrix(as.numeric(v), nrow = 1,
              dimnames = list(label, c("TP", "FP")))
  m
}

add_ckf_row <- function(res_tpfp, ckf, stat = c("mean", "derand", "median"),
                        label = "CKF", also_derand = FALSE) {
  stat <- match.arg(stat)
  if (is.null(res_tpfp$aldex_gold))
    stop("res_tpfp has no $aldex_gold -- add_ckf_row() expects the object ",
         "summarize_tpfp() returns")

  ## Replace rather than append when the row is already there, so re-running the
  ## block in an interactive session cannot stack three CKF rows on the table.
  g <- res_tpfp$aldex_gold
  rows <- c(label, if (also_derand) paste0(label, "_derand"))
  g <- g[!rownames(g) %in% rows, , drop = FALSE]
  g <- rbind(g, ckf_row(ckf, stat, label))
  if (also_derand)
    g <- rbind(g, ckf_row(ckf, "derand", paste0(label, "_derand")))
  res_tpfp$aldex_gold <- g
  res_tpfp
}

################################################################################
## record_ckf -- write the run, in both the new and the legacy shape.
##
## Two files, for two readers:
##
##   <dir>/ckf_record.rds   one row per (dataset, pair), the same idempotent
##                          replace-not-append rule record_scale_model() uses, so
##                          re-running a study updates its row instead of
##                          duplicating it. Plus ckf_record.csv to read outside R.
##   <dir>/<file>           the legacy per-study object, in the legacy list shape,
##                          written only when `file` is given. This is what the
##                          current `real data results.R` lines already read
##                          (ckf_zemb$mean and so on), so nothing downstream has
##                          to change on the same day the analysis scripts do.
##
## `legacy_index` handles Stammler, whose legacy file is a LIST of TP/FP matrices
## indexed by pair number rather than a single object. Passing the pair's index
## slots this run into that list, leaving the other pairs alone.
################################################################################
record_ckf <- function(ckf, dir, file = NULL, legacy_index = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  ## ---- the keyed record --------------------------------------------------
  row <- data.frame(
    dataset = ckf$dataset, pair = ckf$pair,
    TP = unname(ckf$mean[1]),  FP = unname(ckf$mean[2]),
    sd_TP = unname(ckf$sd[1]), sd_FP = unname(ckf$sd[2]),
    se_TP = unname(ckf$se[1]), se_FP = unname(ckf$se[2]),
    median_TP = ckf$median[1], median_FP = ckf$median[2],
    derand_TP = unname(ckf$derand_tpfp[1]), derand_FP = unname(ckf$derand_tpfp[2]),
    n_derand = length(ckf$derand), derand_thresh = ckf$derand_thresh,
    ratio_mean = ckf$ratio_mean,
    n_empty = ckf$n_empty, empty_rate = ckf$empty_rate,
    R = ckf$R, R_used = ckf$R_used, q = ckf$q, offset = ckf$offset,
    penalty = ckf$penalty, seed = ckf$seed,
    D = ckf$D, N = ckf$N,
    n_truth = ckf$n_truth, n_truth_present = ckf$n_truth_present,
    reference_taxon = ckf$reference_taxon,
    reference_is_truth = ckf$reference_is_truth,
    elapsed_sec = ckf$elapsed_sec, run_at = ckf$run_at,
    stringsAsFactors = FALSE)

  f   <- file.path(dir, CKF_RECORD_FILE)
  old <- if (file.exists(f)) readRDS(f) else NULL
  if (!is.null(old)) {
    key <- function(z) paste0(z$dataset, "|",
                              ifelse(is.na(z$pair) | !nzchar(z$pair), "NA", z$pair))
    old <- old[!(key(old) %in% key(row)), , drop = FALSE]
    for (cn in setdiff(names(row), names(old))) old[[cn]] <- NA
    for (cn in setdiff(names(old), names(row))) row[[cn]] <- NA
    old <- old[, names(row), drop = FALSE]
  }
  rec <- rbind(old, row)
  saveRDS(rec, f)
  utils::write.csv(rec, file.path(dir, CKF_RECORD_CSV), row.names = FALSE)

  ## The selection sets, kept out of the summary because they are variable-length
  ## character vectors and two orders of magnitude larger. With these on disk,
  ## scoring against a revised gold standard needs no refit -- see score_ckf().
  fs   <- file.path(dir, "ckf_selected.rds")
  olds <- if (file.exists(fs)) readRDS(fs) else list()
  k    <- paste0(ckf$dataset, "|", if (is.na(ckf$pair) || !nzchar(ckf$pair)) "NA" else ckf$pair)
  olds[[k]] <- list(dataset = ckf$dataset, pair = ckf$pair, selected = ckf$selected,
                    freq = ckf$freq, derand = ckf$derand, q = ckf$q,
                    offset = ckf$offset, seed = ckf$seed, run_at = ckf$run_at)
  saveRDS(olds, fs)

  ## ---- the legacy per-study file -----------------------------------------
  if (!is.null(file)) {
    legacy <- list(tpfp = ckf$tpfp, mean = ckf$mean,
                   median = ckf$median, ratio_mean = ckf$ratio_mean)
    fl <- file.path(dir, file)
    if (is.null(legacy_index)) {
      saveRDS(legacy, fl)
    } else {
      ## Stammler's shape: a bare list of cbind(tp, fp) matrices, indexed by pair.
      prev <- if (file.exists(fl)) readRDS(fl) else list()
      if (!is.list(prev) || !is.null(names(prev))) prev <- list()
      prev[[legacy_index]] <- ckf$tpfp
      saveRDS(prev, fl)
    }
    message(sprintf("  [CKF] %s%s: TP %.2f FP %.2f -> %s",
                    ckf$dataset, if (is.na(ckf$pair)) "" else paste0(" / ", ckf$pair),
                    ckf$mean[1], ckf$mean[2], fl))
  }
  invisible(row)
}

## Read every dataset's CKF record back -- the companion to
## collect_scale_records() / collect_perp_records() in scale_model_record.R.
collect_ckf_records <- function(dirs) {
  rows <- lapply(dirs, function(d) {
    f <- file.path(d, CKF_RECORD_FILE)
    if (!file.exists(f)) return(NULL)
    readRDS(f)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    warning("no ", CKF_RECORD_FILE, " found -- re-run the DA scripts to create them",
            call. = FALSE)
    return(NULL)
  }
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) {
    for (cn in setdiff(cols, names(r))) r[[cn]] <- NA
    r[, cols, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

################################################################################
## ckf_reference_check: which taxon did MRCKF drop, in every recorded run?
##
## MRCKF builds its design as log(Z[,t] / W[,pstar]), so the LAST column of W has
## no design column and can never be selected. Which taxon that is depends on the
## row order of Y_raw, which depends on the study's filtering -- so it is not
## something to reason about, it is something to look up after the fact.
##
##   ckf_reference_check(c(DIR_SELEX, DIR_ZEMB, DIR_STAMMLER, DIR_JIN,
##                         DIR_PRASSE, DIR_VANDEPUTTE))
##
## One row per (dataset, pair): the dropped taxon and whether it was in that
## run's truth set. `in_truth` TRUE means that run's TP is understated by at
## least one and the run should be repeated with `reference = <another taxon>`.
## Reads the records only -- it fits nothing.
################################################################################
ckf_reference_check <- function(dirs) {
  rec <- collect_ckf_records(dirs)
  if (is.null(rec)) return(NULL)
  out <- rec[, c("dataset", "pair", "reference_taxon", "reference_is_truth",
                 "D", "n_truth", "n_truth_present", "run_at")]
  names(out)[names(out) == "reference_is_truth"] <- "in_truth"
  out <- out[order(out$dataset, out$pair), , drop = FALSE]
  rownames(out) <- NULL
  bad <- which(out$in_truth %in% TRUE)
  if (length(bad))
    warning("the CKF reference taxon is in the truth set for ", length(bad),
            " run(s): ",
            paste0(out$dataset[bad], "/", out$pair[bad], " (", out$reference_taxon[bad], ")",
                   collapse = ", "),
            " -- those TP counts are understated. Re-run with `reference =`.",
            call. = FALSE)
  out
}

################################################################################
## ckf_stability: how many replicates would have been enough?
##
## Answers the question R = 100 was chosen to avoid asking. For each candidate R,
## the Monte Carlo SE of the mean TP and FP implied by the replicates ALREADY
## run -- so it costs nothing and says whether the next run can be shorter.
##
##   ckf_stability(ckf)                 # at 10, 20, 30, 50, 100
##   ckf_stability(ckf, se_target=0.1)  # and the smallest R that reaches it
################################################################################
ckf_stability <- function(ckf, at = c(10, 20, 30, 50, 100), se_target = NULL) {
  R <- ckf$R_used
  at <- at[at <= R]
  if (!length(at)) return(NULL)
  d <- do.call(rbind, lapply(at, function(r) data.frame(
    R = r,
    mean_TP = mean(ckf$tpfp[1:r, 1]), se_TP = stats::sd(ckf$tpfp[1:r, 1]) / sqrt(r),
    mean_FP = mean(ckf$tpfp[1:r, 2]), se_FP = stats::sd(ckf$tpfp[1:r, 2]) / sqrt(r),
    empty   = sum(lengths(ckf$selected[1:r]) == 0) / r,
    stringsAsFactors = FALSE)))
  if (!is.null(se_target)) {
    ok <- which(d$se_TP <= se_target & d$se_FP <= se_target)
    attr(d, "R_enough") <- if (length(ok)) d$R[ok[1]] else NA_integer_
  }
  rownames(d) <- NULL
  d
}
