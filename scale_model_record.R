################################################################################
## scale_model_record.R -- one place where every real-data DA script records
## which scale model sparse SSRV selected, and how many taxa EACH arm would have
## called significant.
##
## Why this exists: a finished analysis used to leave no trace of what the
## unselected arm would have done. The two arms can differ by a lot
## (McMurrough2014: 16 vs 72 taxa), which makes the selection itself a
## reportable result rather than an implementation detail.
##
## sparse_auto() now returns all three results -- the hybrid rule's choice, the
## laplace scale model on its own and the bootstrap scale model on its own -- so
## this file records all three: the COUNTS in scale_model_record.rds, and the
## per-feature tables (which taxon, at what effect size, at what q) in
## sparse_arm_sig.rds and its CSV twin.
##
## Usage, immediately after the sparse_auto() / sparse_auto_sens() call:
##
##   record_scale_model(sparse_out, dataset = "Zemb2020", dir = DIR_ZEMB)
##   record_scale_model(sparse_out, dataset = "Stammler2016", dir = DIR_STAMMLER,
##                      pair = target_pool)          # multi-pair analyses
##
## Writes <dir>/scale_model_record.rds, a data.frame with one row per
## (dataset, pair). Re-running an analysis REPLACES its row rather than
## appending a duplicate, so the file always reflects the latest run of each.
##
## Alongside it, <dir>/sparse_arm_sig.rds keeps each run's per-arm feature
## tables, and <dir>/sparse_arm_features.csv is the same thing readable outside
## R -- one row per feature called by EITHER scale model, labelled shared /
## laplace_only / bootstrap_only.
##
## `real data results.R` reads these back with collect_scale_records() for the
## counts and collect_arm_features() for the identities.
################################################################################

RECORD_FILE <- "scale_model_record.rds"

## Companion to RECORD_FILE. scale_model_record.rds stores the arm COUNTS, which
## is enough to say whether the choice mattered but not enough to score either
## arm: TP/FP needs to know WHICH features an arm called, so it can be
## intersected with that study's truth set. This file keeps those identities,
## one entry per (dataset, pair), so the reporting scripts can put the bootstrap
## and laplace arms next to the selected one without refitting anything.
ARM_SIG_FILE <- "sparse_arm_sig.rds"

## Pull the scale-model fields off a sparse_auto() / sparse_auto_sens() result.
## Returns NA for anything the attributes do not carry, so a result produced by
## an older sparse_SSRV.R still yields a row instead of erroring.
scale_model_summary <- function(out) {
  at   <- function(nm, default = NA) if (is.null(attr(out, nm))) default else attr(out, nm)
  stat <- at("scale_stat", c(bootstrap = NA_real_, laplace = NA_real_))
  nsig <- at("n_sig",      c(bootstrap = NA_integer_, laplace = NA_integer_))
  sba  <- at("sig_by_arm", NULL)

  ## Overlap between the two arms' significant sets -- the quantity that says
  ## whether the choice actually mattered, as opposed to the counts alone.
  shared <- only_b <- only_l <- NA_integer_
  if (!is.null(sba) && !all(is.na(sba$bootstrap)) && !all(is.na(sba$laplace))) {
    sb <- names(which(sba$bootstrap)); sl <- names(which(sba$laplace))
    shared <- length(intersect(sb, sl))
    only_b <- length(setdiff(sb, sl))
    only_l <- length(setdiff(sl, sb))
  }

  data.frame(
    scale_model      = at("scale_model", NA_character_),
    scale_rule       = at("scale_rule",  NA_character_),
    stat_bootstrap   = unname(stat[1]),
    stat_laplace     = unname(stat[2]),
    n_sig_selected   = sum(out$sig.q, na.rm = TRUE),
    n_sig_bootstrap  = unname(nsig["bootstrap"]),
    n_sig_laplace    = unname(nsig["laplace"]),
    n_shared         = shared,
    only_bootstrap   = only_b,
    only_laplace     = only_l,
    n_laplace_fail   = at("n_laplace_fail", NA_integer_),
    n_var_reject     = at("n_var_reject",   NA_integer_),
    D                = nrow(out),
    stringsAsFactors = FALSE)
}

## One key spelling, used by both the writer and the readers. NA and the string
## "NA" have to collapse to the same key or a run recorded without a pair could
## never be looked up again.
arm_sig_key <- function(dataset, pair = NA_character_) {
  paste0(dataset, "|", if (is.na(pair) || !nzchar(pair)) "NA" else as.character(pair))
}

## Human-readable twin of ARM_SIG_FILE. The .rds is what the reporting scripts
## read; this is what gets opened in a spreadsheet when the question is "so which
## taxa does laplace call that bootstrap does not?".
ARM_FEATURE_CSV <- "sparse_arm_features.csv"

################################################################################
## sparse_arm_features: one row per feature EITHER scale model called, with the
## effect size and q from each arm side by side.
##
## Reads the per-arm tables sparse_auto() / sparse_auto_sens() attach
## (attr(out, "arm_out")), so it needs nothing from sparse_SSRV.R and works on a
## saved fit. `set` is the column that answers the question the counts could not:
##
##   shared          both scale models called it
##   laplace_only    only the laplace scale model did
##   bootstrap_only  only the bootstrap scale model did
##
## and it is NA -- not "the other arm only" -- whenever one arm is unusable, so a
## withheld arm is never read as a disagreement.
##
## all_features = TRUE returns every tested feature instead of the union of the
## two significant sets, which is what you want when comparing intervals rather
## than calls.
##
## An older result with no per-arm tables returns NULL with a warning; the
## reporting scripts treat that as "not re-run yet" rather than "found nothing".
################################################################################
sparse_arm_features <- function(out, dataset = NA_character_, pair = NA_character_,
                                all_features = FALSE) {
  arms <- attr(out, "arm_out")
  if (is.null(arms)) {
    warning("result for ", dataset, " carries no per-arm tables -- re-run the fit",
            call. = FALSE)
    return(NULL)
  }
  L <- arms$laplace; B <- arms$bootstrap
  ids <- rownames(L)
  usableL <- !all(is.na(L$sig.q)); usableB <- !all(is.na(B$sig.q))

  sigL <- !is.na(L$sig.q) & L$sig.q
  sigB <- !is.na(B$sig.q) & B$sig.q
  keep <- if (all_features) rep(TRUE, length(ids)) else (sigL | sigB)
  if (!any(keep)) keep <- logical(0)

  chosen <- attr(out, "scale_model")
  hyb    <- if (identical(chosen, "bootstrap")) sigB else sigL

  ## Both arms usable -> a real three-way split. One arm withheld -> NA, because
  ## "bootstrap_only" would otherwise be recorded for a comparison never made.
  set <- rep(NA_character_, sum(keep))
  if (usableL && usableB) {
    sl <- sigL[keep]; sb <- sigB[keep]
    set <- ifelse(sl & sb, "shared", ifelse(sl, "laplace_only", "bootstrap_only"))
  }

  ## rep(, n) rather than letting data.frame() recycle: when NEITHER arm called
  ## anything, n is 0 and a length-1 `dataset` against a length-0 `feature` is an
  ## error, not an empty table. A study where both scale models found nothing is a
  ## real result and has to come back as 0 rows.
  n <- sum(keep)
  d <- data.frame(
    dataset         = rep(dataset,            n),
    pair            = rep(as.character(pair), n),
    feature         = ids[keep],
    set             = set,
    selected_arm    = rep(chosen,             n),
    in_hybrid       = hyb[keep],
    in_laplace      = if (usableL) sigL[keep] else NA,
    in_bootstrap    = if (usableB) sigB[keep] else NA,
    mean_laplace    = L$mean[keep],
    lower_laplace   = L[[1]][keep],
    upper_laplace   = L[[2]][keep],
    q_laplace       = L$q.val[keep],
    mean_bootstrap  = B$mean[keep],
    lower_bootstrap = B[[1]][keep],
    upper_bootstrap = B[[2]][keep],
    q_bootstrap     = B$q.val[keep],
    stringsAsFactors = FALSE)

  ## shared first, then each arm's own calls, most significant first inside each
  ## group -- the order you want to read the disagreement in.
  ord <- order(match(d$set, c("shared", "laplace_only", "bootstrap_only")),
               pmin(d$q_laplace, d$q_bootstrap, na.rm = TRUE), d$feature)
  d <- d[ord, , drop = FALSE]
  rownames(d) <- NULL
  d
}

record_scale_model <- function(out, dataset, dir, pair = NA_character_,
                               S = NA_integer_, N = NA_integer_, notes = NA_character_) {
  row <- cbind(data.frame(dataset = dataset, pair = as.character(pair),
                          stringsAsFactors = FALSE),
               scale_model_summary(out),
               data.frame(N = N, S = S, notes = notes,
                          run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                          stringsAsFactors = FALSE))

  f   <- file.path(dir, RECORD_FILE)
  old <- if (file.exists(f)) readRDS(f) else NULL
  ## Replace the row for this (dataset, pair) so re-running is idempotent.
  if (!is.null(old)) {
    dup <- old$dataset == dataset &
           ifelse(is.na(old$pair), "NA", old$pair) == ifelse(is.na(row$pair), "NA", row$pair)
    old <- old[!dup, , drop = FALSE]
    ## If that removed the LAST row, `old` is a 0-row data.frame, and the column
    ## alignment below then does `old[[cn]] <- NA` on it -- which errors with
    ## "replacement has 1 row, data has 0" rather than adding an empty column.
    ## That fires exactly when a study is re-run into a record file that holds
    ## only its own rows AND the schema has grown since (Vandeputte, 2026-09-01:
    ## the August record predates obs_S / bias_S / covered_S / dlog2_r). Dropping
    ## to NULL is the same thing the no-file branch above does, and rbind(NULL, x)
    ## is x, so the re-run simply rewrites the file from scratch.
    if (nrow(old) == 0L) old <- NULL
  }
  if (!is.null(old)) {
    ## A column added to this file later would otherwise fail to rbind against
    ## rows written by an earlier version; align on the union of both.
    for (cn in setdiff(names(row), names(old))) old[[cn]]  <- NA
    for (cn in setdiff(names(old), names(row))) row[[cn]]  <- NA
    old <- old[, names(row), drop = FALSE]
  }
  saveRDS(rbind(old, row), f)

  ## Arm identities, keyed the same way. Written as a named list rather than
  ## appended to the data.frame above because the entries are variable-length
  ## character vectors, which a data.frame column cannot hold cleanly.
  ##
  ## An arm whose draws all failed yields all-NA from sparse_auto() rather than
  ## a false empty set; that is stored as NULL, which the reader turns into NA
  ## TP/FP instead of a misleading 0/0.
  sba <- attr(out, "sig_by_arm")
  if (!is.null(sba)) {
    arm_ids <- function(v) if (is.null(v) || all(is.na(v))) NULL else names(which(v))
    fa  <- file.path(dir, ARM_SIG_FILE)
    olda <- if (file.exists(fa)) readRDS(fa) else list()
    ## `features` is the per-feature table -- ids plus each arm's effect size and
    ## q -- and is NULL for a fit produced before sparse_auto() returned all three
    ## arms. The bare id vectors above are kept alongside it: arm_tpfp() and
    ## all_pairs_results.R read those, and they still work on old records.
    olda[[arm_sig_key(dataset, pair)]] <- list(
      dataset   = dataset,
      pair      = as.character(pair),
      selected  = row$scale_model,
      bootstrap = arm_ids(sba$bootstrap),
      laplace   = arm_ids(sba$laplace),
      features  = suppressWarnings(sparse_arm_features(out, dataset, pair)),
      run_at    = row$run_at)
    saveRDS(olda, fa)

    ## Rewritten in full each time from the .rds, so the CSV can never describe a
    ## run the .rds has already replaced.
    feats <- do.call(rbind, lapply(olda, `[[`, "features"))
    if (!is.null(feats) && nrow(feats))
      write.csv(feats, file.path(dir, ARM_FEATURE_CSV), row.names = FALSE)
  }

  message(sprintf("  [scale model] %s%s: chose %s | n_sig L=%s B=%s (shared %s) -> %s",
                  dataset, if (is.na(pair)) "" else paste0(" / ", pair),
                  row$scale_model, row$n_sig_laplace, row$n_sig_bootstrap,
                  row$n_shared, f))
  invisible(row)
}

## Read every dataset's record back. Used by `real data results.R`.
collect_scale_records <- function(dirs) {
  rows <- lapply(dirs, function(d) {
    f <- file.path(d, RECORD_FILE)
    if (!file.exists(f)) return(NULL)
    readRDS(f)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    warning("no ", RECORD_FILE, " found -- re-run the DA scripts to create them",
            call. = FALSE)
    return(NULL)
  }
  ## Union the columns so records written by different versions still stack.
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) { for (cn in setdiff(cols, names(r))) r[[cn]] <- NA; r[, cols, drop = FALSE] })
  do.call(rbind, rows)
}

## Read back one run's per-arm significant sets. `pair` may be omitted when the
## study recorded exactly one run; when several exist the caller has to say
## which, because silently picking one would attribute another comparison's
## numbers to this row.
read_arm_sig <- function(dir, dataset, pair = NA_character_) {
  f <- file.path(dir, ARM_SIG_FILE)
  if (!file.exists(f)) {
    warning("no ", ARM_SIG_FILE, " in ", dir,
            " -- re-run that DA script to record the per-arm sets", call. = FALSE)
    return(NULL)
  }
  all <- readRDS(f)
  if (is.na(pair)) {
    hit <- Filter(function(e) identical(e$dataset, dataset), all)
    if (length(hit) == 1) return(hit[[1]])
    if (length(hit) == 0) {
      warning("no arm record for dataset ", dataset, " in ", f, call. = FALSE)
      return(NULL)
    }
    warning(dataset, " has ", length(hit), " recorded runs (",
            paste(vapply(hit, `[[`, "", "pair"), collapse = ", "),
            ") -- name one with `pair`", call. = FALSE)
    return(NULL)
  }
  e <- all[[arm_sig_key(dataset, pair)]]
  if (is.null(e))
    warning("no arm record for ", arm_sig_key(dataset, pair), " in ", f, call. = FALSE)
  e
}

## TP/FP for both arms against a study's truth set, shaped exactly like the
## rows the reporting scripts rbind onto aldex_gold: two columns, TP then FP,
## rownames sparse_SSRV_B / sparse_SSRV_L. Returns NA rows rather than nothing
## when the record is missing, so a table still renders before the DA scripts
## have been re-run -- a withheld number, not a zero.
arm_tpfp <- function(dir, dataset, truth, pair = NA_character_,
                     labels = c(bootstrap = "sparse_SSRV_B", laplace = "sparse_SSRV_L")) {
  e <- suppressWarnings(read_arm_sig(dir, dataset, pair))
  one <- function(ids) {
    if (is.null(e) || is.null(ids)) return(c(TP = NA_integer_, FP = NA_integer_))
    c(TP = sum(ids %in% truth), FP = sum(!ids %in% truth))
  }
  m <- rbind(one(if (is.null(e)) NULL else e$bootstrap),
             one(if (is.null(e)) NULL else e$laplace))
  rownames(m) <- unname(labels[c("bootstrap", "laplace")])
  colnames(m) <- c("TP", "FP")
  m
}

################################################################################
## read_arm_features / collect_arm_features -- the identity side of the record.
##
## collect_scale_records() answers "how many did each scale model call?".
## These answer "WHICH ones?", which is the question the counts raise and cannot
## settle: two arms can agree on a number and disagree on every feature in it.
##
## read_arm_features()    one run, resolved the same way read_arm_sig() resolves
##                        a missing `pair` (unique match, or say which).
## collect_arm_features() every run in `dirs`, stacked, with dataset and pair
##                        columns -- what `real data results.R` reports.
##
## Both return NULL rather than erroring when nothing has been recorded yet, so
## a reporting script still runs before every DA script has been re-run.
################################################################################
read_arm_features <- function(dir, dataset, pair = NA_character_) {
  e <- read_arm_sig(dir, dataset, pair)
  if (is.null(e)) return(NULL)
  if (is.null(e$features)) {
    warning("the record for ", arm_sig_key(dataset, pair), " has arm ids but no ",
            "per-feature table -- it predates them; re-run that DA script",
            call. = FALSE)
    return(NULL)
  }
  e$features
}

collect_arm_features <- function(dirs) {
  rows <- unlist(lapply(dirs, function(d) {
    f <- file.path(d, ARM_SIG_FILE)
    if (!file.exists(f)) return(NULL)
    lapply(readRDS(f), `[[`, "features")
  }), recursive = FALSE)
  ## NULL means "that run predates the per-feature tables". A 0-row table means
  ## "both scale models called nothing", which is a result and is kept -- dropping
  ## it here would make an empty study indistinguishable from an un-re-run one.
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    warning("no per-feature arm records found in any of the given dirs -- re-run ",
            "the DA scripts to write them", call. = FALSE)
    return(NULL)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

## Counts implied by an arm-feature table, one row per (dataset, pair). Same
## three numbers scale_model_record.rds already stores, recomputed here from the
## identities -- so a disagreement between this and collect_scale_records() means
## the two files came from different runs.
arm_feature_counts <- function(feats) {
  if (is.null(feats) || !nrow(feats)) return(NULL)
  key <- paste(feats$dataset, feats$pair, sep = "|")
  do.call(rbind, lapply(split(feats, key), function(d) data.frame(
    dataset        = d$dataset[1],
    pair           = d$pair[1],
    selected_arm   = d$selected_arm[1],
    n_union        = nrow(d),
    n_shared       = sum(d$set %in% "shared"),
    only_laplace   = sum(d$set %in% "laplace_only"),
    only_bootstrap = sum(d$set %in% "bootstrap_only"),
    stringsAsFactors = FALSE)))
}

################################################################################
## report_arm_features / arm_sig_tax -- the two lines a DA script needs.
##
## Every real-data script wants the same thing after its sparse_auto() call: see
## which features each scale model called, and be able to carry the two arms into
## whatever table it builds next. These are that, spelled once instead of
## copy-pasted seven times.
##
##   arm_report <- report_arm_features(sparse_out, "Prasse2015")
##   tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")     # ready for
##   tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")   # sig_features_long()
##
## report_arm_features() only reads attributes -- it fits nothing and writes
## nothing -- so it is safe to re-run on a saved fit.
################################################################################
report_arm_features <- function(out, dataset, pair = NA_character_, max_print = 40) {
  feats <- suppressWarnings(sparse_arm_features(out, dataset, pair))
  arms  <- attr(out, "arm_out")
  lab   <- paste0(dataset, if (is.na(pair)) "" else paste0(" / ", pair))
  if (is.null(feats)) {
    message("  [arms] ", lab, ": this fit carries no per-arm tables -- re-run it")
    return(invisible(NULL))
  }
  cat("\n== features called by each scale model (", lab, ") ==\n", sep = "")
  ## From the attribute, not from feats[1, ]: a run where NEITHER arm called
  ## anything has 0 rows, and the selected arm is still a fact about that run.
  cat("   selected arm: ", attr(out, "scale_model"),
      "  |  shared ",         sum(feats$set %in% "shared"),
      ", laplace only ",      sum(feats$set %in% "laplace_only"),
      ", bootstrap only ",    sum(feats$set %in% "bootstrap_only"), "\n", sep = "")
  show <- c("feature", "set", "in_hybrid", "mean_laplace", "q_laplace",
            "mean_bootstrap", "q_bootstrap")
  print(head(feats[, show, drop = FALSE], max_print), row.names = FALSE)
  if (nrow(feats) > max_print)
    cat("   ... ", nrow(feats) - max_print, " more; full table in ", ARM_FEATURE_CSV,
        "\n", sep = "")
  sig <- function(tab) if (is.null(tab)) NULL else tab[tab$sig.q %in% TRUE, , drop = FALSE]
  invisible(list(features   = feats,
                 laplace    = sig(arms$laplace),
                 bootstrap  = sig(arms$bootstrap)))
}

## One arm's significant rows, shaped like the tax_* tables the analysis scripts
## hand to sig_features_long(): a `category` column carrying the feature id.
## Returns a 0-row frame rather than NULL when an arm called nothing, so an arm
## that found nothing still appears as an all-zero column instead of vanishing
## from the incidence matrix.
arm_sig_tax <- function(arm_report, arm = c("laplace", "bootstrap")) {
  arm <- match.arg(arm)
  tab <- if (is.null(arm_report)) NULL else arm_report[[arm]]
  if (is.null(tab))
    tab <- data.frame(`2.5%` = numeric(0), `97.5%` = numeric(0), mean = numeric(0),
                      p.val = numeric(0), q.val = numeric(0), sig.q = logical(0),
                      check.names = FALSE)
  tab$category  <- rownames(tab)
  rownames(tab) <- NULL
  tab
}

################################################################################
## theta^perp -- recording the estimated SCALE alongside the selected arm.
##
## scale_model_record.rds says which scale model was chosen and what each arm
## called. It says nothing about the number the scale model actually produces:
## theta^perp, the common log2 shift that turns relative log fold changes into
## absolute ones. That number is the method's central claim, and on the three
## studies with external quantification -- Zemb (qPCR), Jin (absolute totals),
## Vandeputte (flow cytometry) -- it is directly measurable, so it can be scored
## rather than only asserted.
##
## Written, one row per (dataset, pair, arm):
##
##   theta_perp_record.rds / .csv   the estimates -- posterior mean, sd, median,
##                                  95% interval, the same interval as a fold
##                                  change on the total load, and the Monte Carlo
##                                  error of the mean
##   theta_perp_draws.rds           the posterior DRAWS, so coverage at any level
##                                  and any plot of estimate-vs-measurement can be
##                                  made later without refitting anything
##
## Usage, right after record_scale_model():
##
##   record_theta_perp(sparse_out, dataset = "Vandeputte2017", dir = DIR_VAN,
##                     S = S, N = N,
##                     obs = mean(log2(cc[X == 1])) - mean(log2(cc[X == 0])))
##
## `obs` is optional and is the EXTERNAL measurement of the same quantity, on the
## same log2 scale and in the same direction (treatment minus control). Supply it
## and the row carries bias, absolute error and whether the credible interval
## covered it; omit it and the estimate is still recorded, just unscored.
##
## `matched` is the list from perp_matched_obs(), and adds obs_S -- that same
## measurement projected onto the FEATURES THE MODEL WAS FITTED ON, which is the
## quantity theta^perp actually estimates. Both references are recorded side by
## side and neither replaces the other: obs is the physical measurement, obs_S is
## the like-for-like comparison. Read the block above perp_matched_obs() before
## quoting either -- in particular, what is treated as outside the retained set
## changes obs_S by more than the filter itself does.
##
## IT MUST BE A DENSITY -- total abundance per unit of sample mass or volume, not
## total abundance per tube. The sequenced composition is invariant to how much
## material was loaded (the mass cancels out of w_dn = d_dn*mg / (mg*sum d_dn)),
## so theta^perp is the log2 ratio of total DENSITY and nothing else. Handing this
## a per-tube total scores the estimate against the pipetting: on Zemb, where the
## spiked tubes hold 11-14 mg and the unspiked ones 38-55 mg, the two differ by
## 1.9 log2 and flip the sign of the answer. Jin's copies/mg and Vandeputte's
## cells-per-gram are already densities; a per-tube figure has to be divided by
## the sample mass before it comes in here.
################################################################################

PERP_FILE       <- "theta_perp_record.rds"
PERP_CSV        <- "theta_perp_record.csv"
PERP_DRAWS_FILE <- "theta_perp_draws.rds"

################################################################################
## perp_matched_obs: the external measurement, projected onto the features the
## model was actually fitted on.
##
## `obs` above measures the WHOLE community. qPCR, flow cytometry and spike-ins
## all quantify what went into the tube, before anyone filtered a count table.
## theta^perp as ESTIMATED is the log2 fold change of the total load of the
## RETAINED features only: sparse_auto() is handed log2(rdirichlet(Y0 + 0.5)) --
## proportions renormalised WITHIN the kept set -- so
## log2 W_dn = log2 p^S_dn + log2 T_S,n and the shift the mode identifies is
## Delta log2 T_S, not Delta log2 T. The two are different numbers, and every
## bias / coverage figure recorded before this function existed compared them as
## though they were the same one.
##
## With S the retained features, A the community the assay measures, and
## r_n = T_S,n / T_n the retained fraction of the load,
##
##     T_S,n = r_n T_n    =>    theta^perp_S = Delta log2 T_n + Delta log2 r_n
##
## That is an IDENTITY, not a first-order correction, because Delta is a
## difference of GROUP MEANS OF LOGS and log T_S,n = log T_n + log r_n holds
## sample by sample. Define Delta the other way round -- log of a ratio of means
## -- and it stops being exact.
##
## This function does not evaluate that formula. It forms W_dn = p_dn * T^obs_n
## across ALL of A and sums over S, which is algebraically the same thing but is
## the definition applied as written, so it cannot drift out of step with the
## quantity it is supposed to compute.
##
## WHAT COUNTS AS "OUTSIDE S" IS A JUDGEMENT AND IT DOMINATES THE ANSWER.
## Measured on the studies here:
##
##   Zemb        abundance filter, no agglomeration   r_n 0.984-0.991   -0.005
##   Vandeputte  glom to genus + prevalence filter    r_n 0.167-0.986   -0.135
##   Prasse      glom to genus + prevalence filter    r_n 0.117-0.342   -0.328
##
## The rare-taxon filters are second order BY CONSTRUCTION -- they drop features
## carrying almost no reads. What actually moves the number is tax_glom(), which
## silently discards every OTU with no assignment at the target rank: 22381 OTUs
## become 561 genera on Prasse and seven-eighths of the reads leave with them.
## So pass `Y_all` as the table BEFORE agglomeration, since those reads are
## community that the assay counted -- and know that the matched reference then
## inherits the annotation database's coverage, which is an uncomfortable
## dependency for a value being used as ground truth. Passing the post-glom table
## instead scores against the prevalence filter alone; on Prasse that is +0.080
## rather than -0.328, so the choice is not cosmetic and belongs in the caption.
##
## r_n IS A READ FRACTION. For a 16S qPCR assay (Zemb, Prasse, Jin) that is close
## to the copy fraction it needs to be. For Vandeputte's flow cytometry the
## target is a CELL fraction, so 16S copy-number variation and extraction bias
## enter the correction -- the one study where this plug-in is on the weakest
## footing is also the one where it makes the bias slightly worse.
##
## THE RETAINED SET IS GIVEN ONE OF TWO WAYS, and agglomeration decides which.
##
##   keep    ids / logical / positions INTO `Y_all`. Correct only when the fitted
##           features are literally rows of `Y_all` -- i.e. no tax_glom().
##   Y_kept  the fitted count table itself. The right form after agglomeration:
##           tax_glom() names each genus after its most abundant constituent OTU,
##           so taxa_names(red) ARE valid ids in the pre-glom table but select
##           only those archetype OTUs. Summing them would compute the retained
##           fraction over 29 OTUs where the model was fitted on 29 GENERA, and
##           would do it silently, since every id resolves.
##
##   m <- perp_matched_obs(total = cc, X = X_0, Y_all = otu_before_glom,
##                         Y_kept = Y_0)
##   record_theta_perp(out, ..., matched = m)
##
## `total` must be the DENSITY, on the same footing as `obs` -- see the note
## above. Samples whose total is missing or non-positive are dropped from BOTH
## halves together, so obs and obs_S always describe the same sample set; that
## subsumes the is.finite() guard the Zemb script used to apply by hand.
##
## Returns a list: obs (whole community), obs_S (matched), dlog2_r, the range of
## r_n, the two feature counts, and how many samples were dropped.
################################################################################
perp_matched_obs <- function(total, X, Y_all, keep = NULL, Y_kept = NULL,
                             taxa_are_rows = TRUE,
                             kept_taxa_are_rows = taxa_are_rows) {
  M <- as.matrix(Y_all)
  if (!taxa_are_rows) M <- t(M)
  storage.mode(M) <- "double"

  if (is.null(keep) == is.null(Y_kept))
    stop("perp_matched_obs: give exactly one of `keep` (ids into Y_all) or ",
         "`Y_kept` (the fitted count table). After tax_glom() it must be Y_kept ",
         "-- see the note above.")

  if (!is.null(Y_kept)) {
    K <- as.matrix(Y_kept)
    if (!kept_taxa_are_rows) K <- t(K)
    storage.mode(K) <- "double"
    if (ncol(K) != ncol(M))
      stop("perp_matched_obs: `Y_kept` has ", ncol(K), " samples and `Y_all` has ",
           ncol(M), ".")
    if (!is.null(colnames(K)) && !is.null(colnames(M)) &&
        !identical(as.character(colnames(K)), as.character(colnames(M))))
      stop("perp_matched_obs: `Y_kept` and `Y_all` are in different sample orders.")
    n_kept <- nrow(K)
    tot_S_raw <- colSums(K)
  } else {
    ## `keep` as ids, as a logical over features, or as positions. Ids are the
    ## form to prefer: they are the only one that fails loudly when `Y_all` is not
    ## the table the filter was applied to.
    idx <- if (is.character(keep)) {
      miss <- setdiff(keep, rownames(M))
      if (length(miss))
        stop("perp_matched_obs: ", length(miss), " retained feature(s) are not rows of ",
             "`Y_all` (", paste(utils::head(miss, 5), collapse = ", "),
             if (length(miss) > 5) ", ..." else "", ").\n",
             "  `Y_all` must be the count table BEFORE filtering, carrying the same ",
             "feature ids as the fitted set. If the fit was agglomerated, pass ",
             "`Y_kept =` instead of `keep =`.")
      match(keep, rownames(M))
    } else if (is.logical(keep)) {
      if (length(keep) != nrow(M))
        stop("perp_matched_obs: logical `keep` has length ", length(keep),
             " but `Y_all` has ", nrow(M), " features.")
      which(keep)
    } else {
      as.integer(keep)
    }
    if (!length(idx)) stop("perp_matched_obs: `keep` selects no features.")
    n_kept <- length(idx)
    tot_S_raw <- colSums(M[idx, , drop = FALSE])
  }

  if (length(total) != ncol(M) || length(X) != ncol(M))
    stop("perp_matched_obs: `total` (", length(total), ") and `X` (", length(X),
         ") must both be one per SAMPLE of `Y_all` (", ncol(M), " columns). ",
         "If `Y_all` is samples x taxa, pass taxa_are_rows = FALSE.")
  ## Positional alignment is the failure this cannot detect on its own, so check
  ## names whenever both sides carry them.
  if (!is.null(names(total)) && !is.null(colnames(M)) &&
      !identical(as.character(names(total)), as.character(colnames(M))))
    stop("perp_matched_obs: names(total) and colnames(Y_all) differ -- the ",
         "measurement and the counts are in different sample orders.")

  X <- as.numeric(X)
  if (!setequal(unique(X), c(0, 1)))
    stop("perp_matched_obs: `X` must be 0/1 with both groups present.")

  ## Drop samples with no usable measurement from BOTH halves at once. Doing it
  ## here rather than at the call site is what keeps obs and obs_S over the same
  ## samples; scoring a whole-community value on 9 samples against a matched one
  ## on 8 would be exactly the mismatch this function was written to remove.
  ok <- is.finite(total) & total > 0
  n_dropped <- sum(!ok)
  if (n_dropped) {
    message("  [theta^perp] ", n_dropped, " sample(s) have no measured density -- ",
            "dropped from BOTH obs and obs_S")
    M <- M[, ok, drop = FALSE]; total <- total[ok]; X <- X[ok]
    tot_S_raw <- tot_S_raw[ok]
    if (!setequal(unique(X), c(0, 1)))
      stop("perp_matched_obs: dropping unmeasured samples left only one group.")
  }

  tot_all <- colSums(M)
  tot_S   <- tot_S_raw
  if (any(tot_all <= 0))
    stop("perp_matched_obs: sample(s) with zero total reads in `Y_all`.")
  if (any(tot_S <= 0))
    stop("perp_matched_obs: sample(s) with zero reads among the retained ",
         "features -- log2 r_n is undefined there.")
  ## r_n > 1 is not a rounding problem, it is the wrong `Y_all`: the retained
  ## features are a SUBSET of the community, so their reads cannot outnumber it.
  ## This is what catches a post-glom table passed as `Y_all`, or `Y_kept` and
  ## `Y_all` coming from different rarefactions.
  if (any(tot_S > tot_all * (1 + 1e-9)))
    stop("perp_matched_obs: retained reads exceed total reads in ",
         sum(tot_S > tot_all * (1 + 1e-9)), " sample(s) -- `Y_all` is not a ",
         "superset of the fitted table.")

  r   <- tot_S / tot_all
  W_S <- r * total                 # sum_{d in S} p_dn * T^obs_n  ==  r_n * T_n
  dl  <- function(v) mean(log2(v[X == 1])) - mean(log2(v[X == 0]))

  structure(list(obs       = dl(total),
                 obs_S     = dl(W_S),
                 dlog2_r   = dl(r),
                 r         = r,
                 r_min     = min(r),
                 r_max     = max(r),
                 D_all     = nrow(M),
                 D_kept    = n_kept,
                 n_dropped = n_dropped),
            class = "perp_matched")
}

## Three rows -- hybrid, laplace, bootstrap -- labelled with the run they came
## from, and scored against `obs` when one is given.
##
## Returns NULL for a fit produced before theta^perp was recorded, so an old
## saved result still loads instead of erroring.
## Both references are carried, never one in place of the other:
##
##   obs    the measurement as made -- the whole community, what the assay
##          actually quantified. Kept because it is the physical number, and
##          because every figure recorded before obs_S existed is scored on it.
##   obs_S  the same measurement projected onto the retained features, which is
##          the quantity theta^perp actually estimates. See perp_matched_obs().
##
## The gap between them (dlog2_r) is a filter diagnostic worth reading on its
## own: near zero means the filter took only features the assay barely counted.
perp_record_rows <- function(out, dataset = NA_character_, pair = NA_character_,
                             obs = NA_real_, matched = NULL) {
  tab <- suppressWarnings(sparse_perp(out))
  if (is.null(tab)) {
    warning("result for ", dataset, " carries no theta^perp -- re-run the fit",
            call. = FALSE)
    return(NULL)
  }
  scal <- function(v) if (length(v) != 1 || !is.finite(v)) NA_real_ else as.numeric(v)
  obs <- scal(obs)

  obs_S <- dlog2_r <- r_min <- r_max <- NA_real_
  D_all <- D_kept <- NA_integer_
  if (!is.null(matched)) {
    if (!is.list(matched) || !all(c("obs", "obs_S", "dlog2_r") %in% names(matched)))
      stop("record_theta_perp: `matched` must be the list returned by ",
           "perp_matched_obs(), not a bare number.")
    obs_S   <- scal(matched$obs_S)
    dlog2_r <- scal(matched$dlog2_r)
    r_min   <- scal(matched$r_min)
    r_max   <- scal(matched$r_max)
    D_all   <- as.integer(matched$D_all)
    D_kept  <- as.integer(matched$D_kept)
    ## The cross-check that earns this argument its keep. The caller's own `obs`
    ## line and the one perp_matched_obs() derives from `total` are computed from
    ## the same measurement, so they can only disagree if the two are looking at
    ## different samples -- which is precisely the alignment bug that would make
    ## obs_S wrong in a way no summary would show.
    mo <- scal(matched$obs)
    if (is.na(obs)) {
      obs <- mo
    } else if (!is.na(mo) && abs(obs - mo) > 1e-8) {
      stop(sprintf(paste0("record_theta_perp: obs = %+0.6f but perp_matched_obs() ",
                          "computes %+0.6f from the same measurement.\n  The two are ",
                          "over different samples -- check the sample subsetting on ",
                          "either side before trusting obs_S."), obs, mo))
    }
  }

  d <- cbind(data.frame(dataset = rep(dataset, nrow(tab)),
                        pair    = rep(as.character(pair), nrow(tab)),
                        stringsAsFactors = FALSE),
             tab)
  d$obs     <- obs
  d$bias    <- d$mean - obs
  d$abs_err <- abs(d$bias)
  ## NA, not FALSE, when nothing was measured -- an unscored row must never read
  ## as a miss.
  d$covered <- if (is.na(obs)) NA else (obs >= d$lower & obs <= d$upper)
  ## How many posterior sds the measurement sits from the estimate. Scale-free, so
  ## it is the one column that compares across datasets whose loads differ by
  ## orders of magnitude.
  d$z       <- d$bias / d$sd

  ## The same five columns against the filter-matched reference. Suffixed rather
  ## than replacing the originals, so a table can report both and nothing that
  ## already reads `bias` silently changes meaning.
  d$obs_S     <- obs_S
  d$bias_S    <- d$mean - obs_S
  d$abs_err_S <- abs(d$bias_S)
  d$covered_S <- if (is.na(obs_S)) NA else (obs_S >= d$lower & obs_S <= d$upper)
  d$z_S       <- d$bias_S / d$sd

  ## What the projection was worth, so the two references can be told apart
  ## without re-deriving them: obs_S = obs + dlog2_r, exactly.
  d$dlog2_r <- dlog2_r
  d$r_min   <- r_min
  d$r_max   <- r_max
  d$D_all   <- D_all
  d$D_kept  <- D_kept
  rownames(d) <- NULL
  d
}

## Print the three arms as a DA script's run log should show them, and return the
## rows invisibly. Reads attributes only -- fits nothing, writes nothing -- so it
## is safe to call on a saved fit.
report_theta_perp <- function(out, dataset, pair = NA_character_, obs = NA_real_,
                              matched = NULL) {
  d <- perp_record_rows(out, dataset, pair, obs, matched)
  lab <- paste0(dataset, if (is.na(pair)) "" else paste0(" / ", pair))
  if (is.null(d)) {
    message("  [theta^perp] ", lab, ": not recorded in this fit -- re-run it")
    return(invisible(NULL))
  }
  cat("\n== estimated theta^perp (", lab, ") ==\n", sep = "")
  cat("   log2 fold change of TOTAL load, treatment - control; ",
      "selected arm: ", d$selected[1], "\n", sep = "")
  show <- c("arm", "mean", "sd", "lower", "upper", "fold_change", "mc_se")
  if (!is.na(d$obs[1]))   show <- c(show, "obs", "bias", "covered")
  if (!is.na(d$obs_S[1])) show <- c(show, "obs_S", "bias_S", "covered_S")
  print(d[, show, drop = FALSE], row.names = FALSE, digits = 4)
  ## Spelled out rather than left to be inferred from two columns: obs is the
  ## whole community the assay measured, obs_S the part of it the model was
  ## fitted on, and the line between them says how much of the load the filter
  ## removed. A dlog2_r near zero means the two references are interchangeable
  ## and the distinction can be dropped from the caption; a large one means it
  ## cannot.
  if (!is.na(d$obs_S[1]))
    cat(sprintf(paste0("   obs is the whole community (D = %d); obs_S is the %d ",
                       "retained features:\n   r_n in [%.4f, %.4f], ",
                       "Delta log2 r = %+0.4f, so obs %+0.4f -> obs_S %+0.4f\n"),
                d$D_all[1], d$D_kept[1], d$r_min[1], d$r_max[1],
                d$dlog2_r[1], d$obs[1], d$obs_S[1]))
  invisible(d)
}

record_theta_perp <- function(out, dataset, dir, pair = NA_character_,
                              S = NA_integer_, N = NA_integer_, obs = NA_real_,
                              matched = NULL, notes = NA_character_, draws = TRUE) {
  d <- report_theta_perp(out, dataset, pair, obs, matched)
  if (is.null(d)) return(invisible(NULL))

  d$N      <- N
  d$S      <- S
  d$notes  <- notes
  d$run_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  f   <- file.path(dir, PERP_FILE)
  old <- if (file.exists(f)) readRDS(f) else NULL
  ## Replace this (dataset, pair)'s three rows, so re-running is idempotent --
  ## the same rule record_scale_model() uses, and the same NA/"NA" collapse, so a
  ## run recorded without a pair can still be found.
  if (!is.null(old)) {
    key_of <- function(z) paste0(z$dataset, "|", ifelse(is.na(z$pair) | !nzchar(z$pair), "NA", z$pair))
    old <- old[!(key_of(old) %in% key_of(d)[1]), , drop = FALSE]
    ## If that removed the LAST row, `old` is a 0-row data.frame, and the column
    ## alignment below then does `old[[cn]] <- NA` on it -- which errors with
    ## "replacement has 1 row, data has 0" rather than adding an empty column.
    ## That fires exactly when a study is re-run into a record file that holds
    ## only its own rows AND the schema has grown since (Vandeputte, 2026-09-01:
    ## the August record predates obs_S / bias_S / covered_S / dlog2_r). Dropping
    ## to NULL is the same thing the no-file branch above does, and rbind(NULL, x)
    ## is x, so the re-run simply rewrites the file from scratch.
    if (nrow(old) == 0L) old <- NULL
  }
  if (!is.null(old)) {
    for (cn in setdiff(names(d),   names(old))) old[[cn]] <- NA
    for (cn in setdiff(names(old), names(d)))   d[[cn]]   <- NA
    old <- old[, names(d), drop = FALSE]
  }
  rec <- rbind(old, d)
  saveRDS(rec, f)
  write.csv(rec, file.path(dir, PERP_CSV), row.names = FALSE)

  ## The draws themselves. Kept in their own file rather than in the record
  ## because they are two orders of magnitude larger than the summary, and the
  ## summary is the thing that gets opened. With these on disk, coverage at a
  ## level other than 95% -- or a plot of the posterior against the measured
  ## value -- needs no refit.
  if (isTRUE(draws)) {
    perp <- attr(out, "theta_perp")
    fd   <- file.path(dir, PERP_DRAWS_FILE)
    oldd <- if (file.exists(fd)) readRDS(fd) else list()
    ## Both references travel with the draws, or coverage recomputed at another
    ## level from this file could only ever be scored against the whole-community
    ## number -- the one the draws are NOT comparable with.
    oldd[[arm_sig_key(dataset, pair)]] <- list(
      dataset = dataset, pair = as.character(pair),
      selected = perp$scale_model, obs = d$obs[1], obs_S = d$obs_S[1],
      dlog2_r = d$dlog2_r[1],
      laplace = perp$laplace, bootstrap = perp$bootstrap,
      S = S, N = N, run_at = d$run_at[1])
    saveRDS(oldd, fd)
  }

  h <- d[d$arm == "hybrid", ]
  message(sprintf("  [theta^perp] %s%s: %+0.3f (95%% %+0.3f, %+0.3f) on the %s arm%s%s -> %s",
                  dataset, if (is.na(pair)) "" else paste0(" / ", pair),
                  h$mean, h$lower, h$upper, h$selected,
                  if (is.na(h$obs)) "" else sprintf("; measured %+0.3f, bias %+0.3f, covered %s",
                                                    h$obs, h$bias, h$covered),
                  if (is.na(h$obs_S)) "" else sprintf("; matched %+0.3f, bias %+0.3f, covered %s",
                                                      h$obs_S, h$bias_S, h$covered_S),
                  f))
  invisible(d)
}

## Read every dataset's theta^perp record back, stacked -- the companion to
## collect_scale_records(). Returns NULL rather than erroring when nothing has
## been recorded yet, so a reporting script still runs before the DA scripts have
## been re-run.
collect_perp_records <- function(dirs) {
  rows <- lapply(dirs, function(d) {
    f <- file.path(d, PERP_FILE)
    if (!file.exists(f)) return(NULL)
    readRDS(f)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    warning("no ", PERP_FILE, " found -- re-run the DA scripts to create them",
            call. = FALSE)
    return(NULL)
  }
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) { for (cn in setdiff(cols, names(r))) r[[cn]] <- NA; r[, cols, drop = FALSE] })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

## One run's theta^perp draws, resolved the way read_arm_sig() resolves a missing
## `pair`: unique match, or say which.
read_perp_draws <- function(dir, dataset, pair = NA_character_) {
  f <- file.path(dir, PERP_DRAWS_FILE)
  if (!file.exists(f)) {
    warning("no ", PERP_DRAWS_FILE, " in ", dir,
            " -- re-run that DA script to record the theta^perp draws", call. = FALSE)
    return(NULL)
  }
  all <- readRDS(f)
  if (is.na(pair)) {
    hit <- Filter(function(e) identical(e$dataset, dataset), all)
    if (length(hit) == 1) return(hit[[1]])
    if (length(hit) == 0) {
      warning("no theta^perp draws for ", dataset, " in ", f, call. = FALSE)
      return(NULL)
    }
    warning(dataset, " has ", length(hit), " recorded runs (",
            paste(vapply(hit, `[[`, "", "pair"), collapse = ", "),
            ") -- name one with `pair`", call. = FALSE)
    return(NULL)
  }
  e <- all[[arm_sig_key(dataset, pair)]]
  if (is.null(e))
    warning("no theta^perp draws for ", arm_sig_key(dataset, pair), " in ", f, call. = FALSE)
  e
}

################################################################################
## score_theta_perp: bias, RMSE and interval coverage for theta^perp, by arm.
##
## The three numbers a reviewer asks for when the claim is "we recover the
## absolute shift", computed over every recorded run that has an external
## measurement attached. Rows without one are dropped -- and counted, so a table
## built mostly from unscored runs cannot be mistaken for a strong result.
##
##   score_theta_perp(collect_perp_records(dirs))
##   score_theta_perp(rec, by = c("dataset", "arm"))   # per study
##
## RMSE over a single run is just |bias|; it is the aggregate over Jin's six
## (level, location) runs, or over the three studies, that makes it informative.
## `coverage` is the fraction of runs whose 95% credible interval contained the
## measured value, so it is only interpretable once there are several runs --
## `n` is reported next to it for exactly that reason.
################################################################################
score_theta_perp <- function(rec, by = "arm") {
  if (is.null(rec) || !nrow(rec)) return(NULL)
  n_total <- nrow(rec)
  d <- rec[is.finite(rec$obs) & is.finite(rec$mean), , drop = FALSE]
  if (!nrow(d)) {
    warning("no recorded run carries an external theta^perp measurement -- ",
            "pass `obs =` to record_theta_perp()", call. = FALSE)
    return(NULL)
  }
  ## Records written before obs_S existed have no such column at all; give them
  ## one so a mixed set of study directories still scores instead of erroring.
  if (is.null(d$obs_S)) d$obs_S <- NA_real_
  key <- do.call(paste, c(unname(d[, by, drop = FALSE]), sep = "|"))
  out <- do.call(rbind, lapply(split(d, key), function(g) {
    r  <- g[1, by, drop = FALSE]
    ## The matched half is scored over its OWN subset, not over `g`. A run with
    ## no obs_S must not drag the RMSE of the ones that have it, and `n_S`
    ## sitting next to the figures is what says how many runs it rests on.
    gs <- g[is.finite(g$obs_S), , drop = FALSE]
    data.frame(r,
               n          = nrow(g),
               bias       = mean(g$bias),
               mean_abs   = mean(g$abs_err),
               rmse       = sqrt(mean(g$bias^2)),
               coverage   = mean(g$covered),
               mean_width = mean(g$width),
               n_S        = nrow(gs),
               bias_S     = if (nrow(gs)) mean(gs$bias_S)             else NA_real_,
               mean_abs_S = if (nrow(gs)) mean(gs$abs_err_S)          else NA_real_,
               rmse_S     = if (nrow(gs)) sqrt(mean(gs$bias_S^2))     else NA_real_,
               coverage_S = if (nrow(gs)) mean(gs$covered_S)          else NA_real_,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  attr(out, "n_unscored")   <- n_total - nrow(d)
  attr(out, "n_unmatched")  <- sum(!is.finite(d$obs_S))
  out
}
