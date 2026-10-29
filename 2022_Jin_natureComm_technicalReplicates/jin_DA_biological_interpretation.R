################################################################################
## Jin et al. (2022): biological interpretation of the Sparse SSRV real-data
## analysis, and the LaTeX tables built from it.
##
## Every number that reaches the LaTeX output is read out of a saved analysis
## object at run time. Nothing is transcribed by hand, so re-running the analysis
## and then re-running this file is enough to bring the manuscript tables back
## in sync with the results.
##
## SOURCES
##   jin_results_all_levels_<DATESTAMP>.RData
##       The Sparse SSRV side. One entry per (level, location), each carrying
##         $D, $n            dimensions of that subset,
##         $tax_lookup       taxonomy for every feature, plus the ASV/OTU -> cOTU
##                           map at the ASV and OTU levels,
##         $sparse_out       whose attributes hold the FULL per-arm tables
##                           (arm_out: interval, posterior mean, p, q, sig.q for
##                           every feature, not only the called ones), the
##                           per-arm significance vectors (sig_by_arm), the
##                           selected scale model and the selection rule,
##         $ckf              the CKF run for that subset (NULL at the ASV and OTU
##                           levels, which have no measured reference). Carries
##                           `reference_taxon`: MRCKF builds its design as
##                           log(Z[,t] / W[,pstar]), so the LAST column of W has
##                           no design column and can never be selected. Which
##                           taxon that is depends on the row order after
##                           filtering, so it is looked up rather than reasoned
##                           about -- see check_ckf_reference() below, which
##                           cross-checks it against the features these tables
##                           report.
##       Reading arm_out rather than the ROOT-level real_data_arm_features.RData
##       is what lets the discordance table report an LFC and a q-value for the
##       cOTUs that only Jin et al. called: that record keeps significant
##       features only, arm_out keeps all of them. The two agree exactly on the
##       features they share (checked in check_against_arm_features() below).
##
##   Supplementary Data 9.xlsx
##       The Jin et al. side: the DESeq2 volcano source data behind
##       Supplementary Figure 18a-f. Which sheet is which location is NOT taken
##       on faith -- see jin_published() -- it is decided by matching feature
##       IDs against each location's own feature set.
##
## LOG BASE -- both sides are log2, and they are directly comparable.
##   Sparse SSRV: sparse_DA_09012026.R builds the resampled array as
##   `resam.array[, , n] <- log2(dirichlet_samples)`, so theta, its posterior
##   mean and the credible interval are all log2 fold changes.
##   Jin et al.: DESeq2 log2FoldChange, as the workbook column name says.
##   Both use VA-deficient versus VA-sufficient, so signs match too. An earlier
##   draft of the supplement described the Sparse SSRV estimates as natural-log;
##   that was wrong, and the captions emitted below say log2 for both.
##   The one caveat left is the estimand, not the units: Sparse SSRV targets the
##   ABSOLUTE (scale-corrected) LFC while DESeq2 targets the compositional one,
##   which at the proximal location differ by the measured load shift of
##   0.86 log2 units. verify_log_base() prints the evidence on each run.
################################################################################

################################################################################
## PATHS -- nothing below this block hardcodes an absolute path.
##
## ROOT is auto-detected rather than toggled by hand, so this same file runs
## unedited on the cluster and on the laptop. The first candidate that actually
## exists wins; to point somewhere else, just assign ROOT yourself afterwards.
################################################################################
ROOT_CANDIDATES <- local({  # manuscript_code: the folder holding .manuscript_code_root (README.md)
  p <- normalizePath(getwd(), mustWork = FALSE)
  while (!file.exists(file.path(p, ".manuscript_code_root"))) {
    if (dirname(p) == p) stop("run from inside the manuscript_code folder (see README.md)")
    p <- dirname(p)
  }
  p
})
ROOT <- ROOT_CANDIDATES[dir.exists(ROOT_CANDIDATES)][1]
if (is.na(ROOT)) {
  stop("None of the candidate project roots exist:\n  ",
       paste(ROOT_CANDIDATES, collapse = "\n  "))
}
message("ROOT = ", ROOT)

DIR_JIN   <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates")
DATESTAMP <- "09012026"
DIR_TAB   <- file.path(DIR_JIN, "tables")     # one .tex per table lands here
dir.create(DIR_TAB, showWarnings = FALSE)

## Appended to every file this script writes. Empty for the hybrid rule, which is
## the headline result; a single-arm run is marked so it cannot overwrite it.
arm_suffix <- function(arm = ARM) if (identical(arm, "hybrid")) "" else paste0("_", arm)

## Under the hybrid rule the reported arm is chosen per subset, and since the
## rule moved to median credible-interval width the six Jin subsets no longer
## agree -- so "hybrid Sparse SSRV" in a caption is not one estimator applied six
## times. Said once, wherever a caption or a marker depends on the arm.
arm_clause <- function() {
  if (!identical(ARM, "hybrid")) return("")
  paste0(" Under the hybrid rule the arm is chosen separately for each subset, so it ",
         "is not necessarily the same one throughout; Table~\\ref{tab:jin2022_summary} ",
         "names the arm selected for each.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

################################################################################
## CONFIGURATION
################################################################################

## Which result the manuscript reports. Everything downstream keys off this.
##
##   "hybrid"     the arm the selection rule chose FOR THAT SUBSET -- so the
##                choice can differ between, say, cOTU proximal and ASV distal.
##                This is what sparse_auto() returns and what the cross-dataset
##                results tables report for every other study, which is the
##                reason it is the default here too.
##   "bootstrap"  the bootstrap arm everywhere, whatever the rule chose.
##   "laplace"    likewise for Laplace.
##
## Outputs are named after this: "hybrid" writes the plain filenames and the
## single-arm settings append "_bootstrap" / "_laplace", so switching ARM cannot
## silently overwrite a table built from a different arm.
##
## Settable from the caller -- `Rscript -e 'ARM <- "bootstrap"; source(...)'` --
## the same `if (!exists())` pattern task_common.R uses for DATESTAMP. It used to
## be an unconditional assignment, so building the _bootstrap variant of every
## table meant editing this line between the two passes and remembering to put it
## back. The default is unchanged.
if (!exists("ARM")) ARM <- "hybrid"

## Jin et al.'s own significance rule, stated in their Fig. 7 / Supplementary
## Fig. 18 legends: FDR < 0.05 AND absolute fold change > 2, the latter being
## |log2FoldChange| > 1.
JIN_FDR_MAX        <- 0.05
JIN_ABS_LOG2FC_MIN <- 1

## Human-readable location labels, and the order rows appear in.
##
## This file builds the LOCATION-STRATIFIED analysis -- the SUPPLEMENTARY half of
## the Jin results. LOC_KEYS is deliberately prox/dist and deliberately excludes
## the pooled ("both") subset that sparse_DA_09012026.R also writes: pooled is the
## PRIMARY result and is reported by `real data results.R` (see JIN_PAIR there),
## against the measured-scale gold reference, which is empty at 4 vs 4 and so
## non-empty only when the two locations are pooled.
##
## Everything downstream -- build_runs(), the T4.x tables in manuscript_tables.R,
## the cross-representation builders -- iterates LEVEL_KEYS x LOC_KEYS and looks
## subsets up BY NAME, so a pooled entry sitting in the same results object is
## simply never read. Adding "both" here would put the primary result into the
## supplementary tables as if it were a third cecal location; it is not one -- it
## is the same eight mice again.
LOC_KEYS   <- c("prox", "dist")
LOC_LABEL  <- c(prox = "Proximal", dist = "Distal")
LEVEL_KEYS <- c("cOTU", "ASV", "OTU")

## The workbook sheets holding the published volcano data. The (level, location)
## each sheet belongs to is resolved at run time, not assumed from this order.
JIN_SOURCE_BOOK   <- "Supplementary Data 9.xlsx"
JIN_SOURCE_SHEETS <- list(
  cOTU = c("Supplementary Figure 18a", "Supplementary Figure 18b"),
  ASV  = c("Supplementary Figure 18c", "Supplementary Figure 18d"),
  OTU  = c("Supplementary Figure 18e", "Supplementary Figure 18f")
)

## Documented departures from the saved mapping object. The BarBIQ mapping
## leaves Otu032 as "Not detected", but Jin et al. name Otu032 together with
## Otu029 as the OTU representation of COTU-CM-2002 (Fig. 7 and the
## accompanying text). Features listed here are re-mapped and get a footnote
## marker wherever they appear, so the override is visible in the typeset table
## rather than buried in the code.
MANUAL_COTU_MAP <- c("OTU-CM-032" = "COTU-CM-2002")

## Rounding used in the emitted tables.
DIG_LFC_LONG  <- 3    # supplementary long tables
DIG_LFC_SHORT <- 2    # main-text comparison tables
DIG_Q         <- 3    # Sparse SSRV q-values

################################################################################
## LOAD -- Sparse SSRV
################################################################################

results_all <- readRDS(file.path(DIR_JIN,
                                 paste0("jin_results_all_levels_", DATESTAMP, ".RData")))

## `pair` is the key the ROOT-level records use: the bare location at cOTU
## level, "<level>_<location>" at ASV and OTU level.
pair_key <- function(level, loc) if (identical(level, "cOTU")) loc else paste0(level, "_", loc)

## One tidy row per (level, location, feature): taxonomy, the mapped cOTU, and
## both arms' posterior mean / q-value / call. This is the single frame every
## table below is cut from.
build_ssrv <- function(res) {
  out <- list()
  for (lv in LEVEL_KEYS) for (loc in LOC_KEYS) {
    o   <- res[[lv]][[loc]]
    tl  <- o$tax_lookup
    arm <- attr(o$sparse_out, "arm_out")
    sig <- attr(o$sparse_out, "sig_by_arm")

    ## At cOTU level a feature IS its own cOTU; the saved lookup has no `cOTU`
    ## column there. Filling it in makes the mapped-cOTU logic uniform across
    ## the three levels instead of special-casing cOTU at every use.
    mapped <- if (is.null(tl$cOTU)) tl$category else tl$cOTU
    overridden <- tl$category %in% names(MANUAL_COTU_MAP)
    mapped[overridden] <- MANUAL_COTU_MAP[tl$category[overridden]]

    i <- match(tl$category, rownames(arm$bootstrap))
    stopifnot(!anyNA(i))

    out[[length(out) + 1]] <- data.frame(
      level          = lv,
      location       = loc,
      pair           = pair_key(lv, loc),
      selected_arm   = attr(o$sparse_out, "scale_model"),
      feature        = tl$category,
      cOTU           = mapped,
      cOTU_override  = overridden,
      Kingdom = tl$Kingdom, Phylum = tl$Phylum, Class = tl$Class,
      Order   = tl$Order,   Family = tl$Family, Genus = tl$Genus,
      lfc_prop       = tl$lfc_prop,
      lfc_abs        = tl$lfc_abs,
      mean_bootstrap = arm$bootstrap$mean[i],
      ## Columns 1 and 2 of an arm table are the 2.5% and 97.5% quantiles.
      ## Positional rather than by name: the names carry a literal `%`.
      lower_bootstrap = arm$bootstrap[[1]][i],
      upper_bootstrap = arm$bootstrap[[2]][i],
      q_bootstrap    = arm$bootstrap$q.val[i],
      sig_bootstrap  = sig$bootstrap[match(tl$category, names(sig$bootstrap))] %in% TRUE,
      mean_laplace   = arm$laplace$mean[i],
      lower_laplace  = arm$laplace[[1]][i],
      upper_laplace  = arm$laplace[[2]][i],
      q_laplace      = arm$laplace$q.val[i],
      sig_laplace    = sig$laplace[match(tl$category, names(sig$laplace))] %in% TRUE,
      stringsAsFactors = FALSE
    )
  }
  d <- do.call(rbind, out)
  ## The reported result, under neutral names, so ARM is the only switch to flip.
  ## Under "hybrid" the pick is PER SUBSET -- `selected_arm` carries which one --
  ## because the selection rule runs separately for every (level, location).
  pick <- function(boot, lap) {
    if (ARM == "bootstrap") return(boot)
    if (ARM == "laplace")   return(lap)
    ifelse(d$selected_arm == "bootstrap", boot, lap)
  }
  d$lfc   <- pick(d$mean_bootstrap,  d$mean_laplace)
  d$lower <- pick(d$lower_bootstrap, d$lower_laplace)
  d$upper <- pick(d$upper_bootstrap, d$upper_laplace)
  d$q     <- pick(d$q_bootstrap,     d$q_laplace)
  d$sig   <- pick(d$sig_bootstrap,   d$sig_laplace)
  rownames(d) <- NULL
  d
}

ssrv <- build_ssrv(results_all)

## Under "hybrid" the reported values must reproduce sparse_out, which IS the
## selected arm's table. Checked rather than assumed: an off-by-one in the
## per-subset pick would silently report the arm that was NOT chosen.
check_hybrid_matches_fit <- function(res, ssrv_tbl) {
  if (!identical(ARM, "hybrid")) return(invisible(NULL))
  worst <- 0
  for (lv in LEVEL_KEYS) for (loc in LOC_KEYS) {
    so <- res[[lv]][[loc]]$sparse_out
    d  <- ssrv_tbl[ssrv_tbl$level == lv & ssrv_tbl$location == loc, ]
    i  <- match(d$feature, rownames(so))
    worst <- max(worst, max(abs(d$lfc - so$mean[i])), max(abs(d$q - so$q.val[i])))
    if (!identical(d$sig, so$sig.q[i] %in% TRUE))
      stop("hybrid calls disagree with sparse_out at ", lv, "/", loc)
  }
  if (worst > 1e-10) stop("hybrid estimates disagree with sparse_out by ", worst)
  message("Hybrid reporting reproduces sparse_out exactly at all six subsets")
  invisible(NULL)
}
check_hybrid_matches_fit(results_all, ssrv)

## Per-(level, location) run metadata: dimensions, which scale model was
## selected, under which rule, and how many features each arm called.
build_runs <- function(res) {
  rows <- list()
  for (lv in LEVEL_KEYS) for (loc in LOC_KEYS) {
    o   <- res[[lv]][[loc]]
    sig <- attr(o$sparse_out, "sig_by_arm")
    b   <- names(sig$bootstrap)[sig$bootstrap %in% TRUE]
    l   <- names(sig$laplace)[sig$laplace %in% TRUE]
    chosen <- attr(o$sparse_out, "scale_model")
    rows[[length(rows) + 1]] <- data.frame(
      dataset     = "Jin2022",
      level       = lv,
      location    = loc,
      pair        = pair_key(lv, loc),
      D           = o$D,
      N           = o$n,
      scale_model = chosen,
      scale_rule  = attr(o$sparse_out, "scale_rule"),
      n_selected  = if (identical(chosen, "bootstrap")) length(b) else length(l),
      n_bootstrap = length(b),
      n_laplace   = length(l),
      n_shared    = length(intersect(b, l)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

runs <- build_runs(results_all)

################################################################################
## LOAD -- Jin et al.'s published DESeq2 results
##
## The workbook labels the volcano sheets a-f without saying which of each pair
## is proximal, so the assignment is derived instead of assumed: a sheet belongs
## to the location whose feature set contains all of its IDs. The margin is not
## subtle (152/152 versus 138/152 for the cOTU pair, and similarly for the other
## two), and the run aborts rather than guessing if it ever narrows.
################################################################################

read_volcano <- function(path, sheet) {
  ## Row 1 is the figure label, row 2 is blank, row 3 is the header, and every
  ## value is stored as text.
  d <- suppressMessages(as.data.frame(read_xlsx(path, sheet = sheet, skip = 2)))
  names(d)[1:3] <- c("FDR", "log2FC", "feature")
  d$FDR    <- as.numeric(d$FDR)
  d$log2FC <- as.numeric(d$log2FC)
  d <- d[!is.na(d$feature) & !is.na(d$log2FC), c("feature", "log2FC", "FDR")]
  rownames(d) <- NULL
  d
}

jin_published <- function(ssrv_tbl, path = file.path(DIR_JIN, JIN_SOURCE_BOOK)) {
  if (!file.exists(path)) {
    stop("Jin et al. source-data workbook not found:\n  ", path,
         "\nThe published-comparison tables cannot be built without it.")
  }
  out <- list()
  for (lv in names(JIN_SOURCE_SHEETS)) {
    feats <- split(ssrv_tbl$feature[ssrv_tbl$level == lv],
                   ssrv_tbl$location[ssrv_tbl$level == lv])
    for (sh in JIN_SOURCE_SHEETS[[lv]]) {
      d     <- read_volcano(path, sh)
      ## Fraction of the sheet's IDs present in each location's feature set.
      ## The right location contains all of them; the wrong one is short by the
      ## features that location filtered out. The winner has to match
      ## essentially everything AND clearly beat the runner-up -- anything less
      ## means the workbook layout changed and the assignment is a guess.
      share  <- vapply(LOC_KEYS, function(k) mean(d$feature %in% feats[[k]]), numeric(1))
      ranked <- sort(share, decreasing = TRUE)
      best   <- LOC_KEYS[which.max(share)]
      if (ranked[1] < 0.99 || ranked[1] - ranked[2] < 0.02) {
        stop("cannot assign sheet '", sh, "' (", lv, ") to a location; ",
             "ID overlap = ", paste(sprintf("%s %.3f", LOC_KEYS, share), collapse = ", "))
      }
      message(sprintf("  %-26s -> %-5s %-4s  (ID overlap %.3f vs %.3f)",
                      sh, lv, best, max(share), min(share)))
      d$level    <- lv
      d$location <- best
      d$sheet    <- sh
      d$sig      <- d$FDR < JIN_FDR_MAX & abs(d$log2FC) > JIN_ABS_LOG2FC_MIN
      out[[length(out) + 1]] <- d
    }
  }
  d <- do.call(rbind, out)
  ## One sheet per (level, location): if two sheets of a pair were sent to the
  ## same location the assignment is wrong, however confident each looked.
  assigned <- unique(d[, c("level", "location", "sheet")])
  if (anyDuplicated(paste(assigned$level, assigned$location))) {
    stop("two sheets were assigned to the same (level, location):\n",
         paste(sprintf("  %s %s <- %s", assigned$level, assigned$location, assigned$sheet),
               collapse = "\n"))
  }
  rownames(d) <- NULL
  d[, c("level", "location", "feature", "log2FC", "FDR", "sig", "sheet")]
}

message("Assigning published volcano sheets to locations:")
jin <- jin_published(ssrv)

## Convenience lookups keyed on level|location|feature.
jin_key   <- with(jin,  paste(level, location, feature, sep = "|"))
jin_lfc   <- setNames(jin$log2FC, jin_key)
jin_fdr   <- setNames(jin$FDR,    jin_key)
jin_sig   <- setNames(jin$sig,    jin_key)
kof       <- function(level, location, feature) paste(level, location, feature, sep = "|")

################################################################################
## CONSISTENCY CHECKS -- run every time, so a stale object is caught here rather
## than in the typeset table.
################################################################################

## The ROOT record covers every Jin run recorded on disk; this file builds only
## the location-stratified ones (see LOC_KEYS). Since 2026-08-28 the DA script
## also writes the POOLED cOTU run under pair = "both", which `real data results.R`
## reports as the primary Jin result -- so the ROOT record legitimately carries
## rows this file has no counterpart for, and matching every row unconditionally
## turned that into a hard failure.
##
## The scope is taken from LEVEL_KEYS x LOC_KEYS, NOT from ssrv_tbl's own pairs:
## deriving it from the table would make a genuinely missing subset look
## out-of-scope and silently disarm the check, which is the one thing it exists
## to catch. In-scope rows must still match exactly; out-of-scope rows are named
## and skipped.
check_against_arm_features <- function(ssrv_tbl) {
  f <- file.path(ROOT, "real_data_arm_features.RData")
  if (!file.exists(f)) { message("real_data_arm_features.RData absent -- skipping cross-check"); return(invisible(NULL)) }
  ft <- readRDS(f)$features
  ft <- ft[ft$dataset == "Jin2022", ]

  scope <- unlist(lapply(LEVEL_KEYS, function(lv) vapply(LOC_KEYS, function(lo)
    pair_key(lv, lo), character(1))), use.names = FALSE)
  outside <- ft[!(ft$pair %in% scope), ]
  if (nrow(outside)) {
    message(sprintf("Cross-check: skipping %d row(s) from run(s) this file does not build: %s",
                    nrow(outside), paste(sort(unique(outside$pair)), collapse = ", ")))
  }
  ft <- ft[ft$pair %in% scope, ]
  if (!nrow(ft)) { message("Cross-check vs real_data_arm_features: no in-scope rows recorded"); return(invisible(NULL)) }

  i  <- match(paste(ft$pair, ft$feature), paste(ssrv_tbl$pair, ssrv_tbl$feature))
  if (anyNA(i)) {
    miss <- unique(paste(ft$pair[is.na(i)], ft$feature[is.na(i)]))
    stop("features recorded at ROOT are missing from the per-location objects: ",
         paste(head(miss, 10), collapse = ", "),
         if (length(miss) > 10) sprintf(" (and %d more)", length(miss) - 10) else "")
  }
  dm <- max(abs(ft$mean_bootstrap - ssrv_tbl$mean_bootstrap[i]))
  dq <- max(abs(ft$q_bootstrap    - ssrv_tbl$q_bootstrap[i]))
  message(sprintf("Cross-check vs real_data_arm_features: max |dmean| = %.3g, max |dq| = %.3g", dm, dq))
  if (dm > 1e-8 || dq > 1e-8) warning("ROOT-level record disagrees with the per-location objects")
  invisible(NULL)
}

check_scale_model_record <- function(runs_tbl) {
  f <- file.path(ROOT, "real_data_scale_model.RData")
  if (!file.exists(f)) { message("real_data_scale_model.RData absent -- skipping cross-check"); return(invisible(NULL)) }
  sm <- readRDS(f)
  sm <- sm[sm$dataset == "Jin2022" & !is.na(sm$pair), ]
  i  <- match(runs_tbl$pair, sm$pair)
  ok <- !is.na(i)
  bad <- runs_tbl$pair[ok][ runs_tbl$n_bootstrap[ok] != sm$n_sig_bootstrap[i[ok]] |
                            runs_tbl$n_laplace[ok]   != sm$n_sig_laplace[i[ok]]   |
                            runs_tbl$D[ok]           != sm$D[i[ok]] ]
  if (length(bad)) {
    warning("real_data_scale_model.RData disagrees with the analysis objects for: ",
            paste(bad, collapse = ", "))
  } else {
    message("Cross-check vs real_data_scale_model: counts and dimensions agree")
  }
  invisible(NULL)
}

## Prints the evidence for the log-base claim made in the captions: Jin et al.'s
## log2FoldChange against this pipeline's own log2 empirical LFC and against the
## Sparse SSRV posterior mean. If the posterior mean were natural-log, dividing
## it by ln 2 would land near the empirical log2 value; it does not.
verify_log_base <- function(ssrv_tbl) {
  d <- ssrv_tbl[ssrv_tbl$level == "cOTU" & ssrv_tbl$sig & is.finite(ssrv_tbl$lfc_prop), ]
  if (!nrow(d)) { message("log-base check: no finite-LFC significant cOTU to compare"); return(invisible(NULL)) }
  d$jin <- jin_lfc[kof(d$level, d$location, d$feature)]
  message("Log-base check (all columns log2 unless noted):")
  print(data.frame(location = d$location, feature = d$feature,
                   jin_log2        = round(d$jin, 3),
                   empirical_log2  = round(d$lfc_prop, 3),
                   ssrv_posterior  = round(d$lfc, 3),
                   if_it_were_ln   = round(d$lfc / log(2), 3),
                   row.names = NULL))
  invisible(NULL)
}

## The taxon MRCKF dropped, per subset, checked against what these tables report.
##
## CKF is not one of the methods compared here -- that comparison lives in the
## TP/FP table in `real data results.R` -- but the dropped taxon is a fact about
## the fit that cannot be recovered from its output, and the question it raises
## is answered with THIS script's objects: was the untestable taxon one of the
## features we go on to report as a real finding? If so, CKF's failure to call it
## says nothing about CKF, and the TP/FP table understates it by at least one.
##
## Writes jin_ckf_reference_<DATESTAMP>.csv so the answer sits next to the
## tables rather than only in a console message.
check_ckf_reference <- function(ssrv_tbl, res) {
  rows <- list()
  for (lv in LEVEL_KEYS) for (loc in LOC_KEYS) {
    ckf <- res[[lv]][[loc]]$ckf
    if (is.null(ckf)) next                    # not run at this level
    f  <- ckf$reference_taxon
    i  <- which(ssrv_tbl$level == lv & ssrv_tbl$location == loc &
                ssrv_tbl$feature == f)
    ## The dropped taxon should be a feature of this subset; if it is not, the
    ## CKF record and the analysis object came from different runs.
    if (!length(i)) {
      warning("CKF reference taxon '", f, "' for ", lv, "/", loc,
              " is not a feature of that subset -- the ckf record and the ",
              "analysis object are out of sync", call. = FALSE)
      next
    }
    rows[[length(rows) + 1]] <- data.frame(
      level = lv, location = loc, pair = pair_key(lv, loc),
      dropped_feature = f,
      mapped_cOTU     = ssrv_tbl$cOTU[i],
      Phylum = ssrv_tbl$Phylum[i], Family = ssrv_tbl$Family[i],
      Genus  = ssrv_tbl$Genus[i],
      lfc = ssrv_tbl$lfc[i], q = ssrv_tbl$q[i],
      sig_in_sparse_SSRV = ssrv_tbl$sig[i],
      in_ckf_truth_set   = ckf$reference_is_truth,
      D = ckf$D, run_at = ckf$run_at,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    message("No CKF run recorded on the analysis objects -- skipping reference check")
    return(invisible(NULL))
  }
  d <- do.call(rbind, rows); rownames(d) <- NULL
  message("CKF reference taxon (dropped by MRCKF, never selectable):")
  print(d[, c("level", "location", "dropped_feature", "mapped_cOTU", "Genus",
              "sig_in_sparse_SSRV", "in_ckf_truth_set")], row.names = FALSE)

  f_out <- file.path(DIR_JIN, paste0("jin_ckf_reference_", DATESTAMP, arm_suffix(), ".csv"))
  write.csv(d, f_out, row.names = FALSE)
  message("Wrote ", f_out)

  bad <- d[d$sig_in_sparse_SSRV %in% TRUE | d$in_ckf_truth_set %in% TRUE, ]
  if (nrow(bad)) {
    warning("the CKF reference taxon is a reported finding for: ",
            paste0(bad$pair, " (", bad$dropped_feature, ")", collapse = ", "),
            " -- CKF could not have selected it, so its TP is understated. ",
            "Re-run that subset with `reference = <another taxon>`.",
            call. = FALSE, immediate. = TRUE)
  } else {
    message("  none of the dropped taxa is significant here or in a CKF truth set")
  }
  invisible(d)
}

check_against_arm_features(ssrv)
check_scale_model_record(runs)
verify_log_base(ssrv)
ckf_reference <- check_ckf_reference(ssrv, results_all)

################################################################################
## LATEX HELPERS
################################################################################

DASH <- "--"

tex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  for (ch in c("&", "%", "$", "#", "_", "{", "}")) {
    x <- gsub(ch, paste0("\\", ch), x, fixed = TRUE)
  }
  x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("^", "\\textasciicircum{}", x, fixed = TRUE)
  x
}

## Math minus, so negatives typeset as $-5.37$ rather than a hyphen.
tex_num <- function(x, digits = DIG_LFC_SHORT, na = DASH) {
  vapply(x, function(v) {
    if (is.na(v))       return(na)
    if (!is.finite(v))  return(if (v > 0) "$\\infty$" else "$-\\infty$")
    s <- formatC(v, format = "f", digits = digits)
    if (v < 0) paste0("$", s, "$") else s
  }, character(1), USE.NAMES = FALSE)
}

## Three significant digits, switching to scientific notation below 1e-3 --
## the convention the published FDR column already used.
tex_pval <- function(x, na = DASH) {
  vapply(x, function(v) {
    if (is.na(v)) return(na)
    if (v == 0)   return("$<10^{-4}$")           # below the posterior-draw resolution
    if (v < 1e-3) {
      e <- floor(log10(v)); m <- v / 10^e
      return(sprintf("$%.2f\\times10^{%d}$", m, e))
    }
    formatC(v, format = "fg", digits = 3, flag = "#")
  }, character(1), USE.NAMES = FALSE)
}

## Sparse SSRV q-values: plain decimals, since BH on Monte-Carlo p-values gives
## a coarse grid; an exact zero means the p-value fell below the draw resolution.
tex_q <- function(x, digits = DIG_Q, na = DASH) {
  vapply(x, function(v) {
    if (is.na(v))  return(na)
    if (v == 0)    return("$<10^{-4}$")
    if (v < 10^(-digits)) return(tex_pval(v))
    formatC(v, format = "f", digits = digits)
  }, character(1), USE.NAMES = FALSE)
}

tex_genus <- function(x) {
  ifelse(is.na(x) | x == "", DASH, paste0("\\textit{", tex_escape(x), "}"))
}

tex_tax <- function(x) ifelse(is.na(x) | x == "", DASH, tex_escape(x))

## Several features stacked in one cell, left-aligned, as in the existing tables.
tex_stack <- function(items) {
  if (!length(items)) return(DASH)
  if (length(items) == 1) return(items)
  paste0("\\shortstack[l]{", paste(items, collapse = " \\\\ "), "}")
}

row_tex <- function(...) paste0("        ", paste(c(...), collapse = " & "), " \\\\")

## A 95% credible interval, set in math mode so the minus signs are real minus
## signs rather than hyphens.
tex_interval <- function(lo, hi, digits = DIG_LFC_SHORT) {
  if (length(lo) != 1L || length(hi) != 1L || is.na(lo) || is.na(hi)) return("")
  sprintf("$[%s,\\, %s]$",
          formatC(lo, format = "f", digits = digits),
          formatC(hi, format = "f", digits = digits))
}

## Significance markers: * called by Jin et al., \dag called by us, and the
## trailing letter flags a feature whose cOTU mapping we overrode. \dag is
## written \dag{} so a following letter is not swallowed into the macro name.
sig_marks <- function(level, location, feature, our_sig, override = FALSE) {
  m <- character(0)
  if (isTRUE(unname(jin_sig[kof(level, location, feature)]))) m <- c(m, "*")
  if (isTRUE(our_sig))  m <- c(m, "\\dag{}")
  if (isTRUE(override)) m <- c(m, "a")
  if (!length(m)) return("")
  paste0("\\textsuperscript{", paste(m, collapse = ""), "}")
}

## A stack of entries in one cell, TOP-aligned: a nested tabular with [t] puts
## its first line on the outer row's baseline, so the cOTU in column 2 sits
## level with the first mapped feature rather than floating to the middle the
## way \shortstack would leave it.
tex_block <- function(entries, sep = " \\\\[3pt] ", align = "l") {
  if (!length(entries)) return(DASH)
  paste0("\\begin{tabular}[t]{@{}", align, "@{}}",
         paste(entries, collapse = sep),
         "\\end{tabular}")
}

## Feature IDs sort by their trailing number, so ASV-CM-106 follows ASV-CM-93.
## nat_num is the key on its own, for callers that need it inside a wider order().
nat_num <- function(x) {
  num <- suppressWarnings(as.numeric(sub("^.*?([0-9]+)$", "\\1", x)))
  ifelse(is.na(num), Inf, num)
}
nat_order <- function(x) order(sub("([0-9]+)$", "", x), nat_num(x), x)

################################################################################
## TABLE BUILDERS -- each returns a character vector of LaTeX lines.
################################################################################

## ---------------------------------------------------------------------------
## Dimensions and scale-model selection, one row per analysed subset.
## ---------------------------------------------------------------------------
tab_scale_model <- function(runs_tbl) {
  d <- runs_tbl[order(match(runs_tbl$level, LEVEL_KEYS),
                      match(runs_tbl$location, LOC_KEYS)), ]
  body <- mapply(function(...) row_tex(...),
                 tex_escape(d$dataset), tex_escape(d$pair),
                 sprintf("$%d \\times %d$", d$D, d$N),
                 tex_escape(d$scale_model), tex_escape(d$scale_rule),
                 d$n_selected, d$n_bootstrap, d$n_laplace, d$n_shared,
                 USE.NAMES = FALSE)
  c("\\begin{table}[ht]",
    "\\centering",
    "\\caption{Dimensions of each Jin et al. subset and the number of features called by each scale model. $D$ is the number of features retained after filtering and $N$ the number of samples. The scale model column reports which arm the selection rule chose for that subset; the last four columns count the features called by the selected arm, by the bootstrap arm, by the Laplace arm, and by both arms.}",
    "\\label{tab:jin2022_summary}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lllcccccc}",
    "\\toprule",
    "Dataset & Pair & Dimension & Scale Model & Scale Rule &",
    "$n_{\\text{sig, selected}}$ &",
    "$n_{\\text{sig, bootstrap}}$ &",
    "$n_{\\text{sig, Laplace}}$ &",
    "$n_{\\text{shared}}$ \\\\",
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\end{table}")
}

## ---------------------------------------------------------------------------
## Main-text Table 1: discovery counts, published versus ours.
## ---------------------------------------------------------------------------
tab_discovery_counts <- function(ssrv_tbl, jin_tbl) {
  n_ssrv <- function(lv, loc) sum(ssrv_tbl$sig & ssrv_tbl$level == lv & ssrv_tbl$location == loc)
  n_jin  <- function(lv, loc) sum(jin_tbl$sig  & jin_tbl$level  == lv & jin_tbl$location  == loc)
  body <- vapply(LEVEL_KEYS, function(lv) {
    sp <- n_ssrv(lv, "prox"); sd <- n_ssrv(lv, "dist")
    row_tex(sprintf("%-4s", lv), n_jin(lv, "prox"), n_jin(lv, "dist"), sp, sd,
            if (sd > 0) formatC(sp / sd, format = "f", digits = 2) else DASH)
  }, character(1), USE.NAMES = FALSE)
  arm_lab <- paste0(toupper(substring(ARM, 1, 1)), substring(ARM, 2))
  c("\\begin{table}[!tbp]",
    "    \\centering",
    "    \\small",
    "    \\begin{tabular}{@{}lrrrrr@{}}",
    "        \\toprule",
    paste0("        & \\multicolumn{2}{c}{Jin et al.} & \\multicolumn{3}{c}{", arm_lab, " Sparse SSRV} \\\\"),
    "        \\cmidrule(lr){2-3}\\cmidrule(l){4-6}",
    "        Feature representation & Proximal & Distal & Proximal & Distal & Proximal/distal \\\\",
    "        \\midrule",
    body,
    "        \\bottomrule",
    "    \\end{tabular}",
    paste0("    \\caption{Numbers of significant features at the proximal and distal cecal locations reported by \\citet{jin2022quantsinglecell} and identified by ",
           tolower(arm_lab), " Sparse SSRV. The final column reports the proximal-to-distal ratio for ",
           tolower(arm_lab), " Sparse SSRV. Published counts are recomputed from the Supplementary Figure 18 source data under the rule used in that paper, FDR $<",
           JIN_FDR_MAX, "$ and absolute fold change $>", 2^JIN_ABS_LOG2FC_MIN, "$.",
           arm_clause(), "}"),
    "    \\label{tab:jin_discovery_counts}",
    "\\end{table}")
}

## ---------------------------------------------------------------------------
## Main-text Table 2: each significant cOTU beside the published estimate and
## beside whichever ASVs and OTUs mapping to it were themselves called.
## ---------------------------------------------------------------------------
tab_cross_representation <- function(ssrv_tbl, jin_tbl) {
  used_override <- FALSE
  rows <- character(0)
  for (loc in LOC_KEYS) {
    sig_c <- ssrv_tbl[ssrv_tbl$sig & ssrv_tbl$level == "cOTU" & ssrv_tbl$location == loc, ]
    if (!nrow(sig_c)) next
    ## Strongest evidence first: q-value, then magnitude.
    sig_c <- sig_c[order(sig_c$q, -abs(sig_c$lfc)), ]
    if (length(rows)) rows <- c(rows, "        \\addlinespace")
    for (i in seq_len(nrow(sig_c))) {
      f  <- sig_c$feature[i]
      k  <- kof("cOTU", loc, f)
      pub <- if (k %in% names(jin_lfc)) {
        paste0(tex_num(jin_lfc[[k]], DIG_LFC_SHORT),
               if (isTRUE(jin_sig[[k]])) "\\textsuperscript{*}" else "")
      } else DASH
      cell <- function(lv) {
        m <- ssrv_tbl[ssrv_tbl$sig & ssrv_tbl$level == lv &
                      ssrv_tbl$location == loc & ssrv_tbl$cOTU == f, ]
        if (!nrow(m)) return(DASH)
        m <- m[nat_order(m$feature), ]
        used_override <<- used_override || any(m$cOTU_override)
        tex_stack(paste0(tex_escape(m$feature),
                         ifelse(m$cOTU_override, "\\textsuperscript{a}", ""),
                         " (", tex_num(m$lfc, DIG_LFC_SHORT), ")"))
      }
      rows <- c(rows, row_tex(LOC_LABEL[[loc]], tex_escape(f), pub,
                              tex_num(sig_c$lfc[i], DIG_LFC_SHORT),
                              cell("ASV"), cell("OTU")))
    }
  }
  note <- paste0("\\footnotesize \\textsuperscript{*}Significant in \\citet{jin2022quantsinglecell}, defined by FDR $<",
                 JIN_FDR_MAX, "$ and absolute fold change $>", 2^JIN_ABS_LOG2FC_MIN, "$.")
  if (used_override) {
    ov <- names(MANUAL_COTU_MAP)
    note <- paste0(note, " \\textsuperscript{a}The saved mapping object labels ",
                   paste(tex_escape(ov), collapse = " and "),
                   " as ``Not detected.'' However, \\citet{jin2022quantsinglecell} explicitly identify ",
                   paste(tex_escape(ov), collapse = " and "),
                   " as an OTU representation of ",
                   paste(tex_escape(unname(MANUAL_COTU_MAP)), collapse = " and "), ".")
  }
  arm_lab <- paste0(toupper(substring(ARM, 1, 1)), substring(ARM, 2))
  c("\\begin{table}[!tbp]",
    "    \\centering",
    "    \\footnotesize",
    "    \\setlength{\\tabcolsep}{3pt}",
    "    \\begin{tabular}{@{}llrrll@{}}",
    "        \\toprule",
    "        Location & Significant cOTU & \\shortstack{Jin et al.\\\\$\\log_2$ LFC} & \\shortstack{Sparse SSRV\\\\$\\log_2$ LFC} & \\shortstack[l]{Significant mapped ASV(s)\\\\($\\log_2$ LFC)} & \\shortstack[l]{Significant mapped OTU(s)\\\\($\\log_2$ LFC)} \\\\",
    "        \\midrule",
    rows,
    "        \\bottomrule",
    "    \\end{tabular}",
    paste0("    \\caption{Comparison with the cOTU results reported by \\citet{jin2022quantsinglecell} and agreement among significant mapped feature representations. Published cOTU estimates are DESeq2 $\\log_2$ fold changes; Sparse SSRV estimates are ",
           tolower(arm_lab),
           " posterior means, also on the $\\log_2$ scale. All estimates use the vitamin A-deficient versus vitamin A-sufficient contrast, so signs and magnitudes are on a common footing. The two nevertheless estimate different quantities: Sparse SSRV targets the absolute fold change, which incorporates the measured difference in total microbial load between the diet groups, whereas DESeq2 targets the compositional fold change. A dash indicates that no significant mapped feature was identified for that representation. cOTUs are ordered by q-value and then by decreasing magnitude within each location.}"),
    "    \\label{tab:jin_cross_representation}",
    "    \\vspace{2pt}",
    "    \\begin{minipage}{0.96\\textwidth}",
    paste0("        ", note),
    "    \\end{minipage}",
    "\\end{table}")
}

## ---------------------------------------------------------------------------
## The comprehensive cross-representation tables. One builder, two scopes.
##
## Ten columns in three groups of three (feature, their LFC, our LFC) behind the
## location. A cOTU with several mapped ASVs or OTUs spans several physical
## rows, with the cOTU's own four cells filled on the first of them and blank
## afterwards, so the cOTU always heads its block. Stacking the mapped features
## inside a single cell instead would not work: the feature name and the
## published LFC are one line while our estimate is two, and nested stacks of
## unequal height stop lining up across the group.
##
##   scope = "cotu"     Rows are the UNION of the cOTUs the two analyses called,
##                      Jin et al.'s first, then ours -- so a cOTU we missed is
##                      visible rather than absent. Every mapped ASV and OTU is
##                      listed whether or not it was called, because a mapped
##                      feature nobody called is evidence about the cOTU too.
##                      The cOTU carries its own significance markers.
##
##   scope = "feature"  The selection runs the other way: every ASV and OTU
##                      called by either analysis is listed, EXCEPT those whose
##                      mapped cOTU already heads a block in the "cotu" table --
##                      the two are complementary rather than overlapping, so a
##                      called feature appears in exactly one of them. The cOTU
##                      is context here, not a finding, so it carries NO marker.
##                      Most of the cOTUs reached this way were called by
##                      neither analysis; some were called by one or both, and
##                      the table does not say which -- see the note printed by
##                      report_marked_cotus() after the tables are built.
##
## Unmapped features stand alone in this scope. "Not detected" is not a cOTU,
## and collapsing every unmapped ASV and OTU into one block would claim they
## share one; each gets its own block with the cOTU columns dashed.
## ---------------------------------------------------------------------------
## The cOTUs the cOTU-driven table heads a block with, at one location. Defined
## once because the feature-driven table has to exclude exactly this set, and a
## second copy of the rule is how the two tables start double-counting.
cotu_scope_keys <- function(ssrv_tbl, jin_tbl, loc) {
  ours   <- ssrv_tbl[ssrv_tbl$level == "cOTU" & ssrv_tbl$location == loc, ]
  theirs <- jin_tbl[jin_tbl$level  == "cOTU" & jin_tbl$location  == loc, ]
  union(theirs$feature[theirs$sig], ours$feature[ours$sig])
}

tab_cross_representation_full <- function(ssrv_tbl, jin_tbl, scope = c("cotu", "feature")) {
  scope <- match.arg(scope)
  used_override <- FALSE
  missing_pub   <- FALSE
  unfitted_cotu <- character(0)
  N_GROUP <- 3L                       # feature, Jin et al. LFC, our LFC

  ## Our estimate over its credible interval, top-aligned so the first line
  ## sits on the row baseline with the single-line cells beside it.
  ##
  ## Centred, not left-aligned, and the column that holds it is `c` rather than
  ## `r`. The cell is as wide as its interval and intervals differ in width, so
  ## a flush-left estimate inside a right-aligned column lands in a different
  ## place on every row. Centring both the cell and the column puts every point
  ## estimate on the column's centre line.
  est_cell <- function(lfc, lo, hi) {
    tex_block(paste0(tex_num(lfc, DIG_LFC_SHORT), " \\\\ ",
                     "{\\scriptsize ", tex_interval(lo, hi), "}"),
              align = "c")
  }

  ## `lv` may be one level for the whole vector (mapped_of) or one per feature
  ## (location_groups); rep_len covers both. Getting this wrong is silent: a
  ## vector `lv` makes kof() return a vector, isTRUE() is then always FALSE, and
  ## the published calls drop out of the selection without any error.
  called_by_either <- function(lv, loc, feats, our_sig) {
    lv <- rep_len(lv, length(feats))
    our_sig | vapply(seq_along(feats), function(i)
      isTRUE(unname(jin_sig[kof(lv[i], loc, feats[i])])), logical(1))
  }

  ## Block key: the mapped cOTU, or the feature itself when unmapped.
  gkey <- function(d) ifelse(d$cOTU == "Not detected",
                             paste0("ND|", d$level, "|", d$feature), d$cOTU)

  ## Every feature of one level in one block, as a three-column frame, or NULL.
  mapped_of <- function(lv, loc, key) {
    m <- ssrv_tbl[ssrv_tbl$level == lv & ssrv_tbl$location == loc, ]
    m <- m[gkey(m) == key, ]
    if (scope == "feature" && nrow(m)) m <- m[called_by_either(lv, loc, m$feature, m$sig), ]
    if (!nrow(m)) return(NULL)
    m <- m[nat_order(m$feature), ]
    used_override <<- used_override || any(m$cOTU_override)
    data.frame(
      name = paste0(tex_escape(m$feature),
                    vapply(seq_len(nrow(m)), function(i)
                      sig_marks(lv, loc, m$feature[i], m$sig[i], m$cOTU_override[i]),
                      character(1))),
      ## NA -- the feature was dropped by DESeq2's independent filtering and is
      ## absent from the published volcano data -- prints as a dash.
      jin  = { v <- unname(jin_lfc[kof(lv, loc, m$feature)])
               missing_pub <<- missing_pub || anyNA(v)
               tex_num(v, DIG_LFC_SHORT) },
      ours = vapply(seq_len(nrow(m)), function(i)
               est_cell(m$lfc[i], m$lower[i], m$upper[i]), character(1)),
      stringsAsFactors = FALSE)
  }

  ## One group's three cells on physical row `r` of a block.
  group_cells <- function(d, r) {
    n <- if (is.null(d)) 0L else nrow(d)
    if (n == 0L) {
      if (r == 1L) paste0("\\multicolumn{", N_GROUP, "}{c}{", DASH, "}")
      else rep("", N_GROUP)
    } else if (r <= n) {
      c(d$name[r], d$jin[r], d$ours[r])
    } else rep("", N_GROUP)
  }

  ## All physical rows for one block.
  block_rows <- function(loc, key) {
    a <- mapped_of("ASV", loc, key)
    o <- mapped_of("OTU", loc, key)
    nb <- max(1L, if (is.null(a)) 0L else nrow(a), if (is.null(o)) 0L else nrow(o))
    unmapped <- startsWith(key, "ND|")
    ours <- if (unmapped) NULL else
      ssrv_tbl[ssrv_tbl$level == "cOTU" & ssrv_tbl$location == loc &
               ssrv_tbl$feature == key, ]
    ## The ASV/OTU mapping points into the full BarBIQ cOTU catalogue, which is
    ## wider than the cOTU table analysed at this location: prevalence filtering
    ## drops some of them. Such a cOTU is named but has no estimate on either
    ## side, so its two numeric cells are dashed.
    no_cotu_fit <- !unmapped && (is.null(ours) || !nrow(ours))
    if (no_cotu_fit) unfitted_cotu <<- c(unfitted_cotu, paste(loc, key))
    vapply(seq_len(nb), function(r) {
      head_cells <- if (r == 1L) {
        if (unmapped) c(LOC_LABEL[[loc]], "Not detected", DASH, DASH)
        else if (no_cotu_fit) c(LOC_LABEL[[loc]], tex_escape(key), DASH, DASH)
        else c(LOC_LABEL[[loc]],
               paste0(tex_escape(key),
                      if (scope == "cotu") sig_marks("cOTU", loc, key, ours$sig) else ""),
               tex_num(unname(jin_lfc[kof("cOTU", loc, key)]), DIG_LFC_SHORT),
               est_cell(ours$lfc, ours$lower, ours$upper))
      } else rep("", 4)
      row_tex(c(head_cells, group_cells(a, r), group_cells(o, r)))
    }, character(1))
  }

  ## The blocks of one location, in order, as a list of key vectors -- one
  ## element per group that gets its own 6pt separator.
  location_groups <- function(loc) {
    if (scope == "cotu") {
      ours   <- ssrv_tbl[ssrv_tbl$level == "cOTU" & ssrv_tbl$location == loc, ]
      theirs <- jin_tbl[jin_tbl$level  == "cOTU" & jin_tbl$location  == loc, ]
      st <- theirs$feature[theirs$sig]
      so <- setdiff(ours$feature[ours$sig], st)
      list(st[nat_order(st)], so[nat_order(so)])
    } else {
      sel <- ssrv_tbl[ssrv_tbl$level %in% c("ASV", "OTU") & ssrv_tbl$location == loc, ]
      sel <- sel[called_by_either(sel$level, loc, sel$feature, sel$sig), ]
      ## Drop what the cOTU-driven table already shows, so the two partition the
      ## called features rather than repeating a third of them.
      sel <- sel[!(sel$cOTU %in% cotu_scope_keys(ssrv_tbl, jin_tbl, loc)), ]
      if (!nrow(sel)) return(list())
      k  <- unique(gkey(sel))
      nd <- startsWith(k, "ND|")
      ## Mapped blocks by cOTU identifier, then the unmapped ones by level and
      ## identifier, so the "Not detected" rows collect at the end.
      lv <- sub("^ND\\|([^|]+)\\|.*$", "\\1", k[nd])
      ft <- sub("^ND\\|[^|]+\\|", "", k[nd])
      list(sort(k[!nd]), k[nd][order(match(lv, LEVEL_KEYS), nat_num(ft), ft)])
    }
  }

  rows <- character(0)
  for (loc in LOC_KEYS) {
    groups <- location_groups(loc)
    if (!length(unlist(groups))) next
    if (length(rows)) rows <- c(rows, "        \\addlinespace[6pt]")
    for (bi in seq_along(groups)) {
      keys <- groups[[bi]]
      if (!length(keys)) next
      if (bi > 1) rows <- c(rows, "        \\addlinespace[6pt]")
      ## 2pt between blocks inside a group; the 6pt separators above already
      ## space the group and location boundaries.
      for (k in seq_along(keys)) {
        if (k > 1) rows <- c(rows, "        \\addlinespace[2pt]")
        rows <- c(rows, block_rows(loc, keys[k]))
      }
    }
  }

  note <- paste0("\\scriptsize \\textsuperscript{*}Significant in \\citet{jin2022quantsinglecell} ",
                 "(FDR $<", JIN_FDR_MAX, "$, absolute fold change $>", 2^JIN_ABS_LOG2FC_MIN,
                 "$). \\textsuperscript{\\dag}Significant in ", tolower(ARM),
                 " Sparse SSRV", arm_clause(),
                 " Significance there means a Benjamini-Hochberg adjusted posterior predictive p-value ",
                 "below 0.05 rather than a credible interval excluding zero -- an interval here may ",
                 "exclude zero for a feature carrying no dagger.",
                 if (scope == "cotu")
                   " A feature with neither marker maps to the cOTU but was called by neither analysis."
                 else
                   paste0(" Markers are carried by the ASVs and OTUs only. Every cOTU either analysis ",
                          "called heads a block in the preceding table, so no cOTU appearing here was ",
                          "called in its own right; the cOTU column is context, not a finding, and is ",
                          "left unmarked."),
                 if (length(unfitted_cotu))
                   paste0(" A cOTU whose own two columns are dashed was not retained in that ",
                          "location's cOTU-level analysis, so no estimate exists for it on ",
                          "either side; the ASVs and OTUs mapping to it were still tested.")
                 else "",
                 if (missing_pub)
                   paste0(" A dash in a published-LFC column marks a feature absent from the ",
                          "corresponding volcano source data.")
                 else "")
  if (used_override) {
    ov <- names(MANUAL_COTU_MAP)
    note <- paste0(note, " \\textsuperscript{a}The saved mapping object labels ",
                   paste(tex_escape(ov), collapse = " and "),
                   " as ``Not detected''; it is placed here because \\citet{jin2022quantsinglecell} ",
                   "identify it as an OTU representation of ",
                   paste(tex_escape(unname(MANUAL_COTU_MAP)), collapse = " and "), ".")
  }

  arm_lab <- paste0(toupper(substring(ARM, 1, 1)), substring(ARM, 2))
  caption <- if (scope == "cotu") {
    paste0("Every cOTU called by either analysis, with all ASVs and OTUs mapping to it, ",
           "and both analyses' estimates for all three feature representations. Within each location the ",
           "cOTUs reported by \\citet{jin2022quantsinglecell} are listed first, then those identified only by ",
           tolower(arm_lab), " Sparse SSRV; both groups are in identifier order. A cOTU with more than one ",
           "mapped ASV or OTU spans several rows and heads its own block. Mapped features are listed whether ",
           "or not either analysis called them, so agreement and disagreement across feature representations ",
           "can be read directly. A dash spanning a group indicates that no feature of that representation ",
           "maps to the cOTU.")
  } else {
    paste0("ASVs and OTUs called by either analysis whose mapped cOTU is not itself listed in ",
           "Table~\\ref{tab:jin_cross_representation_full}, grouped by that cOTU, with both analyses' ",
           "estimates for all three feature representations. Unlike the cOTU-driven table, membership ",
           "here is decided at the ASV and OTU level: a feature appears if ",
           "\\citet{jin2022quantsinglecell} called it, if ", tolower(arm_lab),
           " Sparse SSRV called it, or both. Features whose mapped cOTU is itself listed in ",
           "Table~\\ref{tab:jin_cross_representation_full} are shown there instead, so the two tables ",
           "partition the called features rather than overlapping. ",
           "The cOTU columns are therefore context and carry no significance markers. Blocks run in cOTU ",
           "identifier order within each location, followed by features the mapping leaves unassigned, ",
           "which stand alone because ``Not detected'' is not a shared cOTU. A dash spanning a group ",
           "indicates that no feature of that representation was called for that cOTU.")
  }
  caption <- paste0(caption,
    " Published estimates are DESeq2 $\\log_2$ fold changes, taken from the cOTU, ASV and OTU volcano ",
    "source data of Supplementary Figure 18; Sparse SSRV estimates are ", tolower(arm_lab),
    " posterior means with 95\\% credible intervals beneath, also on the $\\log_2$ scale. All estimates ",
    "use the vitamin A-deficient versus vitamin A-sufficient contrast. The two nevertheless estimate ",
    "different quantities: Sparse SSRV targets the absolute fold change, which incorporates the measured ",
    "difference in total microbial load between the diet groups, whereas DESeq2 targets the compositional ",
    "fold change.")

  ## Repeated three times, once per feature group.
  hdr_pair <- paste0("\\shortstack{Jin et al.\\\\$\\log_2$ LFC} & ",
                     "\\shortstack{Sparse SSRV\\\\$\\log_2$ LFC}")
  head_block <- c(
    "    \\toprule",
    "    & \\multicolumn{3}{c}{cOTU} & \\multicolumn{3}{c}{Mapped ASV} & \\multicolumn{3}{c}{Mapped OTU} \\\\",
    "    \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(l){8-10}",
    paste0("    Location & Feature & ", hdr_pair,
           " & Feature & ", hdr_pair, " & Feature & ", hdr_pair, " \\\\"),
    "    \\midrule")
  lab <- if (scope == "cotu") "tab:jin_cross_representation_full"
         else                 "tab:jin_cross_representation_sig_features"

  ## A longtable rather than a float: ten columns two lines deep overflowed a
  ## single landscape page, and the row count grows with every feature either
  ## analysis adds. No \resizebox -- longtable cannot be scaled, so the column
  ## widths have to fit the landscape text block on their own, which they do.
  c("\\clearpage",
    "\\begin{landscape}",
    "{\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\setlength{\\LTcapwidth}{\\linewidth}",
    "\\begin{longtable}{@{}llcclcclcc@{}}",
    paste0("    \\caption{", caption, "}"),
    paste0("    \\label{", lab, "}\\\\"),
    head_block,
    "    \\endfirsthead",
    "    \\multicolumn{10}{l}{\\tablename~\\thetable\\ continued}\\\\",
    head_block,
    "    \\endhead",
    "    \\midrule",
    "    \\multicolumn{10}{r}{Continued on next page}\\\\",
    "    \\endfoot",
    "    \\bottomrule",
    "    \\endlastfoot",
    rows,
    "\\end{longtable}",
    ## A plain paragraph, not a minipage: an unbreakable \linewidth-wide box
    ## placed after a longtable gets pushed onto a page of its own, orphaning
    ## the symbol key from the table it explains.
    "\\vspace{4pt}\\noindent",
    paste0("    ", note),
    "\\par",
    "}",
    "\\end{landscape}")
}

## The ASVs and OTUs the feature-scoped table must contain: every one called by
## either analysis. Computed independently of the table builder so it can be
## used to check it.
sig_features_either <- function(ssrv_tbl) {
  d <- ssrv_tbl[ssrv_tbl$level %in% c("ASV", "OTU"), ]
  keep <- d$sig | vapply(seq_len(nrow(d)), function(i)
    isTRUE(unname(jin_sig[kof(d$level[i], d$location[i], d$feature[i])])), logical(1))
  d[keep, c("level", "location", "feature")]
}

## Does the emitted table hold exactly those features, once each? The selection
## bug this guards against is silent -- blocks reachable only through a feature
## Jin et al. called simply vanished -- so it is checked rather than trusted.
check_feature_table <- function(lines, ssrv_tbl, exclude_cotu_scope = FALSE,
                                jin_tbl = NULL) {
  want <- sig_features_either(ssrv_tbl)
  if (isTRUE(exclude_cotu_scope)) {
    keep <- vapply(seq_len(nrow(want)), function(i) {
      cot <- ssrv_tbl$cOTU[ssrv_tbl$level == want$level[i] &
                           ssrv_tbl$location == want$location[i] &
                           ssrv_tbl$feature == want$feature[i]][1]
      !(cot %in% cotu_scope_keys(ssrv_tbl, jin_tbl, want$location[i]))
    }, logical(1))
    want <- want[keep, , drop = FALSE]
  }
  body <- grep("^        .+ & ", lines, value = TRUE)
  ## The lookbehind matters: without it "OTU-CM-..." also matches inside every
  ## "COTU-CM-..." in the first column.
  got  <- unlist(regmatches(body, gregexpr("(?<![A-Za-z])(ASV|OTU)-CM-[0-9]+",
                                           body, perl = TRUE)))
  miss  <- setdiff(unique(want$feature), got)
  extra <- setdiff(got, unique(want$feature))
  if (length(got) != nrow(want) || length(miss) || length(extra)) {
    stop("feature-scoped table is not the expected set of ASVs/OTUs: ",
         length(got), " cells emitted for ", nrow(want), " expected",
         if (length(miss))  paste0("; missing: ",    paste(miss,  collapse = ", ")) else "",
         if (length(extra)) paste0("; unexpected: ", paste(extra, collapse = ", ")) else "")
  }
  message("Feature-scoped table: ", nrow(want),
          " ASVs/OTUs called by either analysis, all present")
  invisible(NULL)
}

## Since the feature-scoped table excludes every cOTU the cOTU-scoped table
## heads, none of the cOTUs it reaches should have been called by either
## analysis. That is an invariant of the partition, so it is checked rather than
## assumed: anything reported here means the exclusion did not take.
report_marked_cotus <- function(ssrv_tbl, jin_tbl) {
  rows <- list()
  for (loc in LOC_KEYS) {
    sel <- ssrv_tbl[ssrv_tbl$level %in% c("ASV", "OTU") & ssrv_tbl$location == loc, ]
    keep <- sel$sig | vapply(seq_len(nrow(sel)), function(i)
      isTRUE(unname(jin_sig[kof(sel$level[i], loc, sel$feature[i])])), logical(1))
    sel <- sel[keep, , drop = FALSE]
    sel <- sel[!(sel$cOTU %in% cotu_scope_keys(ssrv_tbl, jin_tbl, loc)), , drop = FALSE]
    keep <- rep(TRUE, nrow(sel))
    cot <- setdiff(unique(sel$cOTU[keep]), "Not detected")
    ours   <- ssrv_tbl[ssrv_tbl$level == "cOTU" & ssrv_tbl$location == loc, ]
    theirs <- jin_tbl[jin_tbl$level  == "cOTU" & jin_tbl$location  == loc, ]
    for (f in cot[nat_order(cot)]) {
      in_ours   <- isTRUE(ours$sig[match(f, ours$feature)])
      in_theirs <- isTRUE(theirs$sig[match(f, theirs$feature)])
      if (in_ours || in_theirs)
        rows[[length(rows) + 1]] <- data.frame(
          location = loc, cOTU = f, jin_et_al = in_theirs, sparse_SSRV = in_ours)
    }
  }
  if (!length(rows)) {
    message("Feature-scoped table: none of its cOTUs was called by either analysis, as the ",
            "partition with the cOTU-scoped table requires")
    return(invisible(NULL))
  }
  d <- do.call(rbind, rows)
  warning("the feature-scoped table still reaches ", nrow(d), " cOTU(s) that were called, ",
          "so it overlaps the cOTU-scoped table instead of partitioning with it:\n",
          paste(sprintf("  %s %s", d$location, d$cOTU), collapse = "\n"),
          call. = FALSE, immediate. = TRUE)
  print(d, row.names = FALSE)
  invisible(d)
}

## ---------------------------------------------------------------------------
## Supplementary: every feature the reported arm called, with its LFC.
## ---------------------------------------------------------------------------
tab_all_significant <- function(ssrv_tbl) {
  body <- character(0)
  first <- TRUE
  for (loc in LOC_KEYS) for (lv in LEVEL_KEYS) {
    d <- ssrv_tbl[ssrv_tbl$sig & ssrv_tbl$level == lv & ssrv_tbl$location == loc, ]
    if (!nrow(d)) next
    d <- d[nat_order(d$feature), ]
    if (!first) body <- c(body, "    \\addlinespace[2pt]")
    first <- FALSE
    body <- c(body, paste0("    ", LOC_LABEL[[loc]], " & ", lv, " & ",
                           tex_escape(d$feature), " & ",
                           tex_num(d$lfc, DIG_LFC_LONG), " \\\\"))
  }
  arm_lab <- paste0(toupper(substring(ARM, 1, 1)), substring(ARM, 2))
  hdr <- c("    \\toprule",
           "    Location & Feature representation & Feature & Estimated LFC \\\\",
           "    \\midrule")
  c("{\\small",
    "\\setlength{\\LTcapwidth}{0.95\\textwidth}",
    "\\begin{longtable}{@{}lllr@{}}",
    paste0("    \\caption{All features identified as significant by ", tolower(arm_lab),
           " Sparse SSRV in the Jin et al. data. Estimated LFCs are posterior means for the vitamin A-deficient versus vitamin A-sufficient contrast on the $\\log_2$ scale. Features are ordered by location, representation, and feature identifier.}"),
    "    \\label{tab:supp_jin_all_significant}\\\\",
    hdr,
    "    \\endfirsthead",
    "    \\multicolumn{4}{l}{\\tablename~\\thetable\\ continued}\\\\",
    hdr,
    "    \\endhead",
    "    \\midrule",
    "    \\multicolumn{4}{r}{Continued on next page}\\\\",
    "    \\endfoot",
    "    \\bottomrule",
    "    \\endlastfoot",
    body,
    "\\end{longtable}",
    "}")
}

## ---------------------------------------------------------------------------
## Supplementary: significant ASVs (or OTUs) whose mapped cOTU was NOT itself
## called at the same location -- the discoveries the cOTU analysis missed.
## ---------------------------------------------------------------------------
tab_without_sig_cotu <- function(ssrv_tbl, lv) {
  rows <- character(0)
  excluded <- character(0)
  for (loc in LOC_KEYS) {
    sig_c <- ssrv_tbl$feature[ssrv_tbl$sig & ssrv_tbl$level == "cOTU" & ssrv_tbl$location == loc]
    d <- ssrv_tbl[ssrv_tbl$sig & ssrv_tbl$level == lv & ssrv_tbl$location == loc, ]
    if (!nrow(d)) next
    drop <- d$cOTU %in% sig_c
    excluded <- c(excluded, d$feature[drop & d$cOTU_override])
    d <- d[!drop, ]
    if (!nrow(d)) next
    d <- d[nat_order(d$feature), ]
    rows <- c(rows, mapply(function(...) row_tex(...),
                           LOC_LABEL[[loc]], tex_escape(d$feature), tex_tax(d$cOTU),
                           tex_tax(d$Kingdom), tex_tax(d$Phylum), tex_tax(d$Class),
                           tex_tax(d$Order),   tex_tax(d$Family), tex_genus(d$Genus),
                           tex_num(d$lfc, DIG_LFC_LONG), USE.NAMES = FALSE))
  }
  arm_lab <- tolower(ARM)
  note <- paste0("\\footnotesize ``Not detected'' denotes a feature for which the supplied mapping did not identify a cOTU; dashes indicate unavailable taxonomic ranks. A feature mapped to a cOTU can still appear in this table when that cOTU was not significant at the same location.")
  if (length(excluded)) {
    note <- paste0(note, " ", paste(tex_escape(unique(excluded)), collapse = " and "),
                   if (length(unique(excluded)) > 1) " are" else " is",
                   " excluded because \\citet{jin2022quantsinglecell} identify ",
                   if (length(unique(excluded)) > 1) "them" else "it",
                   " as ", if (length(unique(excluded)) > 1) paste0(lv, " representations")
                            else paste0("an ", lv, " representation"),
                   " of the significant ",
                   paste(tex_escape(unique(unname(MANUAL_COTU_MAP))), collapse = " and "),
                   ", although the supplied mapping labels ",
                   if (length(unique(excluded)) > 1) "them" else "it", " ``Not detected.''")
  }
  c("\\clearpage",
    "\\begin{landscape}",
    "\\begin{table}[p]",
    "    \\centering",
    "    \\scriptsize",
    "    \\setlength{\\tabcolsep}{3pt}",
    "    \\resizebox{\\linewidth}{!}{%",
    "    \\begin{tabular}{@{}lllllllllr@{}}",
    "        \\toprule",
    "        Location & Feature & Mapped cOTU & Kingdom & Phylum & Class & Order & Family & Genus & Estimated LFC \\\\",
    "        \\midrule",
    rows,
    "        \\bottomrule",
    "    \\end{tabular}%",
    "    }",
    paste0("    \\caption{Significant ", lv, "s without a significant cOTU counterpart at the same cecal location. The mapped cOTU and its taxonomic lineage are reported when available. Estimated LFCs are ",
           arm_lab, " Sparse SSRV posterior means for the vitamin A-deficient versus vitamin A-sufficient contrast on the $\\log_2$ scale.}"),
    paste0("    \\label{tab:supp_jin_", tolower(lv), "_without_sig_cotu}"),
    "    \\vspace{2pt}",
    "    \\begin{minipage}{0.98\\linewidth}",
    paste0("        ", note),
    "    \\end{minipage}",
    "\\end{table}",
    "\\end{landscape}")
}

## ---------------------------------------------------------------------------
## Supplementary: cOTUs called by exactly one of the two analyses.
##
## Because arm_out keeps every feature, the Sparse SSRV LFC and q-value are
## reported for the Jin-only cOTUs as well, rather than left as dashes.
## ---------------------------------------------------------------------------
tab_discordant_cotus <- function(ssrv_tbl, jin_tbl) {
  rows      <- character(0)
  agreed    <- character(0)
  for (loc in LOC_KEYS) {
    ours  <- ssrv_tbl[ssrv_tbl$level == "cOTU" & ssrv_tbl$location == loc, ]
    theirs <- jin_tbl[jin_tbl$level == "cOTU" & jin_tbl$location == loc, ]
    sig_ours   <- ours$feature[ours$sig]
    sig_theirs <- theirs$feature[theirs$sig]
    if (!length(setdiff(sig_ours, sig_theirs)) && !length(setdiff(sig_theirs, sig_ours))) {
        agreed <- c(agreed, sprintf("%s location, where only %s %s significant",
                                  tolower(LOC_LABEL[[loc]]),
                                  paste(tex_escape(sig_ours), collapse = " and "),
                                  if (length(sig_ours) > 1) "were" else "was"))
      next
    }
    blocks <- list("Jin et al. only"    = setdiff(sig_theirs, sig_ours),
                   "Sparse SSRV only"   = setdiff(sig_ours,   sig_theirs))
    for (b in names(blocks)) {
      f <- blocks[[b]]
      if (!length(f)) next
      i <- match(f, ours$feature)
      j <- match(f, theirs$feature)
      o <- nat_order(f)
      i <- i[o]; j <- j[o]; f <- f[o]
      if (length(rows)) rows <- c(rows, "        \\addlinespace")
      rows <- c(rows, mapply(function(...) row_tex(...),
                             LOC_LABEL[[loc]], tex_escape(f), b,
                             tex_tax(ours$Kingdom[i]), tex_tax(ours$Phylum[i]),
                             tex_tax(ours$Class[i]),   tex_tax(ours$Order[i]),
                             tex_tax(ours$Family[i]),  tex_genus(ours$Genus[i]),
                             tex_num(theirs$log2FC[j], DIG_LFC_SHORT),
                             tex_pval(theirs$FDR[j]),
                             tex_num(ours$lfc[i], DIG_LFC_SHORT),
                             tex_q(ours$q[i]),
                             USE.NAMES = FALSE))
    }
  }
  agree_sentence <- if (length(agreed)) {
    paste0(" The two analyses agreed completely at the ", paste(agreed, collapse = ", and at the "), ".")
  } else ""
  arm_lab <- tolower(ARM)
  c("\\clearpage",
    "\\begin{landscape}",
    "\\begin{table}[p]",
    "    \\centering",
    "    \\scriptsize",
    "    \\setlength{\\tabcolsep}{2.5pt}",
    "    \\resizebox{\\linewidth}{!}{%",
    "    \\begin{tabular}{@{}lllllllllrrrr@{}}",
    "        \\toprule",
    "        Location & cOTU & Significant in & Kingdom & Phylum & Class & Order & Family & Genus & \\shortstack{Jin et al.\\\\$\\log_2$ LFC} & \\shortstack{Jin et al.\\\\FDR} & \\shortstack{Sparse SSRV\\\\$\\log_2$ LFC} & \\shortstack{Sparse SSRV\\\\q-value} \\\\",
    "        \\midrule",
    rows,
    "        \\bottomrule",
    "    \\end{tabular}%",
    "    }",
    paste0("    \\caption{Location-specific cOTU discoveries that differed between \\citet{jin2022quantsinglecell} and ",
           arm_lab, " Sparse SSRV. Published DESeq2 LFCs and Sparse SSRV posterior means are both on the $\\log_2$ scale; the corresponding adjusted significance measures are the FDR and the ",
           arm_lab, " q-value respectively. Estimates and q-values are reported for every cOTU listed, including those only \\citet{jin2022quantsinglecell} called, so the two analyses can be compared on the features where they disagree.",
           agree_sentence, "}"),
    "    \\label{tab:supp_jin_discordant_cotus}",
    "    \\vspace{2pt}",
    "    \\begin{minipage}{0.98\\linewidth}",
    paste0("        \\footnotesize Jin et al. significance was defined by FDR $<", JIN_FDR_MAX,
           "$ and absolute fold change $>", 2^JIN_ABS_LOG2FC_MIN,
           "$. Raw p-values were not included in the published source-data table, so FDR values are reported. A Sparse SSRV q-value shown as $<10^{-4}$ was zero to the resolution of the ",
           format(attr(results_all$cOTU$prox$sparse_out, "n_draws"), big.mark = ","),
           " posterior draws."),
    "    \\end{minipage}",
    "\\end{table}",
    "\\end{landscape}")
}

################################################################################
## EMIT
################################################################################

TABLES <- list(
  jin_scale_model        = tab_scale_model(runs),
  jin_discovery_counts   = tab_discovery_counts(ssrv, jin),
  jin_cross_representation = tab_cross_representation(ssrv, jin),
  jin_cross_representation_full = tab_cross_representation_full(ssrv, jin, scope = "cotu"),
  jin_cross_representation_sig_features = tab_cross_representation_full(ssrv, jin, scope = "feature"),
  jin_all_significant    = tab_all_significant(ssrv),
  jin_asv_without_sig_cotu = tab_without_sig_cotu(ssrv, "ASV"),
  jin_otu_without_sig_cotu = tab_without_sig_cotu(ssrv, "OTU"),
  jin_discordant_cotus   = tab_discordant_cotus(ssrv, jin)
)

preamble_note <- c(
  "% ---------------------------------------------------------------------------",
  "% Jin et al. (2022) tables -- GENERATED FILE, do not edit by hand.",
  paste0("% Produced by jin_DA_biological_interpretation.R on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "."),
  paste0("% Sparse SSRV result reported: ", ARM,
         if (identical(ARM, "hybrid")) " (selection rule's choice, per subset)" else "", ".  Analysis objects: jin_results_all_levels_", DATESTAMP, ".RData."),
  paste0("% Published comparison: ", JIN_SOURCE_BOOK, ", Supplementary Figure 18a-f."),
  "% All log fold changes, on both sides, are log2.",
  "%",
  "% Requires in the preamble:",
  "%   \\usepackage{booktabs, longtable, pdflscape, graphicx, multirow, array}",
  "%   \\usepackage[authoryear]{natbib}   % for \\citet",
  "% ---------------------------------------------------------------------------",
  ""
)

for (nm in names(TABLES)) {
  writeLines(TABLES[[nm]], file.path(DIR_TAB, paste0(nm, arm_suffix(), ".tex")))
}

out_tex <- file.path(DIR_JIN, paste0("jin_tables_", DATESTAMP, arm_suffix(), ".tex"))
writeLines(c(preamble_note,
             unlist(lapply(names(TABLES), function(nm)
               c(paste0("% ---- ", nm, " ----"), TABLES[[nm]], "")),
               use.names = FALSE)),
           out_tex)
message("Wrote ", out_tex)
message("Wrote ", length(TABLES), " individual table files to ", DIR_TAB)

## Echo the whole thing when run non-interactively, so `Rscript this.R` prints
## paste-ready LaTeX. Suppressed when manuscript_tables.R sources this file:
## that script owns stdout and re-emits the tables it actually uses.
if (!interactive() && !isTRUE(getOption("jin.quiet"))) cat(readLines(out_tex), sep = "\n")

################################################################################
## SIDE OUTPUT -- the significant-feature table with taxonomy and LFC, kept for
## downstream use. Same shape as before; `cOTU` now carries the mapped cOTU at
## every level (the feature itself at cOTU level) instead of the taxon label,
## so joining on it is safe.
################################################################################

sig_tax_est <- lapply(setNames(LOC_KEYS, LOC_KEYS), function(loc) {
  d <- ssrv[ssrv$sig & ssrv$location == loc, ]
  d <- d[order(match(d$level, LEVEL_KEYS), nat_order(d$feature)), ]
  data.frame(level = d$level, location = d$location,
             ## sparse_SSRV is the hybrid in this project's method naming;
             ## the single-arm runs carry the _B / _L suffixes.
             method = switch(ARM, hybrid = "sparse_SSRV",
                             bootstrap = "sparse_SSRV_B", laplace = "sparse_SSRV_L"),
             category = d$feature,
             Kingdom = d$Kingdom, Phylum = d$Phylum, Class = d$Class,
             Order = d$Order, Family = d$Family, Genus = d$Genus,
             cOTU = d$cOTU, estimated_LFC = d$lfc,
             q_value = d$q, stringsAsFactors = FALSE, row.names = NULL)
})

f_tax <- file.path(DIR_JIN, paste0("jin_real_data_tax_estLFC", arm_suffix(), ".RData"))
saveRDS(sig_tax_est, file = f_tax)
message("Wrote ", f_tax)

################################################################################
## CONSOLE SUMMARY
################################################################################
check_feature_table(TABLES$jin_cross_representation_sig_features, ssrv,
                    exclude_cotu_scope = TRUE, jin_tbl = jin)
report_marked_cotus(ssrv, jin)

message("\nSignificant features by level and location (", ARM, " arm):")
print(with(ssrv[ssrv$sig, ], table(level = factor(level, LEVEL_KEYS),
                                   location = factor(location, LOC_KEYS))))
message("\nPublished (Jin et al.) counts under FDR < ", JIN_FDR_MAX,
        " and |log2FC| > ", JIN_ABS_LOG2FC_MIN, ":")
print(with(jin[jin$sig, ], table(level = factor(level, LEVEL_KEYS),
                                 location = factor(location, LOC_KEYS))))
