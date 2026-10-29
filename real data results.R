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

## Study subdirectories are spelled once here rather than repeated on every
## readRDS() line, so a rename only has to be made in one place.
DIR_SELEX      <- file.path(ROOT, "selex")
DIR_ZEMB       <- file.path(ROOT, "2020_zemb_microOpen_spike")
DIR_STAMMLER   <- file.path(ROOT, "2016_stammler_microbiome_micehuman")
DIR_JIN        <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates")
DIR_PRASSE     <- file.path(ROOT, "2015_prasse_appliedEnvMicrobiology_wetlands")
DIR_VANDEPUTTE <- file.path(ROOT, "Vandeputte")

## collect_scale_records(): reads the scale_model_record.rds each DA script
## writes, so the scale-model choice and both arms' significant counts come out
## of the same run that produced the TP/FP table below.
## arm_tpfp(): reads sparse_arm_sig.rds from the same run and scores the
## bootstrap and laplace scale-model arms separately against a study's truth.
## collect_arm_features(): the same records at FEATURE level -- which taxa each
## scale model called, not just how many -- for the last section of this file.
source(file.path(ROOT, "scale_model_record.R"))

################################################################################
## DATESTAMP RESOLUTION
##
## Every study writes its outputs under the datestamp of the run that produced
## them, and this file used to name "03212026" in fourteen places. Once a study
## was re-run under a new stamp, the re-run wrote new files and this file went on
## reading the old ones -- silently, because the old ones still exist. That is
## how Zemb came to report 10 true positives against a truth set of 1: the
## March fit had been scored against a ten-taxon gold reference, the August fit
## scores against the spike-in alone, and the table was still reading March.
##
## stamped() resolves each file to the newest stamp that actually exists, records
## which it used, and report_stamps() prints the result. A study whose files sit
## at different stamps, or a table built from a mixture of runs, is now visible
## rather than something you have to notice in a number.
################################################################################
RESULT_STAMPS <- c("09012026", "08242026", "03212026")
.stamp_log <- new.env(parent = emptyenv())

## `pattern` carries one %s where the datestamp goes.
stamped <- function(dir, pattern, what = pattern) {
  for (st in RESULT_STAMPS) {
    f <- file.path(dir, sprintf(pattern, st))
    if (file.exists(f)) {
      assign(paste(basename(dir), what), st, envir = .stamp_log)
      if (identical(what, "res_tpfp")) assign(paste0("PATH ", basename(dir)), f, envir = .stamp_log)
      return(f)
    }
  }
  stop("no file for ", what, " in ", dir, " under any known datestamp (",
       paste(RESULT_STAMPS, collapse = ", "), ")")
}

## A datestamp is not enough to tell two runs apart. Every DA script calls
## record_scale_model() partway through and saves res_tpfp at the very end, so a
## run that was interrupted between the two leaves fresh per-arm records beside a
## stale method table -- both under the same datestamp, and the resulting row is
## half one run and half another. That is what put Zemb's hybrid column at odds
## with the arm it is supposed to equal. Compared here by clock time, since the
## filename cannot show it.
## `runs` is a named list: directory -> the pair this table reports. Comparing
## against the pair matters -- Jin fits six subsets and Stammler fifteen pool
## pairs in one pass, so the LAST subset's record is legitimately later than the
## reported subset's method table, and comparing against the maximum would flag
## every multi-subset study on every run.
check_run_coherence <- function(runs) {
  bad <- character(0)
  for (i in seq_along(runs)) {
    d <- names(runs)[i]; pr <- runs[[i]]
    f_arm <- file.path(d, "sparse_arm_sig.rds")
    f_res <- get0(paste0("PATH ", basename(d)), envir = .stamp_log)
    if (!file.exists(f_arm) || is.null(f_res)) next
    a <- readRDS(f_arm)
    key <- paste0(basename(d), "|")   # placeholder, resolved below
    hit <- Filter(function(z) {
      zp <- if (is.null(z$pair) || is.na(z$pair)) NA_character_ else as.character(z$pair)
      (is.na(pr) && is.na(zp)) || (!is.na(pr) && !is.na(zp) && zp == pr)
    }, a)
    if (!length(hit)) next
    ra <- suppressWarnings(max(as.POSIXct(vapply(hit, function(z)
      if (is.null(z$run_at)) NA_character_ else z$run_at, character(1))), na.rm = TRUE))
    rr <- file.mtime(f_res)
    if (is.finite(ra) && !is.na(rr) && ra > rr + 60) {
      bad <- c(bad, sprintf("%s%s: per-arm record %s is newer than %s (%s)",
                            basename(d), if (is.na(pr)) "" else paste0(" / ", pr),
                            format(ra, "%Y-%m-%d %H:%M"),
                            basename(f_res), format(rr, "%Y-%m-%d %H:%M")))
    }
  }
  if (length(bad)) {
    warning("a later run wrote per-arm records but not the method table for:\n  ",
            paste(bad, collapse = "\n  "),
            "\nThat run did not reach its res_tpfp save, so these rows mix two fits. ",
            "Re-run those DA scripts to completion.", call. = FALSE, immediate. = TRUE)
  }
  invisible(bad)
}

report_stamps <- function() {
  k <- setdiff(ls(.stamp_log), grep("^PATH ", ls(.stamp_log), value = TRUE))
  if (!length(k)) return(invisible(NULL))
  d <- data.frame(file = k, stamp = unlist(mget(k, envir = .stamp_log)),
                  row.names = NULL, stringsAsFactors = FALSE)
  d <- d[order(d$stamp, d$file), ]
  cat("\n== datestamp actually read, per file ==\n")
  print(d, row.names = FALSE)
  ## A study reading two stamps at once is the dangerous case: its method rows
  ## and its per-arm rows then describe different fits, and the hybrid column
  ## stops matching the arm it is supposed to equal. Called out separately from
  ## the milder case of one study lagging the others.
  d$study <- sub(" [^ ]+$", "", d$file)
  split_stamps <- tapply(d$stamp, d$study, function(z) length(unique(z)))
  mixed_study <- names(split_stamps)[split_stamps > 1]
  if (length(mixed_study)) {
    warning("these studies are read from MORE THAN ONE run at once: ",
            paste(mixed_study, collapse = ", "),
            ". Their method-level and per-arm rows describe different fits, so the ",
            "hybrid column need not match either arm. Re-run them so every file ",
            "carries one stamp.", call. = FALSE, immediate. = TRUE)
  }
  if (length(unique(d$stamp)) > 1) {
    old <- unique(d$stamp)[unique(d$stamp) != RESULT_STAMPS[1]]
    warning("this table mixes runs: ", sum(d$stamp != RESULT_STAMPS[1]), " of ",
            nrow(d), " files are still at ", paste(old, collapse = "/"),
            " because no newer version exists. Re-run those studies to make the ",
            "comparison uniform.", call. = FALSE, immediate. = TRUE)
  }
  invisible(d)
}

## Jin does not produce ONE result. Its script runs the VS-vs-VD contrast per
## feature level (cOTU / ASV / OTU) and per subset -- pooled ("both") plus the two
## cecal locations (prox / dist) -- writing a tagged file for each, so this table
## has to say which one it is reporting.
##
## THE PRIMARY SUBSET IS cOTU / POOLED, and the reason is the gold reference.
## Every figure in this table is a TP/FP count against the measured-scale
## reference, and at 4 vs 4 BOTH stratified references are EMPTY at padj<0.05 --
## not marginal, empty. Scoring against an empty set makes every call a false
## positive and every method score 0 TP, which says something about the
## reference's power at n=4, not about the methods. That is why the Jin block
## below sat commented out. Pooled (8 vs 8) is the smallest subset of this design
## with a non-empty reference, so it is the one the main table reports.
##
## The stratified runs remain the SUPPLEMENTARY analysis and are unaffected by
## this choice: manuscript_tables.R builds them through
## jin_DA_biological_interpretation.R, which enumerates LOC_KEYS = prox/dist and
## never reads the pooled files.
##
## Both halves use the SAME filter rule (>=5 raw counts in >=3 samples, recomputed
## per subset), so the primary tests 189 cOTUs and the supplementary 150 / 153.
##
## JIN_TAG feeds the file names -- "" for the pooled cOTU run, which keeps the
## bare tax_scale_<stamp>.RData / res_tpfp_<stamp>.RData names, and "_prox" /
## "_dist" for the stratified ones. JIN_PAIR is the matching key inside Jin's
## scale_model_record.rds / sparse_arm_sig.rds, which sparse_DA_09012026.R writes
## as the bare subset name for cOTU and "<level>_<subset>" for ASV and OTU.
##
## Set these to "_prox"/"prox" to report a stratified run instead -- but expect
## every TP to be 0 until that run has a non-empty reference.
JIN_TAG  <- ""
JIN_PAIR <- "both"

## Prasse likewise produces more than ONE result, and they are NOT interchangeable.
## Several runs exist, all writing into the same directory under the dataset name
## "Prasse2015", so the run key is the only thing that tells them apart. The first
## three are the CURRENT rowMeans(raw counts) >= 5 family:
##
##   rawVSqseq_QSeq  sparse_DA_09012026_rawVSqseq_mean5.R, the QSeq pass. THE
##                   REPORTED RUN. The methods are fitted on the QSeq-rescaled
##                   counts -- what the single-pass script also does -- while the
##                   FILTER is decided once on the raw counts and shared by both
##                   passes.
##   rawVSqseq_raw   the same script's raw-count pass. Under the August
##                   prevalence filter its gold reference was EMPTY at the
##                   reported gamma (1 / 0 / 0 / 0 across the four), so scoring
##                   against it gave every method TP = 0 and called every
##                   discovery a false positive -- it measured the reference's
##                   power, not the methods. Whether the mean filter changes that
##                   is an OPEN QUESTION until this run exists; the block below
##                   still raises its warning before reporting it, and that
##                   warning is now the thing to read rather than a formality.
##   mean5           sparse_DA_09012026_mean5.R, one pass.
##
## and the rest are SUPERSEDED, kept readable so an August or March fit can be
## pulled up deliberately. Selecting one also means setting
## PRASSE_RECORD_FILTER below, or the record key will not match the files:
##
##   prev10_QSeq     the August rawVSqseq QSeq pass (10% prevalence filter).
##   prev10          the August single-pass run.
##   sum10           the older sum(counts) > 10 filter.
##
## FILENAMES ARE PATTERNS, with %s where the datestamp goes, resolved by
## prasse_file() to the newest stamp that exists. They used to be literals naming
## 03212026, so after the August re-runs this section went on reading the March
## files while the per-arm records it sits beside came from August -- the exact
## split-run failure check_run_coherence() further down exists to catch, made
## invisible because the filename never changed.
##
## `input` is the extra level the rawVSqseq object carries. That script fits both
## passes and saves ONE res_tpfp keyed [[INPUT]][[gamma]], where the single-pass
## scripts save [[gamma]]; non-NA `input` selects that level. Its sens and gold
## files are per-pass and so need no such indexing.
##
## `ckf` is NA everywhere now: every current run folds CKF into res_tpfp itself
## through add_ckf_row(), so it arrives as a row of aldex_gold like any other
## method. It is kept for the sum10 run, whose CKF still lives in its own file.
PRASSE_RUNS <- list(
  rawVSqseq_QSeq = c(sens  = "sparse_sens_%s_mean5_QSeq.RData",
                     gold  = "tax_scale_%s_mean5_QSeq.RData",
                     tpfp  = "res_tpfp_%s_mean5_rawVSqseq.RData",
                     ckf   = NA_character_,
                     input = "QSeq"),
  rawVSqseq_raw  = c(sens  = "sparse_sens_%s_mean5_raw.RData",
                     gold  = "tax_scale_%s_mean5_raw.RData",
                     tpfp  = "res_tpfp_%s_mean5_rawVSqseq.RData",
                     ckf   = NA_character_,
                     input = "raw"),
  mean5          = c(sens  = "sparse_sens_%s_mean5.RData",
                     gold  = "tax_scale_%s_mean5.RData",
                     tpfp  = "res_tpfp_%s_mean5.RData",
                     ckf   = NA_character_,
                     input = NA_character_),
  ## SUPERSEDED runs, kept so an older fit can still be read deliberately. Their
  ## files only exist under the 08242026 / 03212026 stamps.
  prev10_QSeq    = c(sens  = "sparse_sens_%s_prev10_QSeq.RData",
                     gold  = "tax_scale_%s_prev10_QSeq.RData",
                     tpfp  = "res_tpfp_%s_prev10_rawVSqseq.RData",
                     ckf   = NA_character_,
                     input = "QSeq"),
  prev10         = c(sens  = "sparse_sens_%s_prev10.RData",
                     gold  = "tax_scale_%s_prev10.RData",
                     tpfp  = "res_tpfp_%s_prev10.RData",
                     ckf   = NA_character_,
                     input = NA_character_),
  sum10          = c(sens  = "sparse_sens_%s.RData",
                     gold  = "tax_scale_%s_filter0_1.RData",
                     tpfp  = "res_tpfp_%s_filter0_1.RData",
                     ckf   = "prasse_ckf_sensitivity_%s_filter0_1.RData",
                     input = NA_character_))

## The stamps Prasse's runs carry, newest first. 09012026 is the rowMeans >= 5
## re-run; 08262026 is here because the August rawVSqseq script was renamed from
## _08262026_ to _08242026_ AFTER its last run, so its outputs still sit under
## that stamp. Because the current patterns carry a _mean5 suffix and the
## superseded ones a _prev10 suffix, a mean5 pattern CANNOT silently resolve to
## an August file: the name would not exist. The stamp list is the fallback for
## the sum10 run, whose pattern is stamp-only.
PRASSE_STAMPS <- c("09012026", "08242026", "08262026", "03212026")

prasse_file <- function(pattern, what = pattern) {
  if (is.na(pattern)) return(NA_character_)
  for (st in PRASSE_STAMPS) {
    f <- file.path(DIR_PRASSE, sprintf(pattern, st))
    if (file.exists(f)) {
      assign(paste("2015_prasse_appliedEnvMicrobiology_wetlands", what), st, envir = .stamp_log)
      if (identical(what, "res_tpfp"))
        assign("PATH 2015_prasse_appliedEnvMicrobiology_wetlands", f, envir = .stamp_log)
      return(f)
    }
  }
  stop("no Prasse file matching '", pattern, "' under any known stamp (",
       paste(PRASSE_STAMPS, collapse = ", "), ")")
}

## The run this table reports. Doubles as the key into Prasse's
## scale_model_record.rds / sparse_arm_sig.rds, so the scale-model choice, the
## per-arm TP/FP rows and the TP/FP table all come from the same fit -- the
## rawVSqseq script records its two passes as "mean5_QSeq" and "mean5_raw".
PRASSE_PAIR  <- "rawVSqseq_QSeq"
PRASSE_FILES <- PRASSE_RUNS[[PRASSE_PAIR]]
PRASSE_INPUT <- unname(PRASSE_FILES[["input"]])

## The record key, which is not always the run name: the rawVSqseq passes are
## recorded under "<FILTER>_<INPUT>" because that is what the DA script's PAIR is.
## PRASSE_RECORD_FILTER tracks the filter suffix in the DA script's `PAIR <- ...`
## line -- "mean5" since 2026-09-01, "prev10" for the archived August fits. It is
## a separate knob from PRASSE_PAIR so that selecting a superseded run above also
## needs its record key selected here, out loud, rather than reading a mean5
## record against a prev10 result table.
PRASSE_RECORD_FILTER <- "mean5"
PRASSE_KEY <- if (is.na(PRASSE_INPUT)) PRASSE_PAIR else
  paste0(PRASSE_RECORD_FILTER, "_", PRASSE_INPUT)

#------------------------------------------------------------------------------------
# SELEX
#------------------------------------------------------------------------------------
sparse_sens_selex <- readRDS(stamped(DIR_SELEX, "sparse_sens_%s.RData", "sparse_sens"))
res_selex <- readRDS(stamped(DIR_SELEX, "res_tpfp_%s.RData", "res_tpfp"))
ckf_selex <- readRDS(file.path(DIR_SELEX, "selex_mrckf_out.RData"))
list_taxa_selex <- readRDS(stamped(DIR_SELEX, "results_%s.RData", "results"))
selex_truth <- read.table(file.path(DIR_SELEX, "selex-truth.txt"))
res_selex_sens <- lapply(sparse_sens_selex[c("nrd","SJ-ste")], FUN = function(x){c(sum(rownames(x) %in% rownames(selex_truth)), sum(!rownames(x) %in% rownames(selex_truth)))})
res_selex_sens <- do.call(rbind, res_selex_sens)
rownames(res_selex_sens) <- paste0("sparse_SSRV_", rownames(res_selex_sens))
## selex's CKF figure is WITHHELD rather than reported.
##
## selex_mrckf_out.RData predates the fix to selex_mrckf.R, where the FP line was
## missing its "!" and so computed FP with the same expression as TP. The stored
## object still shows it: identical(tpfp[,1], tpfp[,2]) is TRUE across all 100
## replicates, and the mean is 0.74/0.74. That number is the true-positive count
## printed twice, not a false-positive rate.
##
## The taxon-filter work did not change selex's CKF inputs, so the file was left
## alone rather than recomputed -- at p=1600, n=14 a local rerun is ~30h. Rerun
## selex_mrckf.R (now fixed) on Roar and delete these two lines to restore it.
ckf_selex$mean <- c(TP = NA_real_, FP = NA_real_)
res_selex$aldex_gold_full <- rbind(res_selex$aldex_gold, CKF = ckf_selex$mean, res_selex_sens,
                                   arm_tpfp(DIR_SELEX, "McMurrough2014", rownames(selex_truth)))

bws <- names(sparse_sens_selex)

#------------------------------------------------------------------------------------
# Zemb
#------------------------------------------------------------------------------------
sparse_sens_zemb <- readRDS(stamped(DIR_ZEMB, "sparse_sens_%s.RData", "sparse_sens"))
res_zemb <- readRDS(stamped(DIR_ZEMB, "res_tpfp_%s.RData", "res_tpfp"))
ckf_zemb <- readRDS(file.path(DIR_ZEMB, "zemb_mrckf_out.RData"))
zemb_truth <- "Seq1"
res_zemb_sens <- lapply(sparse_sens_zemb[c("nrd","SJ-ste")], FUN = function(x){c(sum(rownames(x) %in% zemb_truth), sum(!rownames(x) %in% zemb_truth))})
res_zemb_sens <- do.call(rbind, res_zemb_sens)
rownames(res_zemb_sens) <- paste0("sparse_SSRV_", rownames(res_zemb_sens))
res_zemb$aldex_gold_full <- rbind(res_zemb$aldex_gold, CKF=ckf_zemb$mean, res_zemb_sens,
                                  arm_tpfp(DIR_ZEMB, "Zemb2020", zemb_truth))
rownames(res_zemb$aldex_gold_full)[c(9,11)] <- c("ANCOM-BC", "ANCOM-BC2")


#------------------------------------------------------------------------------------
# Stammler
#------------------------------------------------------------------------------------
sparse_sens_stammler <- readRDS(stamped(DIR_STAMMLER, "sparse_sens_%s.RData", "sparse_sens"))
res_stammler <- readRDS(stamped(DIR_STAMMLER, "res_tpfp_%s.RData", "res_tpfp"))
ckf_stammler <- readRDS(file.path(DIR_STAMMLER, "stammler_mrckf_out_all_pairs.RDS"))
stammler_truth <- c("AB076660XXXX", "AB247615XXXX")
res_stammler_sens <- lapply(sparse_sens_stammler[c("nrd","SJ-ste")], FUN = function(x){c(sum(rownames(x) %in% stammler_truth), sum(!rownames(x) %in% stammler_truth))})
res_stammler_sens <- do.call(rbind, res_stammler_sens)
rownames(res_stammler_sens) <- paste0("sparse_SSRV_", rownames(res_stammler_sens))
res_stammler$aldex_gold_full <- rbind(res_stammler$aldex_gold, CKF = colMeans(ckf_stammler[[5]]), res_stammler_sens,
                                      arm_tpfp(DIR_STAMMLER, "Stammler2016", stammler_truth, pair = "Pool6"))
rownames(res_stammler$aldex_gold_full)[rownames(res_stammler$aldex_gold_full)=="Sparse_SSRV"] <- "sparse_SSRV"

#------------------------------------------------------------------------------------
# Jin
#------------------------------------------------------------------------------------
message("Jin: reporting the cOTU / '", JIN_PAIR, "' run")

## GUARD -- the pooled run's files carry NO subset tag, and files with those exact
## names already exist at the older datestamp from the pre-split script, which
## agglomerated to genus and tested 132 features rather than 189 cOTUs. stamped()
## falls back through RESULT_STAMPS, so without this check a missing re-run would
## be reported silently as if it were the current pooled analysis -- the one
## failure this whole block is arranged to avoid.
##
## If this stops: run sparse_DA_09012026.R with POOLED_LEVELS <- "cOTU" (its
## default), which writes the bare-tag files at the current datestamp.
f_jin_tpfp <- file.path(DIR_JIN, sprintf(paste0("res_tpfp", JIN_TAG, "_%s.RData"),
                                         RESULT_STAMPS[1]))
if (!file.exists(f_jin_tpfp)) {
  stop("Jin's '", JIN_PAIR, "' run has no res_tpfp at the current datestamp (",
       RESULT_STAMPS[1], "):\n  ", f_jin_tpfp,
       "\nOlder files with this name exist but come from the pre-split, ",
       "genus-agglomerated script and are NOT this analysis. Re-run ",
       "sparse_DA_09012026.R (POOLED_LEVELS <- \"cOTU\") before building this table.")
}

sparse_sens_jin <- readRDS(stamped(DIR_JIN, paste0("sparse_sens", JIN_TAG, "_%s.RData"), "sparse_sens"))
res_jin         <- readRDS(stamped(DIR_JIN, paste0("res_tpfp",    JIN_TAG, "_%s.RData"), "res_tpfp"))

## The gold reference SAVED by that run, the same way Prasse reads its own: this
## is provably the set res_tpfp was scored against, rather than a tax_scale file
## that happens to sit in the directory.
tax_scale_jin <- readRDS(stamped(DIR_JIN, paste0("tax_scale", JIN_TAG, "_%s.RData"), "tax_scale"))
jin_truth     <- tax_scale_jin$category

## An empty reference makes every TP 0 and every call a false positive, which is
## a statement about the reference's power and not about the methods. That is the
## condition that kept this block commented out while the primary subset was a
## single 4-vs-4 location; say so loudly rather than printing a column of zeros.
if (!length(jin_truth)) {
  warning("Jin's '", JIN_PAIR, "' measured-scale reference is EMPTY -- every TP ",
          "in the Jin column will be 0 by construction. Do not read that column ",
          "as a method comparison.", call. = FALSE, immediate. = TRUE)
} else {
  message("  gold (measured-scale) reference: ", length(jin_truth), " cOTUs")
}

## NO `CKF = ...` row is added here, unlike Zemb and selex above.
## sparse_DA_09012026.R already merges it with add_ckf_row() before saving
## res_tpfp, so res_jin$aldex_gold arrives with its CKF row in place. Re-adding it
## from the legacy jin_mrckf_out.RData would stack a second one -- and that legacy
## file is written by ONE subset only, so it need not be the subset reported here.
res_jin_sens <- lapply(sparse_sens_jin[c("nrd","SJ-ste")], FUN = function(x){c(sum(rownames(x) %in% jin_truth), sum(!rownames(x) %in% jin_truth))})
res_jin_sens <- do.call(rbind, res_jin_sens)
rownames(res_jin_sens) <- paste0("sparse_SSRV_", rownames(res_jin_sens))
res_jin$aldex_gold_full <- rbind(res_jin$aldex_gold, res_jin_sens,
                                 arm_tpfp(DIR_JIN, "Jin2022", jin_truth, pair = JIN_PAIR))
## Renamed BY NAME, not by position: this matrix now has the CKF row merged in
## upstream and two arm rows appended, so the positional c(9,11) this line used to
## carry no longer points where it did.
rownames(res_jin$aldex_gold_full)[rownames(res_jin$aldex_gold_full) == "ANCOMBC"]  <- "ANCOM-BC"
rownames(res_jin$aldex_gold_full)[rownames(res_jin$aldex_gold_full) == "ANCOMBC2"] <- "ANCOM-BC2"

#------------------------------------------------------------------------------------
# Prasse
#------------------------------------------------------------------------------------
message("Prasse: reporting the '", PRASSE_PAIR, "' run (record key '", PRASSE_KEY, "')")
sparse_sens_prasse <- readRDS(prasse_file(PRASSE_FILES[["sens"]], "sparse_sens"))
## The gold reference SAVED by that run, not the tax_scale_*.RData it read in:
## the saved copy carries the gold_gamma_* names and is guaranteed to be the set
## res_tpfp was scored against.
tax_scale_prasse <- readRDS(prasse_file(PRASSE_FILES[["gold"]], "tax_scale"))
res_prasse_obj <- readRDS(prasse_file(PRASSE_FILES[["tpfp"]], "res_tpfp"))
## One level deeper for the two-pass rawVSqseq object; see `input` above.
res_prasse <- if (is.na(PRASSE_INPUT)) res_prasse_obj else res_prasse_obj[[PRASSE_INPUT]]
golden_ref_prasse <- tax_scale_prasse[[4]]

## An empty gold reference is not a result. Scored against it every method takes
## TP = 0 and every call becomes a false positive, which measures the reference's
## power at this gamma and not the methods -- the same reason Jin's stratified
## locations are not scored. The raw-count pass is in exactly that position, so
## reporting it is refused here rather than producing a row of zeros that reads
## as "no method found anything".
if (!nrow(golden_ref_prasse)) {
  stop("the gold reference for Prasse run '", PRASSE_PAIR, "' is EMPTY at the ",
       "reported gamma, so no method can score a true positive against it. ",
       "Report a run whose reference is non-empty, or score this one at a gamma ",
       "where it is (gold sizes by gamma: ",
       paste(vapply(tax_scale_prasse, nrow, integer(1)), collapse = ", "), ").")
}
message("  gold (scale-aware ALDEx2) reference: ", nrow(golden_ref_prasse), " genera")

res_prasse_sens <- lapply(sparse_sens_prasse[c("nrd","SJ-ste")], FUN = function(x){c(sum(rownames(x) %in% golden_ref_prasse$category), sum(!rownames(x) %in% golden_ref_prasse$category))})
res_prasse_sens <- do.call(rbind, res_prasse_sens)
rownames(res_prasse_sens) <- paste0("sparse_SSRV_", rownames(res_prasse_sens))

## CKF. Every current Prasse run fits CKF inside the DA script and folds it into
## res_tpfp through add_ckf_row(), so aldex_gold already carries a "CKF" row and
## rbind-ing a second one would leave two rows of that name -- which is not an
## error, and that is the problem: `m[method_list, ]` silently takes the first,
## so which CKF the table reports would depend on row order. Only the older runs,
## whose CKF sits in its own file, need the rbind.
has_ckf_row <- "CKF" %in% rownames(res_prasse[[4]]$aldex_gold)
if (has_ckf_row) {
  if (!is.na(PRASSE_FILES[["ckf"]]))
    warning("Prasse run '", PRASSE_PAIR, "' carries a CKF row in res_tpfp AND names ",
            "a separate CKF file; the row in res_tpfp is used and the file ignored.",
            call. = FALSE, immediate. = TRUE)
  ckf_prasse <- NULL
  ckf_extra  <- NULL                       # nothing to rbind
} else if (is.na(PRASSE_FILES[["ckf"]])) {
  warning("Prasse CKF withheld: run '", PRASSE_PAIR, "' has no CKF row in its ",
          "results object and names no CKF file.", call. = FALSE, immediate. = TRUE)
  ckf_prasse <- NULL
  ckf_extra  <- rbind(CKF = c(TP = NA_real_, FP = NA_real_))
} else {
  ckf_prasse <- readRDS(prasse_file(PRASSE_FILES[["ckf"]], "prasse_ckf"))
  ckf_extra  <- rbind(CKF = ckf_prasse[[4]]$mean)
}

res_prasse[[4]]$aldex_gold_full <- rbind(res_prasse[[4]]$aldex_gold, ckf_extra, res_prasse_sens,
                                         arm_tpfp(DIR_PRASSE, "Prasse2015", golden_ref_prasse$category,
                                                  pair = PRASSE_KEY))

#------------------------------------------------------------------------------------
# Vandeputte
#------------------------------------------------------------------------------------
sparse_sens_vandeputte <- readRDS(stamped(DIR_VANDEPUTTE, "sparse_sens_%s.RData", "sparse_sens"))
res_vandeputte <- readRDS(stamped(DIR_VANDEPUTTE, "res_tpfp_%s.RData", "res_tpfp"))
vandeputte_truths <- readRDS(stamped(DIR_VANDEPUTTE, "results_%s.RData", "results"))$ALDEx2_gold$category
res_vandeputte_sens <- lapply(sparse_sens_vandeputte[c("nrd","SJ-ste")], FUN = function(x){c(sum(rownames(x) %in% vandeputte_truths), sum(!rownames(x) %in% vandeputte_truths))})
res_vandeputte_sens <- do.call(rbind, res_vandeputte_sens)
ckf_vandeputte <- readRDS(file.path(DIR_VANDEPUTTE, "vandeputte_mrckf_out.RData"))
rownames(res_vandeputte_sens) <- paste0("sparse_SSRV_", rownames(res_vandeputte_sens))
res_vandeputte$aldex_gold_full <- rbind(res_vandeputte$aldex_gold, CKF = ckf_vandeputte$mean, res_vandeputte_sens,
                                        arm_tpfp(DIR_VANDEPUTTE, "Vandeputte2017", vandeputte_truths))


## sparse_SSRV_B and sparse_SSRV_L sit immediately after sparse_SSRV: the
## selected arm first, then what each scale model would have given on its own.
## sparse_SSRV always equals whichever of the two was selected, so a row where
## it matches neither means the record and the fit came from different runs.
method_list <- c("sparse_SSRV", "sparse_SSRV_B", "sparse_SSRV_L",
                 "ALDEx2_def", "DESeq2", "edgeR", "limma", "baySeq", "MiMIX",
                 "ANCOM-BC2", "ANCOM-BC", "ANCOM", "LinDA", "CKF",
                 "sparse_SSRV_SJ-ste","sparse_SSRV_nrd")

## Every row `method_list` asks for must exist before the table indexes it, or the
## column build fails with a bare "subscript out of bounds" that says nothing about
## which run is at fault. The row that actually goes missing is CKF: res_tpfp files
## written before add_ckf_row() existed carry twelve method rows and no CKF, so a
## stale file lands here rather than at the datestamp guard above.
jin_missing <- setdiff(method_list, rownames(res_jin$aldex_gold_full))
if (length(jin_missing)) {
  stop("Jin's '", JIN_PAIR, "' res_tpfp is missing table row(s): ",
       paste(jin_missing, collapse = ", "),
       "\nThis is a res_tpfp from an older script version. Re-run ",
       "sparse_DA_09012026.R for this subset.")
}

real_fptp_out <- t(cbind(
selex = paste0(res_selex$aldex_gold_full[method_list,2],"/",paste0("(",res_selex$aldex_gold_full[method_list,1],")")),
zemb = paste0(res_zemb$aldex_gold_full[method_list,2],"/",paste0("(",res_zemb$aldex_gold_full[method_list,1],")")),
stammler = paste0(res_stammler$aldex_gold_full[method_list,2],"/",paste0("(",res_stammler$aldex_gold_full[method_list,1],")")),
jin = paste0(res_jin$aldex_gold_full[method_list,2],"/",paste0("(",res_jin$aldex_gold_full[method_list,1],")")),
prasse = paste0(res_prasse[[4]]$aldex_gold_full[method_list,2],"/",paste0("(",res_prasse[[4]]$aldex_gold_full[method_list,1],")")),
vandeputte = paste0(res_vandeputte$aldex_gold_full[method_list,2],"/",paste0("(",res_vandeputte$aldex_gold_full[method_list,1],")"))
))
colnames(real_fptp_out) <- method_list

## Render withheld cells as "--" rather than "NA/(NA)", so a deliberately absent
## figure cannot be mistaken for a computed one. Only cells that are NA in BOTH
## TP and FP are blanked.
real_fptp_out[real_fptp_out == "NA/(NA)"] <- "--"

real_fptp_out

## Persist it. manuscript_tables.R reads this file to build the results tables,
## and until now it was only ever built in memory here, so the copy on disk was
## whatever had last been saved by hand -- which is exactly how those tables end
## up reporting a run nobody can identify. Saved next to the other records this
## script writes, so one run of this file refreshes all of them together.
## Which run each cell of the table above actually came from. Printed next to
## the table rather than at the end of the log, because a mixed-stamp table is a
## caveat on THESE numbers.
report_stamps()
check_run_coherence(setNames(
  list(NA_character_, NA_character_, "Pool6", JIN_PAIR, PRASSE_KEY, NA_character_),
  c(DIR_SELEX, DIR_ZEMB, DIR_STAMMLER, DIR_JIN, DIR_PRASSE, DIR_VANDEPUTTE)))

saveRDS(real_fptp_out, file = file.path(ROOT, "real_data_results.RData"))
write.csv(real_fptp_out, file = file.path(ROOT, "real_data_results.csv"))


#------------------------------------------------------------------------------------
# Scale model: which arm was selected, and what the OTHER arm would have found
#------------------------------------------------------------------------------------
## Written by record_scale_model() inside each DA script, so these rows always
## describe the run that produced the numbers above. A dataset missing here has
## simply not been re-run since the recording was added -- it is not a failure.
scale_records <- collect_scale_records(c(DIR_SELEX, DIR_ZEMB, DIR_STAMMLER,
                                         DIR_JIN, DIR_PRASSE, DIR_VANDEPUTTE))

if (!is.null(scale_records)) {
  ## The reportable view: the choice, both arms' counts, and how much the choice
  ## actually changed. n_sig_selected must equal whichever arm was chosen; a
  ## mismatch means the record and the fit came from different runs.
  scale_model_out <- scale_records[, c("dataset", "pair", "scale_model",
                                       "n_sig_selected", "n_sig_laplace",
                                       "n_sig_bootstrap", "n_shared",
                                       "only_laplace", "only_bootstrap",
                                       "stat_bootstrap", "stat_laplace",
                                       "n_laplace_fail", "D", "N", "S", "run_at")]
  scale_model_out <- scale_model_out[order(scale_model_out$dataset,
                                           scale_model_out$pair), ]
  print(scale_model_out, row.names = FALSE)

  ## A dataset with BOTH an unkeyed row (pair = NA) and keyed ones has a record
  ## left over from before that study's runs were keyed -- Prasse's two filters
  ## both wrote pair = NA until 2026-08-20, so the second run silently replaced
  ## the first. The stale row is reported here rather than dropped, because
  ## deleting a recorded run is not this script's call; it disappears on its own
  ## once each DA script has been re-run under its own key.
  keyed <- !is.na(scale_records$pair) & nzchar(scale_records$pair)
  stale <- !keyed & scale_records$dataset %in% scale_records$dataset[keyed]
  if (any(stale)) {
    warning("superseded scale-model record(s) with no filter/pair key: ",
            paste0(scale_records$dataset[stale], " (run ",
                   scale_records$run_at[stale], ")", collapse = ", "),
            " -- predates the keyed records and is NOT what this table reports. ",
            "Drop with: r <- readRDS(f); saveRDS(r[!(r$dataset == \"<name>\" & is.na(r$pair)), ], f)",
            call. = FALSE)
  }

  agree <- with(scale_records,
                n_sig_selected == ifelse(scale_model == "laplace",
                                         n_sig_laplace, n_sig_bootstrap))
  if (any(!agree, na.rm = TRUE)) {
    warning("selected-arm count disagrees with the per-arm count on: ",
            paste(scale_records$dataset[which(!agree)], collapse = ", "),
            " -- the record and the fit are out of sync; re-run those scripts.",
            call. = FALSE)
  }

  saveRDS(scale_model_out, file = file.path(ROOT, "real_data_scale_model.RData"))
  write.csv(scale_model_out, file = file.path(ROOT, "real_data_scale_model.csv"),
            row.names = FALSE)
}


#------------------------------------------------------------------------------------
# Scale model, FEATURE by feature: which taxa each arm flagged
#------------------------------------------------------------------------------------
## The section above says how many taxa each scale model called and how many they
## share. That is not enough to say anything about the features themselves: two
## arms can agree on a count and disagree on every taxon in it, and a TP/FP pair
## cannot tell you WHICH true positive an arm found.
##
## sparse_auto() now returns all three results -- hybrid, laplace-only,
## bootstrap-only -- and record_scale_model() writes the per-feature tables, so
## this reads them back and puts the identities next to the counts. Nothing is
## refitted here.
arm_features <- collect_arm_features(c(DIR_SELEX, DIR_ZEMB, DIR_STAMMLER,
                                       DIR_JIN, DIR_PRASSE, DIR_VANDEPUTTE))

## nrow() as well as NULL: collect_arm_features() returns a 0-row table when the
## records exist but no scale model called anything, and the annotation below
## needs at least one row to attach to.
if (!is.null(arm_features) && nrow(arm_features) > 0) {

  ## Each study's truth set, keyed exactly as the records are keyed, so a row can
  ## say TP or FP rather than only "called". A (dataset, pair) with no entry here
  ## keeps is_true = NA -- withheld, not false. That is deliberate for the runs
  ## this file does not score: Jin's ASV/OTU levels have no measured total
  ## abundance, and its two stratified locations have EMPTY gold references at
  ## 4 vs 4, so only the pooled subset named by JIN_PAIR is scored here.
  ## Stammler's fifteen pool pairs are scored in all_pairs_results.R.
  truth_sets <- list(
    `McMurrough2014|NA`  = rownames(selex_truth),
    `Zemb2020|NA`        = zemb_truth,
    `Stammler2016|Pool6` = stammler_truth,
    `Vandeputte2017|NA`  = vandeputte_truths)
  truth_sets[[paste0("Jin2022|",    JIN_PAIR)]]    <- jin_truth
  ## Keyed by the RECORD key, not the run name: collect_arm_features() keys its
  ## rows "<dataset>|<pair>" straight from sparse_arm_sig.rds, and the rawVSqseq
  ## passes are recorded there as "mean5_QSeq" / "mean5_raw". Keying this on
  ## PRASSE_PAIR would leave every Prasse feature unscored (is_true = NA) while
  ## looking like a deliberate withholding.
  truth_sets[[paste0("Prasse2015|", PRASSE_KEY)]] <- golden_ref_prasse$category

  ## Persisted so downstream files do not have to re-derive them. The truth sets
  ## are assembled here from six different study-specific sources; a table that
  ## reports "signal sparsity" as truths/D needs their sizes, and transcribing
  ## those by hand is how that column goes stale when a study is re-filtered.
  saveRDS(truth_sets, file = file.path(ROOT, "real_data_truth_sets.RData"))
  write.csv(data.frame(key = names(truth_sets),
                       n_truth = lengths(truth_sets),
                       row.names = NULL),
            file = file.path(ROOT, "real_data_truth_sizes.csv"), row.names = FALSE)

  key <- paste(arm_features$dataset,
               ifelse(is.na(arm_features$pair), "NA", arm_features$pair), sep = "|")
  arm_features$scored  <- key %in% names(truth_sets)
  arm_features$is_true <- mapply(function(k, f)
    if (!k %in% names(truth_sets)) NA else f %in% truth_sets[[k]],
    key, arm_features$feature, USE.NAMES = FALSE)

  ## Where the arms actually disagree, with the truth label attached. This is the
  ## table to read when asking whether the scale-model choice cost or bought
  ## anything: a laplace_only row with is_true TRUE is a true positive the
  ## bootstrap arm would have missed, and FALSE is a false positive it avoided.
  disagree <- arm_features[arm_features$set %in% c("laplace_only", "bootstrap_only"), ]

  cat("\n== how each scale model's calls split, per study ==\n")
  print(arm_feature_counts(arm_features), row.names = FALSE)

  cat("\n== features the two scale models DISAGREE on ==\n")
  if (nrow(disagree)) {
    print(disagree[, c("dataset", "pair", "feature", "set", "selected_arm",
                       "is_true", "mean_laplace", "q_laplace",
                       "mean_bootstrap", "q_bootstrap")], row.names = FALSE)
  } else {
    cat("  none -- the two arms called the same features in every recorded run\n")
  }

  ## TP/FP for each arm, recomputed from the identities. These must reproduce the
  ## sparse_SSRV_B / sparse_SSRV_L rows of the table above; they are printed
  ## because a mismatch is the cheapest signal that a study's two record files
  ## came from different runs.
  scored <- arm_features[arm_features$scored, ]
  if (nrow(scored)) {
    arm_score <- do.call(rbind, lapply(split(scored, paste(scored$dataset, scored$pair)),
      function(d) data.frame(
        dataset = d$dataset[1], pair = d$pair[1], selected_arm = d$selected_arm[1],
        TP_laplace   = sum(d$in_laplace   & d$is_true,  na.rm = TRUE),
        FP_laplace   = sum(d$in_laplace   & !d$is_true, na.rm = TRUE),
        TP_bootstrap = sum(d$in_bootstrap & d$is_true,  na.rm = TRUE),
        FP_bootstrap = sum(d$in_bootstrap & !d$is_true, na.rm = TRUE),
        stringsAsFactors = FALSE)))
    rownames(arm_score) <- NULL
    cat("\n== TP/FP per scale model, from the feature identities ==\n")
    print(arm_score, row.names = FALSE)
  } else {
    arm_score <- NULL
  }

  saveRDS(list(features = arm_features, disagree = disagree, score = arm_score),
          file = file.path(ROOT, "real_data_arm_features.RData"))
  write.csv(arm_features, file = file.path(ROOT, "real_data_arm_features.csv"),
            row.names = FALSE)
}

## The feature-level companion to real_scale_out: same runs, one row per taxon
## either scale model called, so the counts above can be traced to names.
if (file.exists(file.path(ROOT, "real_data_arm_features.RData"))) {
  real_arm_out <- readRDS(file = file.path(ROOT, "real_data_arm_features.RData"))
  View(real_arm_out$features)
}

readRDS(file = file.path(ROOT, "real_data_scale_model.RData"))

################################################################################
## Estimated theta^perp, against the external quantification --------------------
##
## The scale itself, not its consequences. theta^perp is the common log2 shift
## that turns this method's relative log fold changes into absolute ones, and
## under the model it IS the log2 fold change of the total load between the two
## groups. Three of these studies measured that total independently -- Zemb by
## qPCR, Jin by BarBIQ absolute counts, Vandeputte by flow cytometry -- so on
## those the estimate can be SCORED (bias, RMSE, interval coverage) and not only
## used. The other studies have no external total and appear here unscored, which
## is why `obs` and `covered` are NA for them rather than 0 and FALSE.
##
## Everything below is read back from the theta_perp_record.rds each DA script
## writes; nothing is refitted. The draws are on disk too
## (theta_perp_draws.rds), so coverage at a level other than 95% needs no refit
## -- see read_perp_draws().
##
## All three arms are reported. `hybrid` is the selection rule's answer and
## always duplicates whichever of laplace / bootstrap it chose; the other two are
## what the unselected scale model would have estimated.
################################################################################
perp_records <- collect_perp_records(c(DIR_SELEX, DIR_ZEMB, DIR_STAMMLER,
                                       DIR_JIN, DIR_PRASSE, DIR_VANDEPUTTE))

if (!is.null(perp_records) && nrow(perp_records) > 0) {

  cat("\n== estimated theta^perp (log2 fold change of TOTAL load, treatment - control) ==\n")
  perp_show <- c("dataset", "pair", "arm", "mean", "sd", "lower", "upper",
                 "fold_change", "obs", "bias", "covered")
  print(perp_records[, perp_show], row.names = FALSE, digits = 4)

  ## TWO references, and they are not interchangeable.
  ##
  ##   obs    what the assay measured: the whole community, before any filter.
  ##   obs_S  that measurement projected onto the features the model was fitted
  ##          on, which is the quantity theta^perp actually estimates -- the
  ##          Dirichlet draws are proportions renormalised within the kept set,
  ##          so the shift identified is Delta log2 T_S, not Delta log2 T.
  ##
  ## obs_S = obs + dlog2_r exactly. Both are reported because obs is the physical
  ## number and obs_S is the like-for-like one, and because how much of the
  ## community counts as "outside S" is a judgement (see perp_matched_obs()):
  ## quoting only the corrected value would bury that judgement in one column.
  if (!is.null(perp_records$obs_S) && any(is.finite(perp_records$obs_S))) {
    cat("\n== the same runs against the FILTER-MATCHED reference ==\n")
    print(perp_records[, c("dataset", "pair", "arm", "D_all", "D_kept",
                           "r_min", "r_max", "dlog2_r",
                           "obs_S", "bias_S", "covered_S")],
          row.names = FALSE, digits = 4)
    n_um <- sum(!is.finite(perp_records$obs_S) & is.finite(perp_records$obs))
    if (n_um > 0)
      cat("   (", n_um, " scored row(s) carry no matched reference -- their DA ",
          "script has not been re-run with `matched =`)\n", sep = "")
  }

  ## The three numbers the claim "we recover the absolute shift" has to be
  ## defended with. Computed over the recorded runs that carry a measurement;
  ## `n` says how many those are, which is the context `coverage` is meaningless
  ## without -- a coverage of 1.00 over n = 1 is one interval, not a result.
  perp_score <- score_theta_perp(perp_records)
  if (!is.null(perp_score)) {
    cat("\n== theta^perp accuracy against the external measurement, by scale model ==\n")
    cat("   bias/rmse/coverage are against obs (whole community);\n")
    cat("   the _S columns are against obs_S (matched to each study's filter),\n")
    cat("   over the n_S runs that carry one.\n")
    print(perp_score, row.names = FALSE, digits = 4)
    if (!is.na(attr(perp_score, "n_unscored")) && attr(perp_score, "n_unscored") > 0)
      cat("   (", attr(perp_score, "n_unscored"),
          " rows have no external measurement and are excluded)\n", sep = "")

    cat("\n== the same, per study ==\n")
    print(score_theta_perp(perp_records, by = c("dataset", "arm")),
          row.names = FALSE, digits = 4)
  }

  saveRDS(list(records = perp_records, score = perp_score),
          file = file.path(ROOT, "real_data_theta_perp.RData"))
  write.csv(perp_records, file = file.path(ROOT, "real_data_theta_perp.csv"),
            row.names = FALSE)
}
