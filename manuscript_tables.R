################################################################################
## Manuscript tables -- the whole bundle, in one file.
##
## Running this file prints, and writes, the LaTeX for every table in the
## manuscript and its supplement. Nothing is transcribed by hand: each number is
## read from a saved analysis artifact, so re-running the analyses and then
## re-running this file brings the manuscript back in sync.
##
## This is now an AGGREGATOR. T1, T2, T3 and the S-series are built by the task
## scripts, which it sources; only the Jin tables and T5 are built here, because
## they belong to the Jin objects and the theta^perp record this file loads. To
## work on one of those tables, edit its task script -- not this file.
##
## Emits, in order:
## Emits, in bundle order -- the S-number each one carries is in the left column:
##
##   S1   T1     Summary of all datasets.                                [Task 1]
##   S2   T2     DA results, Sparse SSRV under the hybrid rule at the
##               default bandwidth.                                     [Task 2]
##   S3   S1     All nine Sparse SSRV variants on the scored datasets.   [Task 3]
##   S4   T5     Recovery of the measured change in total load.
##   S5   T4.1   Jin: dimensions and per-scale-model discovery counts.
##   S6   T4.2   Jin: discovery counts against the published analysis.
##   S7   T4.3   Jin: cOTU-driven cross-representation table.
##   S8   T4.4   Jin: feature-driven cross-representation table.
##   S9   T4.5   Jin: taxonomy of every ASV and OTU called by either analysis.
##   S10  S3     The same nine settings on Jin's SUBSETS, discovery counts. [Task 3]
##
## RETIRED: the two-arm results table (its arms are columns of S3 now) and the
## fifteen-pool-pair variant table (the Task 5 figure shows the same thing).
##
## FIGURES are produced by "Task 4 ...", "Task 5 ..." and "Task 6 ...", and by
## the simulation code. Only tables are assembled here.
##
## SOURCES
##   real_data_scale_model.RData    dimensions, selected scale model, per-arm
##                                  discovery counts, one row per (dataset, pair).
##   real_data_results.RData        the FP/(TP) matrix, datasets x methods, that
##                                  `real data results.R` builds and saves.
##   real_data_arm_features.RData   per-arm TP/FP at feature level ($score), used
##                                  to sanity-check the declared truth sizes.
##   jin_DA_biological_interpretation.R
##                                  sourced, not copied. It owns the Jin objects
##                                  and every LaTeX helper, and the Jin tables
##                                  below are built with its functions so the two
##                                  files cannot drift apart. Sourcing it also
##                                  refreshes the Jin-only outputs in
##                                  2022_Jin.../tables, which is intended.
##
##   task_common.R                sourced by every task script: the paths, the
##                                dataset registry, the truth sets, the LaTeX
##                                helpers, the consistency checks, and the reader
##                                that pulls all nine Sparse SSRV variants off a
##                                saved fit.
##
## DEPENDENCY. T2, T3 and S1 rest on `real_data_results.RData` and on the per-arm
## records, which only `real data results.R` and the DA scripts write. Nothing
## here refits anything, so if those have not been re-run since the last analysis
## change these tables describe an earlier fit. check_results_freshness()
## compares timestamps and warns, and Task 3 additionally checks its variant grid
## against T2 cell by cell.
################################################################################

################################################################################
## PATHS -- nothing below this block hardcodes an absolute path.
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

DIR_JIN_ <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates")
DIR_OUT  <- file.path(ROOT, "manuscript_tables")
dir.create(DIR_OUT, showWarnings = FALSE)

## Sourcing the Jin file brings in ROOT, DIR_JIN, DATESTAMP, the tidy Jin frames
## (ssrv, jin, runs), every tex_* helper, and the two cross-representation
## builders. `jin.quiet` stops it echoing its own LaTeX over ours.
options(jin.quiet = TRUE)
source(file.path(DIR_JIN_, "jin_DA_biological_interpretation.R"))
message("")

################################################################################
## THE CROSS-DATASET TABLES ARE BUILT BY THE TASK SCRIPTS
##
## T1, T2 and T3 used to be defined here as well as being the deliverable of
## "Task 1 ..." and "Task 2 ...". Two definitions of one table is how a caption
## gets fixed in one place and printed from the other, so this file now SOURCES
## those scripts and collects what they build. Each writes its own .tex on the
## way through, which is what it would do run on its own; `task.quiet` stops them
## echoing their LaTeX into the middle of this file's bundle.
##
## Sourcing Task 3 also refreshes real_data_ssrv_variants.RData and emits the
## three supplementary variant tables, which belong in the same bundle.
##
## task_common.R -- sourced by each of them -- provides the paths, the dataset
## registry, the truth sets, the LaTeX helpers and the consistency checks that
## this file used to declare inline. It is loaded here too, so anything below
## that reaches for scale_records, fptp, truth_sets or DATASET_META finds them
## whether or not a task script has run yet.
################################################################################
source(file.path(ROOT, "task_common.R"))

options(task.quiet = TRUE)
source(file.path(ROOT, "Task 1 dataset summary table.R"))   # T1
source(file.path(ROOT, "Task 2 DA results table.R"))        # T2
source(file.path(ROOT, "Task 3 sparse SSRV variant tables.R"))  # S1, S3
options(task.quiet = FALSE)

################################################################################
## JIN-ONLY CHECKS
##
## The cross-dataset checks -- truth-set sizes, hybrid-versus-arm agreement, the
## true-positive bound -- run inside "Task 2 DA results table.R", on the same
## numbers it prints. This one is Jin's and belongs here, because `runs` comes
## from the Jin file.
##
## Scope note: `runs` comes from build_runs() in jin_DA_biological_interpretation.R,
## which enumerates LEVEL_KEYS x LOC_KEYS. LOC_KEYS is prox/dist, so every Jin
## table built in this file is the LOCATION-STRATIFIED (supplementary) analysis and
## is unaffected by which subset the main DA results table reports. The pooled
## subset is the primary result and lives in `real data results.R` (JIN_PAIR).
##
## One selection rule should be in force across every Jin subset. The rule became
## "width" -- the larger median 95% credible interval width -- when the simulation
## study moved to it; a subset still carrying an older rule was fitted before that
## change and is not comparable with the rest.
################################################################################
check_scale_rule <- function() {
  r <- unique(runs$scale_rule)
  if (length(r) > 1) {
    warning("Jin subsets were fitted under different scale-selection rules (",
            paste(sort(r), collapse = ", "),
            "). Re-fit the stragglers, and restore the scale-rule column in ",
            "Table~\\ref{tab:jin_dimensions} until they agree.",
            call. = FALSE, immediate. = TRUE)
  } else {
    message("Scale-selection rule in force across every Jin subset: ", r)
  }
  invisible(r)
}
check_scale_rule()


################################################################################
## T4.1 -- Jin dimensions and per-scale-model discovery counts
##
## The pair key is split into the two things it encodes, representation and
## location, because that is how the section discusses the design. The scale-rule
## column is dropped: one rule is in force across every subset, so the column
## repeats a constant. If a future run ever mixes rules, put it back --
## check_scale_rule() below says when that happens.
################################################################################
tab_jin_dimensions <- function() {
  d <- runs[order(match(runs$level, LEVEL_KEYS), match(runs$location, LOC_KEYS)), ]
  rows <- character(0); seen <- character(0)
  for (i in seq_len(nrow(d))) {
    if (d$level[i] %in% seen) lv <- "" else {
      lv <- paste0("\\multirow{", sum(d$level == d$level[i]), "}{*}{", d$level[i], "}")
      if (length(seen)) rows <- c(rows, "        \\addlinespace[2pt]")
    }
    seen <- c(seen, d$level[i])
    rows <- c(rows, row_tex(lv, LOC_LABEL[[d$location[i]]],
                            sprintf("$%d \\times %d$", d$D[i], d$N[i]),
                            d$scale_model[i], d$n_selected[i],
                            d$n_bootstrap[i], d$n_laplace[i], d$n_shared[i]))
  }
  c("\\begin{table}[!tbp]",
    "    \\centering",
    "    \\resizebox{\\textwidth}{!}{%",
    "    \\begin{tabular}{llccccc c}",
    "        \\toprule",
    "        \\shortstack{Feature\\\\representation} & Location & \\shortstack{Dimension\\\\($D \\times N$)} & \\shortstack{Selected\\\\scale model} &",
    "        $n_{\\text{sig, selected}}$ & $n_{\\text{sig, bootstrap}}$ & $n_{\\text{sig, Laplace}}$ & $n_{\\text{shared}}$ \\\\",
    "        \\midrule",
    rows,
    "        \\bottomrule",
    "    \\end{tabular}%",
    "    }",
    paste0("    \\caption{Dimensions of each \\citet{jin2022quantsinglecell} subset and the number of ",
           "features called under each scale model. $D$ is the number of features retained after ",
           "filtering and $N$ the number of samples. The selected column names the arm the selection ",
           "rule chose for that subset; the last four columns count the features called by the selected ",
           "arm, by the bootstrap arm, by the Laplace arm, and by both.}"),
    "    \\label{tab:jin_dimensions}",
    "\\end{table}")
}

################################################################################
## T4.5 -- taxonomy of every ASV and OTU called by either analysis
##
## The layout of the old "without a significant cOTU" supplement tables, but over
## the full set: every ASV and OTU either analysis called, mapped or unmapped,
## with the same * / dagger markers as the feature-driven cross-representation
## table so the two can be read together.
################################################################################
tab_jin_sig_feature_taxonomy <- function(ssrv_tbl) {
  want <- sig_features_either(ssrv_tbl)

  cotu_cell <- function(cotu, loc) {
    if (is.na(cotu) || !nzchar(cotu) || cotu == "Not detected") return(tex_tax(cotu))
    row <- ssrv_tbl[ssrv_tbl$level == "cOTU" & ssrv_tbl$location == loc &
                    ssrv_tbl$feature == cotu, ]
    ## A mapped cOTU need not be among the cOTUs analysed at that location --
    ## the mapping points into the full BarBIQ catalogue, prevalence filtering
    ## trims it -- and an unanalysed cOTU cannot carry a call either way.
    if (!nrow(row)) return(tex_escape(cotu))
    paste0(tex_escape(cotu), sig_marks("cOTU", loc, cotu, row$sig[1]))
  }
  rows <- character(0)
  for (loc in LOC_KEYS) for (lv in c("ASV", "OTU")) {
    f <- want$feature[want$location == loc & want$level == lv]
    if (!length(f)) next
    d <- ssrv_tbl[ssrv_tbl$level == lv & ssrv_tbl$location == loc &
                  ssrv_tbl$feature %in% f, ]
    ## Alphabetical by genus within each representation, unclassified features
    ## last, and by identifier inside a genus. Grouping the lineages together is
    ## the point of the table; identifier order scatters them.
    d <- d[order(is.na(d$Genus), d$Genus, nat_num(d$feature), d$feature), ]
    if (length(rows)) rows <- c(rows, "        \\addlinespace[2pt]")
    rows <- c(rows, vapply(seq_len(nrow(d)), function(i) row_tex(
      LOC_LABEL[[loc]], lv,
      paste0(tex_escape(d$feature[i]),
             sig_marks(lv, loc, d$feature[i], d$sig[i], d$cOTU_override[i])),
      ## Which analyses called the MAPPED cOTU, so the feature can be read
      ## against its cOTU without turning back to the cross-representation
      ## tables. Blank when the mapping is absent, or when the cOTU is not one
      ## of the cOTUs analysed at that location.
      cotu_cell(d$cOTU[i], loc),
      tex_tax(d$Kingdom[i]), tex_tax(d$Phylum[i]), tex_tax(d$Class[i]),
      tex_tax(d$Order[i]),   tex_tax(d$Family[i]), tex_genus(d$Genus[i]),
      tex_num(d$lfc[i], DIG_LFC_LONG)), character(1)))
  }
  note <- paste0("\\scriptsize \\textsuperscript{*}Significant in \\citet{jin2022quantsinglecell} ",
                 "(FDR $<", JIN_FDR_MAX, "$, absolute fold change $>", 2^JIN_ABS_LOG2FC_MIN,
                 "$). \\textsuperscript{\\dag}Significant in ", tolower(ARM), " Sparse SSRV. ",
                 "Markers on a feature say which analysis called that feature; markers on a mapped ",
                 "cOTU say which analysis called the cOTU itself, and a mapped cOTU with no marker was ",
                 "called by neither, or was not among the cOTUs analysed at that location. ",
                 "``Not detected'' denotes a feature the supplied mapping did not assign to a cOTU; ",
                 "dashes mark unavailable taxonomic ranks. ",
                 "\\textsuperscript{a}The saved mapping object labels ",
                 paste(tex_escape(names(MANUAL_COTU_MAP)), collapse = " and "),
                 " as ``Not detected''; it is mapped here because \\citet{jin2022quantsinglecell} ",
                 "identify it as an OTU representation of ",
                 paste(tex_escape(unname(MANUAL_COTU_MAP)), collapse = " and "), ".")
  head_block <- c(
    "    \\toprule",
    "    Location & \\shortstack{Feature\\\\representation} & Feature & \\shortstack{Mapped\\\\cOTU} & Kingdom & Phylum & Class & Order & Family & Genus & \\shortstack{Estimated\\\\$\\log_2$ LFC} \\\\",
    "    \\midrule")
  c("\\clearpage",
    "\\begin{landscape}",
    "{\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\setlength{\\LTcapwidth}{\\linewidth}",
    "\\begin{longtable}{@{}llllllllllr@{}}",
    paste0("    \\caption{Taxonomy of every ASV and OTU called by either analysis in the ",
           "\\citet{jin2022quantsinglecell} data, with its mapped cOTU where the supplied mapping ",
           "provides one. Features are ordered by location and representation, and by genus within each, ",
           "with unclassified features last. Estimated ",
           "LFCs are ", tolower(ARM), " Sparse SSRV posterior means for the vitamin A-deficient versus ",
           "vitamin A-sufficient contrast on the $\\log_2$ scale; the markers, which follow ",
           "Table~\\ref{tab:jin_cross_representation_sig_features}, say which analysis called each ",
           "feature. Both mapped and unmapped features are listed.}"),
    "    \\label{tab:jin_sig_feature_taxonomy}\\\\",
    head_block,
    "    \\endfirsthead",
    "    \\multicolumn{11}{l}{\\tablename~\\thetable\\ continued}\\\\",
    head_block,
    "    \\endhead",
    "    \\midrule",
    "    \\multicolumn{11}{r}{Continued on next page}\\\\",
    "    \\endfoot",
    "    \\bottomrule",
    "    \\endlastfoot",
    rows,
    "\\end{longtable}",
    "\\vspace{4pt}\\noindent",
    paste0("    ", note),
    "\\par",
    "}",
    "\\end{landscape}")
}

################################################################################
## T5 -- does the scale model recover the measured change in total load?
##
## theta^perp is the log2 fold change of the TOTAL load between conditions: the
## one quantity sequence counts cannot identify and the scale model has to
## supply. Three studies measured it directly -- Zemb by qPCR, Vandeputte by flow
## cytometry, Jin by BarBIQ absolute counts -- so for those the posterior can be
## held against an observed value.
##
## The summaries come from theta_perp_check.csv, the record the check writes --
## NOT from the draws file beside it. The two disagree for the Laplace arm by
## about a percent, and it is not a provenance problem: one call to
## record_theta_perp() writes both, from one fit. They are two estimators of the
## same posterior.
##
##   Laplace   theta^perp is a Gaussian mixture, one component per posterior
##             draw, and the record summarises it ANALYTICALLY -- mean(center)
##             and the law of total variance over (center, var), with mixture
##             quantiles. The draws are one Monte Carlo realisation of that same
##             mixture, so they carry sampling noise the record does not.
##   Bootstrap no mixture to exploit, so the record IS the summary of the draws
##             and the two agree to floating point.
##
## The analytic summary is the better estimate of the two, so it is what the
## table reports. check_theta_perp_agreement() verifies the relationship holds:
## bootstrap must match the draws almost exactly, and Laplace must sit within
## Monte Carlo error of them.
##
##   bias     posterior mean minus the observed value.
##   RMSE     sqrt(E[(theta^perp - observed)^2]) = sqrt(bias^2 + sd^2). An exact
##            identity, so with the analytic mean and sd this is exact too --
##            which computing it from the draws would not be.
##   coverage whether the 95% credible interval contains the observed value. One
##            subset contributes one interval, so this is an indicator, not a
##            rate over replicates.
################################################################################

F_PERP_CSV   <- file.path(ROOT, "theta_perp_check", "theta_perp_check.csv")
F_PERP_DRAWS <- file.path(ROOT, "theta_perp_check", "theta_perp_draws.rds")

## Arms reported, in row order. "hybrid" is not a scale model -- it duplicates
## whichever of these two the selection rule chose -- so it is marked on the
## selected row instead of repeated as a third one. Add it here for its own row.
THETA_ARMS  <- c("bootstrap", "laplace")
THETA_LABEL <- c(bootstrap = "Bootstrap", laplace = "Laplace")

## Jin measured total abundance at cOTU level only, so its ASV and OTU subsets
## carry no observed value and drop out of the is.na(obs) filter.
theta_perp_summaries <- function() {
  if (!file.exists(F_PERP_CSV)) {
    stop("missing theta^perp record:\n  ", F_PERP_CSV,
         "\nRun theta_perp_check.R, then re-run this file.")
  }
  d <- read.csv(F_PERP_CSV, stringsAsFactors = FALSE)
  d <- d[d$arm %in% THETA_ARMS & !is.na(d$obs), , drop = FALSE]
  if (!nrow(d)) stop("no theta^perp record carries an observed value")
  d$pair     <- ifelse(is.na(d$pair) | !nzchar(d$pair), NA_character_, d$pair)
  d$selected <- d$arm == d$selected
  d$bias     <- d$mean - d$obs
  ## Exact given the recorded mean and sd; see the note above.
  d$rmse     <- sqrt(d$bias^2 + d$sd^2)
  d$covered  <- d$obs >= d$lower & d$obs <= d$upper

  ## The same three against the FILTER-MATCHED reference, recomputed here from
  ## the recorded obs_S for exactly the reason bias/rmse/covered are recomputed
  ## from obs: the manuscript's numbers should follow from the columns the table
  ## prints, not from a summary written by another script.
  ##
  ## obs_S is `obs` projected onto the features the model was fitted on -- see
  ## perp_matched_obs() in scale_model_record.R. Both are carried and neither
  ## replaces the other: obs is the assay's own value over the whole community,
  ## obs_S is the like-for-like comparison, and obs_S = obs + dlog2_r exactly.
  ## Which of the two a table should lead with is a presentational decision this
  ## function deliberately does not make; it supplies both.
  ##
  ## NA for any row recorded before obs_S existed, so an old CSV still summarises
  ## instead of erroring -- but a table built on the _S columns will then be
  ## empty rather than silently falling back to the uncorrected value.
  if (is.null(d$obs_S)) d$obs_S <- NA_real_
  d$bias_S    <- d$mean - d$obs_S
  d$rmse_S    <- sqrt(d$bias_S^2 + d$sd^2)
  d$covered_S <- d$obs_S >= d$lower & d$obs_S <= d$upper
  rownames(d) <- NULL
  d
}

## Jin reports a replicate standard deviation for every scale measurement, so
## its observed theta^perp comes with a measurement interval. No other study in
## the record does, and those rows keep a bare point estimate.
source(file.path(ROOT, "jin_scale_uncertainty.R"))
theta_perp <- theta_perp_summaries()
obs_ci <- jin_scale_uncertainty(DIR_JIN_)
if (!is.null(obs_ci)) {
  i <- match(paste(theta_perp$dataset, theta_perp$pair),
             paste(obs_ci$dataset, obs_ci$pair))
  theta_perp$obs_lower <- obs_ci$lower[i]
  theta_perp$obs_upper <- obs_ci$upper[i]
  ## The record and the sheet must agree on the point estimate; if they do not,
  ## the interval belongs to a different number from the one being reported.
  d <- abs(theta_perp$obs - obs_ci$obs[i])
  if (any(d > 1e-6, na.rm = TRUE)) {
    warning("Jin's observed theta^perp differs between the theta^perp record and ",
            "Supplementary Data 8 by up to ", signif(max(d, na.rm = TRUE), 3),
            call. = FALSE, immediate. = TRUE)
  }
} else {
  theta_perp$obs_lower <- NA_real_; theta_perp$obs_upper <- NA_real_
}

## Confirms the record and the draws stand in the expected relationship. A
## bootstrap arm that does not match its draws, or a Laplace arm far outside
## Monte Carlo error of them, means something other than sampling noise.
check_theta_perp_agreement <- function(d, tol_boot = 1e-6, tol_mix_se = 6) {
  if (!file.exists(F_PERP_DRAWS)) return(invisible(NULL))
  draws <- readRDS(F_PERP_DRAWS)
  bad <- character(0)
  for (i in seq_len(nrow(d))) {
    e <- draws[[paste0(d$dataset[i], "|", if (is.na(d$pair[i])) "NA" else d$pair[i])]]
    if (is.null(e)) next
    x <- e[[d$arm[i]]]$draws
    if (is.null(x) || !length(x)) next
    gap <- abs(mean(x) - d$mean[i])
    if (d$arm[i] == "bootstrap") {
      if (gap > tol_boot)
        bad <- c(bad, sprintf("%s/%s bootstrap: record %.6f vs draws %.6f",
                              d$dataset[i], d$pair[i], d$mean[i], mean(x)))
    } else {
      ## The draws are a sample of the mixture the record describes exactly, so
      ## the gap should sit inside a few standard errors of that sample mean.
      se <- stats::sd(x) / sqrt(length(x))
      if (gap > tol_mix_se * se)
        bad <- c(bad, sprintf("%s/%s laplace: record %.4f vs draws %.4f (%.1f SE)",
                              d$dataset[i], d$pair[i], d$mean[i], mean(x), gap / se))
    }
  }
  if (length(bad)) {
    warning("theta^perp record and draws disagree beyond Monte Carlo error:\n  ",
            paste(bad, collapse = "\n  "),
            "\nRe-run theta_perp_check.R for those studies.",
            call. = FALSE, immediate. = TRUE)
  } else {
    message("theta^perp: record matches the draws (bootstrap exactly, Laplace ",
            "within Monte Carlo error of the mixture it summarises)")
  }
  invisible(bad)
}
check_theta_perp_agreement(theta_perp)

## THREE VARIANTS of one table, so the shared parts -- row membership, the
## observed column, the arithmetic, the captions' common half -- cannot drift:
##
##   "main"      the five studies with a measurement, one row each, the arm the
##               selection rule chose. Stammler on its headline Pool1-vs-Pool6
##               comparison. No scale-model column: every row is the hybrid
##               result, so naming the arm each time adds nothing the
##               supplementary table does not say better.
##   "stammler"  the same shape, restricted to Stammler and expanded to all
##               fifteen pool pairs -- the theta^perp companion to the all-pairs
##               TP/FP figure in the supplement.
##   "arms"      every run in both arms, the chosen one daggered. Ordered with
##               the five headline runs first and Stammler's remaining pairs
##               after, so the top of it lines up row-for-row with "main".
tab_theta_perp_recovery <- function(d, variant = c("main", "stammler", "arms")) {
  variant          <- match.arg(variant)
  hybrid_only      <- variant %in% c("main", "stammler")
  show_comparison  <- variant %in% c("stammler", "arms")
  show_scale_model <- variant == "arms"

  ## Membership comes from DATASET_META, not from a list kept here, so this table
  ## and T2 cannot describe different fits of the same study. The theta^perp
  ## record accumulates a row for every subset ever fitted and never forgets one,
  ## so without this a superseded filter keeps reporting itself: Prasse's prev10
  ## runs (D = 181) were still rendering after the analysis moved to mean5
  ## (D = 85), giving that study four rows across two filters while every other
  ## table showed it once.
  ##
  ## Stammler is the exception, and deliberately so: it is a dilution series, so
  ## every pair of its six pools is a comparison with a known load ratio. The DA
  ## script records the headline Pool1-vs-Pool6 comparison TWICE, once under the
  ## bare tag "Pool6" and once as "Pool1_vs_Pool6". Which of the two is used
  ## depends on the variant, but never both -- listing both would put an
  ## identical row in the table twice with nothing to say they are one fit.
  reported_pair <- vapply(d$dataset, function(k) {
    p <- if (is.null(DATASET_META[[k]])) NA_character_ else DATASET_META[[k]]$pair
    if (is.null(p) || is.na(p)) "" else as.character(p)
  }, character(1))
  pair_key <- ifelse(is.na(d$pair), "", d$pair)
  stam     <- d$dataset == "Stammler2016"
  is_pair  <- grepl("_vs_", pair_key)

  keep <- switch(variant,
    ## Stammler via its "Pool6" tag, matching DATASET_META and T1/T2.
    main     = pair_key == reported_pair,
    ## Stammler only, and every pool pair of it.
    stammler = stam & is_pair,
    ## The headline run of each study -- Stammler's taken as Pool1_vs_Pool6 so it
    ## sits in the same tag space as the pairs listed after it -- plus the rest
    ## of the pairs.
    arms     = (!stam & pair_key == reported_pair) | (stam & is_pair))

  if (any(!keep)) {
    message("  theta^perp [", variant, "]: dropping ", sum(!keep),
            " row(s) this variant does not report -- ",
            paste(sort(unique(paste0(d$dataset[!keep], "/", pair_key[!keep]))),
                  collapse = ", "))
    d <- d[keep, , drop = FALSE]
    pair_key <- pair_key[keep]; stam <- stam[keep]; is_pair <- is_pair[keep]
    reported_pair <- reported_pair[keep]
  }
  if (!nrow(d)) stop("no theta^perp row survives the '", variant, "' selection")

  ## ORDER. "main" and "stammler" have one run per line and sort naturally.
  ## "arms" is asked to lead with the five headline runs, in dataset order, and
  ## follow them with Stammler's remaining pairs -- so the reader can read the
  ## first five blocks against the main table and the rest as the extra grid.
  headline <- if (variant == "arms")
                ifelse(stam, pair_key == "Pool1_vs_Pool6", TRUE) else rep(TRUE, nrow(d))
  d <- d[order(!headline,
               match(d$dataset, names(DATASET_META)),
               pair_key,
               match(d$arm, THETA_ARMS)), ]
  ## Re-derived from the reordered frame; the old vectors are in the pre-sort order.
  pair_key <- ifelse(is.na(d$pair), "", d$pair)

  ## The comparison column names Stammler's pool pairs. The other studies
  ## contribute one run each and the caption says which, so their cells are a
  ## dash rather than a restatement.
  comparison_label <- function(pair) {
    if (is.na(pair) || !nzchar(pair) || !grepl("_vs_", pair)) return(DASH)
    gsub("_vs_", " vs ", pair)
  }

  ## Dataset, subset and the observed value are properties of the study, not of
  ## a scale model, so they are written once per dataset and left blank on the
  ## second arm's row.
  ##
  ## `first` compares against the PREVIOUS row, not against every key seen so
  ## far. The difference only shows when a subset's rows are not contiguous, and
  ## then it is the difference between a visible fault and a silent one: the
  ## seen-set version blanks the label on a row that has drifted away from its
  ## group, so the row is read under whatever heading precedes it. Comparing with
  ## the previous row re-prints the heading instead, which is wrong-looking
  ## rather than wrong -- and a repeated heading is a bug report you can see.
  ## The sort above should make this unreachable; it is here so that the next
  ## unlisted subset name is caught by the eye rather than by nobody.
  ## MAIN-TEXT VERSION: one row per dataset, the arm the selection rule chose --
  ## which IS the hybrid result, since the hybrid posterior is that arm's, not a
  ## blend of the two. Taken by subsetting the SAME rows the full table renders,
  ## rather than reading the CSV's own "hybrid" rows, so the two tables cannot
  ## report different numbers for the same fit: the main table is exactly the
  ## daggered rows of the supplementary one.
  if (hybrid_only) {
    d <- d[d$selected, , drop = FALSE]
    if (!nrow(d)) stop("no theta^perp row is marked as the selected arm")
    ## RE-DERIVE, because this subset halves `d` and every vector computed from
    ## the pre-subset frame is now the wrong length. pair_key is read by index in
    ## the row loop below, so leaving it stale silently pairs each row's data
    ## with ANOTHER row's label: the Stammler table listed "Pool1 vs Pool2"
    ## twice, the second time carrying Pool1-vs-Pool3's numbers. Nothing in the
    ## output looked wrong -- both are real comparisons with plausible values.
    pair_key <- ifelse(is.na(d$pair), "", d$pair)
  }
  ## One row per run in these variants, so a repeated (dataset, comparison) pair
  ## means two arms survived the filter and the table would be claiming two
  ## different results for one fit.
  if (hybrid_only && anyDuplicated(paste(d$dataset, d$pair)))
    stop("theta^perp [", variant, "]: more than one arm survives for ",
         paste(unique(paste(d$dataset, d$pair)[duplicated(paste(d$dataset, d$pair))]),
               collapse = ", "))

  ## Does any reported run actually carry a measurement interval? The column
  ## subtitle and the caption sentence describing it are both conditional on
  ## this, so neither can promise an interval the table does not show.
  has_obs_ci <- any(!is.na(d$obs_lower))

  keys <- paste(d$dataset, d$pair)
  rows <- character(0)
  for (i in seq_len(nrow(d))) {
    first <- i == 1L || keys[i] != keys[i - 1L]
    ## In the ALL-ARMS table each arm occupies two lines, the estimate over its
    ## interval, and without a gap between arms the interval sits closer to the
    ## NEXT arm's estimate than to its own and reads as belonging to it.
    ##
    ## The one-row-per-run variants need none of that: there is no second arm to
    ## separate, so the gap only spread the table out. They set \arraystretch
    ## instead, exactly as the main results table does, which is what makes their
    ## line spacing match it.
    if (i > 1L)
      rows <- c(rows, sprintf("        \\addlinespace[%dpt]",
                              if (!show_scale_model) 3L else if (first) 7L else 4L))
    rows <- c(rows, row_tex(
      if (first) DATASET_META[[d$dataset[i]]]$label else "",
      if (!show_comparison) NULL
      else if (first) comparison_label(pair_key[i]) else "",
      ## The scale-model column exists only where the arms are shown separately.
      ## In the one-row-per-run variants every row IS the chosen arm, so naming
      ## it on each line marks nothing; the dagger goes with it for the same
      ## reason.
      if (!show_scale_model) NULL
      else paste0(THETA_LABEL[[d$arm[i]]],
                  if (d$selected[i]) "\\textsuperscript{\\dag}" else ""),
      ## The measurement interval sits under the observed value where the study
      ## reports the replicate spread to compute one.
      if (!first) "" else if (is.na(d$obs_lower[i])) tex_num(d$obs[i], DIG_LFC_SHORT)
      else tex_block(paste0(tex_num(d$obs[i], DIG_LFC_SHORT), " \\\\ ",
                            "{\\scriptsize ", tex_interval(d$obs_lower[i], d$obs_upper[i]), "}"),
                     align = "c"),
      ## The interval sits BENEATH the estimate, pulled up tight against it.
      ##
      ## The negative leading is what makes that readable. With the inter-row gap
      ## removed to tighten the table, a plain \\ put the interval equidistant
      ## between its own estimate and the NEXT row's -- 16pt against 17pt in the
      ## rendered table -- so it could be read as belonging to either. Rather
      ## than restore the gap, which is what made the table tall in the first
      ## place, the WITHIN-cell leading is cut: the interval now sits closer to
      ## the estimate it belongs to than to anything else, at no cost in height.
      tex_block(paste0(tex_num(d$mean[i], DIG_LFC_SHORT), " \\\\[-4pt] ",
                       "{\\scriptsize ", tex_interval(d$lower[i], d$upper[i]), "}"),
                align = "c"),
      tex_num(d$bias[i], DIG_LFC_SHORT),
      tex_num(d$rmse[i], DIG_LFC_SHORT),
      if (d$covered[i]) "Yes" else "No"))
  }

  c("\\begin{table}[!tbp]",
    "    \\centering",
    ## 1.2, the value the main results table uses, so the two read at the same
    ## density -- and set OUTSIDE \resizebox, where that table sets it. Applied
    ## only where the inter-row gap has been dropped; the all-arms table keeps
    ## that gap and needs no extra stretch on top of it.
    if (!show_scale_model) "    \\renewcommand{\\arraystretch}{1.2}" else NULL,
    "    \\resizebox{\\textwidth}{!}{%",
    paste0("    \\begin{tabular}{l", if (show_comparison) "l" else "",
           if (show_scale_model) "l" else "", "ccccc}"),
    "        \\toprule",
    ## The "(95% measurement CI)" subtitle appears only if a row actually carries
    ## one. Only Jin reports the replicate spread needed to compute it, and only
    ## for its STRATIFIED locations -- which this table no longer shows, so the
    ## subtitle was promising an interval that appears nowhere beneath it. Driven
    ## off the data rather than deleted, so it returns by itself if a reported run
    ## ever carries the interval again.
    paste0("        Dataset & ", if (show_comparison) "Comparison & " else "",
           if (show_scale_model) "Scale model & " else "",
           if (has_obs_ci)
             "\\shortstack{Observed $\\theta^{\\perp}$\\\\(95\\% measurement CI)}"
           else "Observed $\\theta^{\\perp}$",
           " & \\shortstack{Posterior $\\theta^{\\perp}$\\\\(95\\% CrI)} & Bias & RMSE & \\shortstack{Interval\\\\covers} \\\\"),
    "        \\midrule",
    rows,
    "        \\bottomrule",
    "    \\end{tabular}%",
    "    }",
    paste0("    \\caption{Recovery of the change in total microbial load. $\\theta^{\\perp}$ is the ",
           "$\\log_2$ fold change of the total load between conditions: the component sequence counts ",
           "cannot identify and the scale model must supply. The observed column comes from the external ",
           "measurement each study reports and is a property of the study rather than of a scale model, ",
           "so it is given once per dataset. ",
           ## Conditional for the same reason as the column subtitle: this
           ## described an interval that only Jin's stratified locations carry,
           ## and those are no longer among the reported runs.
           if (has_obs_ci)
             paste0("Where a study reports the replicate spread of its scale ",
                    "measurement -- \\citet{jin2022quantsinglecell} alone does, over four technical ",
                    "replicates per sample -- a 95\\% interval for the observed value is given beneath ",
                    "it, obtained by propagating each sample's standard error through the contrast. ",
                    "That interval covers assay noise only, not variation between mice. ")
           else "",
           "Bias is the posterior mean minus the observed value and ",
           "RMSE is $\\sqrt{\\mathrm{E}[(\\theta^{\\perp}-\\text{observed})^{2}]}=",
           "\\sqrt{\\text{bias}^{2}+\\mathrm{sd}^{2}}$, evaluated at each arm's recorded posterior ",
           ## COUNTED, not written out. This said "the three studies" while the
           ## table listed four -- it was written when Prasse had no theta^perp
           ## row, and nothing tied the sentence to the rows beneath it.
           "moments. Every dataset whose study measured total abundance appears -- ",
           "\\citet{zemb2020absolute}, \\citet{stammler2016adjusting}, ",
           "\\citet{jin2022quantsinglecell}, \\citet{prasse2015soil} and ",
           "\\citet{vandeputte2017quantitative}; only \\citet{mcmurrough2014selex} is absent, having ",
           "no measurement of total abundance at all. \\citet{stammler2016adjusting} reports copies ",
           "per sample rather than per unit mass, so its $\\theta^{\\perp}$ is a per-tube quantity ",
           "while the others are densities. ",
           switch(variant,
             main = paste0(
               "Each study is shown on the single run Table~\\ref{tab:real_studies_result} ",
               "scores -- for \\citet{jin2022quantsinglecell} the pooled cOTU analysis, whose ",
               "measured-scale reference is the only one of its subsets with the power to call ",
               "a feature, for \\citet{prasse2015soil} the QSeq-rescaled counts, and for ",
               "\\citet{stammler2016adjusting} the Pool1-versus-Pool6 comparison. Every entry is ",
               "the hybrid result, that is the arm the selection rule chose for that run; ",
               "Table~\\ref{tab:theta_perp_recovery} gives both arms, and ",
               "Table~\\ref{tab:theta_perp_recovery_stammler} expands ",
               "\\citet{stammler2016adjusting} to all fifteen pool pairs."),
             stammler = paste0(
               "\\citet{stammler2016adjusting} alone, expanded to every pair of its six pools. ",
               "That study is a dilution series, so each pair is a comparison with a known load ",
               "ratio, and these are the same fifteen comparisons the all-pairs TP/FP figure of the ",
               "supplement reports. Every entry is the hybrid result, that is the arm the selection ",
               "rule chose for that pair; Table~\\ref{tab:theta_perp_recovery} gives both arms. ",
               "Pool1 versus Pool6 is the comparison carried in ",
               "Table~\\ref{tab:theta_perp_recovery_hybrid}."),
             arms = paste0(
               "The five headline runs come first, in the order of ",
               "Table~\\ref{tab:theta_perp_recovery_hybrid}, followed by ",
               "\\citet{stammler2016adjusting}'s remaining pool pairs; the comparison column names ",
               "them, and a dash marks a study contributing a single run. Each estimate is followed ",
               "by its 95\\% credible interval. ",
               "Table~\\ref{tab:theta_perp_recovery_hybrid} reports the selected arm alone.")),
           "}"),
    paste0("    \\label{tab:theta_perp_recovery",
           switch(variant, main = "_hybrid", stammler = "_stammler", arms = ""), "}"),
    "    \\vspace{2pt}",
    "    \\begin{minipage}{0.98\\textwidth}",
    paste0("        \\footnotesize ",
           if (show_scale_model)
             paste0("\\textsuperscript{\\dag}The scale model the selection rule chose for ",
                    "that run, and therefore the one the hybrid results of ",
                    "Table~\\ref{tab:real_studies_result} report. ")
           else "",
           "Each run contributes a single interval, so ",
           "the last column indicates whether that interval contains the observed value; it is not a ",
           "coverage rate over replicates. For the Laplace arm the posterior of $\\theta^{\\perp}$ is a ",
           "Gaussian mixture with one component per posterior draw, and its moments and quantiles are ",
           "evaluated analytically rather than from a sample of that mixture."),
    "    \\end{minipage}",
    "\\end{table}")
}

################################################################################
## ASSEMBLE
################################################################################

## T1, T2, T3 and the S-series come from the task scripts sourced at the top of
## this file; the T4 series and T5 are built here, since they belong to the Jin
## objects and the theta^perp record this file already holds.
## BUNDLE ORDER = the order the tables are numbered in. Three deliberate moves
## away from the order these were first written in:
##
##   * the nine-variant table follows the results table directly, in the slot the
##     two-arm table used to hold. Those two arms are two of its nine columns, so
##     the table it replaces was a strict subset of it.
##   * the scale-recovery table moves AHEAD of the Jin block. It is a statement
##     about every study, and sitting after five Jin-only tables made it read as
##     part of that section.
##   * the fifteen-pool-pair variant table is gone -- the Task 5 figure shows the
##     same fifteen comparisons at the same nine settings.
##
## Everything after this list is driven by it: the individual .tex files, the
## bundle, and therefore the S-numbers in the manuscript.
MANUSCRIPT_TABLES <- list(
  T1_dataset_summary        = T1,                                                  # S1
  T2_da_results_hybrid      = T2,                                                  # S2
  S1_ssrv_variants_tpfp     = S1,                                                  # S3
  ## Three views of one table, from one function so they cannot disagree:
  ## the five headline runs for the main text, Stammler's fifteen pool pairs as
  ## the companion to the all-pairs TP/FP figure, and every run in both arms.
  T5_theta_perp_recovery_hybrid   = tab_theta_perp_recovery(theta_perp, "main"),
  T5_theta_perp_recovery_stammler = tab_theta_perp_recovery(theta_perp, "stammler"),
  T5_theta_perp_recovery    = tab_theta_perp_recovery(theta_perp, "arms"),         # S4
  T4_1_jin_dimensions       = tab_jin_dimensions(),                                # S5
  T4_2_jin_discovery_counts = tab_discovery_counts(ssrv, jin),                     # S6
  T4_3_jin_cross_cotu       = tab_cross_representation_full(ssrv, jin, scope = "cotu"),    # S7
  T4_4_jin_cross_features   = tab_cross_representation_full(ssrv, jin, scope = "feature"), # S8
  T4_5_jin_feature_taxonomy = tab_jin_sig_feature_taxonomy(ssrv),                  # S9
  S3_ssrv_variants_jin      = S3                                                   # S10
)

## The taxonomy table covers EVERY called ASV and OTU, unlike the feature-driven
## cross-representation table, which now shows only those whose cOTU is not
## already headed elsewhere. Checked against the full set for that reason.
check_feature_table(MANUSCRIPT_TABLES$T4_5_jin_feature_taxonomy, ssrv)

header <- c(
  "% ---------------------------------------------------------------------------",
  "% Manuscript tables -- GENERATED FILE, do not edit by hand.",
  paste0("% Produced by manuscript_tables.R on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "."),
  paste0("% Jin objects: jin_results_all_levels_", DATESTAMP, ".RData (Jin result reported: ", ARM,
         if (identical(ARM, "hybrid")) ", the selection rule's choice per subset" else "", ")."),
  paste0("% Cross-dataset results: real_data_results.RData, written ",
         format(file.mtime(F_RESULTS), "%Y-%m-%d %H:%M"), "."),
  "% All log fold changes, on both sides, are log2.",
  "%",
  "% Requires in the preamble:",
  "%   \\usepackage{booktabs, longtable, pdflscape, graphicx, multirow, array, rotating}",
  "%   \\usepackage[authoryear]{natbib}   % for \\citet and \\cite",
  "% ---------------------------------------------------------------------------",
  ""
)

for (nm in names(MANUSCRIPT_TABLES)) {
  ## arm_suffix() comes from the sourced Jin file: empty under the hybrid rule,
  ## "_bootstrap" / "_laplace" for a single-arm run, so the two cannot collide.
  writeLines(MANUSCRIPT_TABLES[[nm]], file.path(DIR_OUT, paste0(nm, arm_suffix(), ".tex")))
}
out_all <- file.path(DIR_OUT, paste0("manuscript_tables_", DATESTAMP, arm_suffix(), ".tex"))
writeLines(c(header,
             unlist(lapply(names(MANUSCRIPT_TABLES), function(nm)
               c(paste0("% ---- ", nm, " ----"), MANUSCRIPT_TABLES[[nm]], "")),
               use.names = FALSE)),
           out_all)

message("Wrote ", out_all)
message("Wrote ", length(MANUSCRIPT_TABLES), " individual table files to ", DIR_OUT)

## Echo when run non-interactively, so `Rscript manuscript_tables.R` prints
## paste-ready LaTeX for the whole manuscript.
if (!interactive()) cat(readLines(out_all), sep = "\n")
