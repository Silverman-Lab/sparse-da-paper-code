################################################################################
## TASK 6 -- Figure: the LFC-distribution diagnostic.
##
##   Rscript "Task 6 LFC distribution diagnostic.R"
##
## One panel per analysed subset, showing where the mode of the compositional
## log2 fold changes sits, where the method's fitted scale puts it, and where the
## external measurement says it should be. The gap between those is the quantity
## the sparse scale model has to get right: the method places theta^perp by the
## mode of the full LFC distribution, so when the changed features drag that mode
## away from the truth, the scale is displaced by exactly that much.
##
##   red solid       OBSERVED theta^perp -- the external measurement
##   black longdash  MODE of the LFC distribution, over all features
##   blue twodash    ESTIMATED theta^perp -- posterior mean, hybrid arm
##   pink dotdash    CLR zero
##   grey dashed     zero -- where the mode would sit if the load did not move
##
## The two theta^perp lines are drawn NEGATED, on the compositional LFC axis; see
## the note above observed_perp(). The mode and the CLR line are already on that
## axis. So all four are directly comparable, and the figure reads as: does the
## method's scale (blue) land where the measurement (red) says it should, and how
## far is either from where the naive mode (black) or CLR (pink) would put it?
##
## Coverage: all six datasets, one panel each. Stammler contributes only the pair
## the results tables report, not all fifteen; Jin contributes its POOLED cOTU
## run, the one the manuscript reports. Its six stratified subsets sit behind
## JIN_SUBSET_PANELS, off by default -- see the note above that flag.
##
## This supersedes `real data LFC diagnosis 03212026.R`, which computed the same
## thing but repeated the plotting code once per dataset and wrote to hardcoded
## /storage paths. The per-study data preparation is lifted from it (and, for
## Zemb, Vandeputte and Jin, from theta_perp_check.R, which already had it in a
## ROOT-relative form).
##
## Every subset's feature count is checked against the scale-model record of the
## run the tables report. A panel built on a different feature set would show a
## different mode from the one the method used, which is the one thing this
## figure must not do -- so a mismatch warns loudly rather than passing.
##
## OUTPUTS, in manuscript_figures/
##   lfc_diagnose_<subset>.pdf   one panel per analysed subset. NO title: the
##                               study name belongs in the manuscript caption,
##                               and only the subtitle -- material, contrast --
##                               is drawn.
##   lfc_diagnose_all.pdf        the combined faceted sheet, where the facet
##                               strips DO carry the study name, since that is
##                               the only thing telling the panels apart.
##
## Nothing is refitted; this reads counts only. It is the one task that goes back
## to the raw count tables, because the mode it draws has to be the mode the
## method saw -- which is why every subset's feature count is checked against the
## scale-model record before the panel is written.
##
## This file was manuscript_figures.R until the reporting was split into the six
## numbered tasks; it is the same script under its task name, not a copy, so
## there is no second version of it to fall out of step.
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
if (is.na(ROOT)) stop("no project root found")

## ROOT, DIR_FIG, DATESTAMP, truth_for() and the scale-model record, shared with
## every other task so this figure and the tables describe one run.
source(file.path(ROOT, "task_common.R"))
## sparse_perp(), for the estimated theta^perp line. task_common.R pulls in
## scale_model_record.R, which CALLS sparse_perp() but does not source it, so it
## is not in scope from that chain alone. Reads attributes only -- sourcing this
## fits nothing.
source(file.path(ROOT, "sparse_SSRV.R"))

suppressMessages({
  library(ggplot2); library(dplyr); library(tibble)
  library(phyloseq); library(data.table); library(readxl)
})

################################################################################
## SHARED
################################################################################

## The compositional log2 fold change per feature: mean log2 proportion in the
## treatment group minus the same in the control. The +1 is the pseudocount the
## diagnosis script has always used; it only matters for features that are zero
## in a whole group, which are exactly the ones with no finite LFC otherwise.
## The pseudocount, defined ONCE. lfc_vec() and clr_zero() must use the same one
## or the CLR line stops living on the axis the histogram is drawn against, and
## nothing in the figure would say so.
##
## 0.5, matching the Dirichlet concentration the method uses (Y + 0.5 in the
## resampling that feeds sparse_auto). It was +1, which put the histogram on a
## different LFC construction from the theta^perp drawn over it: on selex, where
## 74.3% of the selected arm is zeros and the unselected arm has none, that
## disagreement was 1.04 of the 2.6 log2 units between the mode and the fitted
## theta^perp.
##
## The remaining gap is NOT a pseudocount effect and this does not close it: the
## method takes the mode PER DIRICHLET DRAW and averages, so it carries a Jensen
## term -- E[log p] < log E[p], large for a near-zero proportion -- that no
## plug-in construction reproduces. On selex that is a further 1.55.
LFC_PSEUDOCOUNT <- 0.5

lfc_vec <- function(Y, X) {
  Y <- as.matrix(Y)
  prop <- sweep(Y + LFC_PSEUDOCOUNT, 2, colSums(Y), "/")
  rowMeans(log2(prop[, X == 1, drop = FALSE])) -
    rowMeans(log2(prop[, X == 0, drop = FALSE]))
}

## Where a CLR normalisation puts zero, ON THE SAME AXIS as the plotted LFC.
##
## CLR subtracts each sample's mean log2 proportion, so
##   clr_lfc_d = lfc_d - (mean_T gm - mean_C gm),   gm_n = mean_d log2(prop_dn).
## A feature CLR calls unchanged therefore has lfc_d = mean_T gm - mean_C gm, and
## that is the line. Equivalently CLR asserts theta^perp = -(that shift), which
## is what makes it comparable with the measured value.
##
## Same pseudocount and same proportions as lfc_vec() -- LFC_PSEUDOCOUNT, shared
## rather than repeated, or the line would not live on the axis it is drawn
## against.
clr_zero <- function(Y, X) {
  Y <- as.matrix(Y)
  prop <- sweep(Y + LFC_PSEUDOCOUNT, 2, colSums(Y), "/")
  gm <- colMeans(log2(prop))
  mean(gm[X == 1]) - mean(gm[X == 0])
}

## Mode of a Sheather-Jones kernel density, which is what the scale model uses.
dens_mode <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  d <- density(x, bw = "SJ")
  d$x[which.max(d$y)]
}

## The three quantities the panel draws, plus the two theta^perp values and the
## CLR shift, which are still computed for the run log but no longer drawn.
##
## WHY theta^perp IS NO LONGER ON THIS FIGURE. Both theta^perp lines were removed
## because neither is commensurable with the histogram they were drawn over. The
## fitted theta^perp is minus the mode of the DIRICHLET-RESAMPLED LFC, and for a
## count of y the resampled log2 proportion sits (digamma(y+0.5) - log(y+0.5))/ln2
## bits below the plug-in one -- -1.83 bits at y = 0, -0.03 by y = 24. Where the
## zeros are split unevenly between arms that offset does not cancel: on selex,
## 74.3% zeros in the selected arm against 0% in the other, it displaced the line
## by 1.56 log2 units against a histogram that has no such term. The observed
## theta^perp went with it because it was the reference the fitted one was being
## judged against, and judging it against a different construction is the whole
## problem.
##
## What replaces it is a comparison that is internally consistent: every line
## below is a mode of THIS panel's LFC vector, or zero. Nothing is converted.
##
##   null   the mode over null features only -- the truth set removed. This is
##          the reference: where the bulk sits once the features that genuinely
##          moved are taken out.
##   all    the mode over every feature, which is the quantity the scale model
##          estimates. Its distance from `null` is how far the changed features
##          drag the estimate.
##   zero   where the mode would sit if the total load had not moved at all.
lfc_modes <- function(lfc, features, truth, clr = NA_real_, obs = NA_real_,
                      est = NA_real_) {
  full <- dens_mode(lfc)
  null <- if (!length(truth) || !any(features %in% truth)) NA_real_
          else dens_mode(lfc[!(features %in% truth)])
  neg <- function(z) if (!is.finite(z)) NA_real_ else -z
  c(null = null, all = full, zero = 0,
    ## kept for the run log only -- absent from LINE_KEYS, so never drawn
    est = neg(est), clr = clr, obs = neg(obs))
}

## Three drawn lines, and colour is never the only thing separating them: the
## advisor reading these has green colour blindness, so the palette avoids green
## and every line carries its own linetype as well. All three linetypes are
## distinct, so the figure survives being read in greyscale or by a deuteranope.
## The null mode is the thickest and the only solid one, because it is the
## reference the all-feature mode is being judged against.
##
## `est`, `clr` and `obs` deliberately have no entry here: an unlisted key is not
## drawn, which is what keeps them out of the panels and out of the legend while
## lfc_modes() still reports them to the log.
LINE_COL <- c(null = "red",   all = "blue",   zero = "black")
LINE_LTY <- c(null = "solid", all = "dashed", zero = "dotted")
LINE_LWD <- c(null = 1.1,     all = 0.9,      zero = 0.6)
## Plain strings, because these are the legend KEYS as well: setNames() and
## `breaks` have to match them exactly. They also go into the run log, which is
## why they stay readable as ASCII.
LINE_LAB <- c(null = "mode, null features only",
              all  = "mode, all features",
              zero = "zero")

## The order lines are drawn and listed in, everywhere. One definition, so a
## panel and the combined legend cannot fall out of step.
LINE_KEYS <- c("null", "all", "zero")

## `band` is accepted and ignored: it was the measurement interval drawn behind
## the observed theta^perp, and that line is gone. Kept in the signature so the
## build loop and the `built` list still carry it -- reinstating the line means
## restoring the annotate() below, not re-plumbing the caller.
##
## `title` is accepted and NOT DRAWN, for the same reason: the single-panel PDFs
## carry no study name, because each is placed in the manuscript under its own
## caption and the study/year printed inside the plot box only repeats it. The
## caller still passes it and `built` still stores it -- it is what labels the
## facets of the combined sheet, where the panels DO have to be told apart.
lfc_panel <- function(d, m, title, sub, band = c(NA_real_, NA_real_)) {
  p <- ggplot(d, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), colour = "grey30",
                   fill = "white", bins = 60, linewidth = 0.2) +
    geom_density(alpha = 0.2, bw = 0.5, fill = "#FF6666", colour = "grey40") +
    labs(subtitle = sub,
         x = expression("LFC (" * hat(theta)^"||" * ")"), y = "Density") +
    theme_bw(base_size = 10) +
    theme(plot.subtitle = element_text(size = 8, colour = "grey30"))
  ## Zero is drawn from LINE_KEYS like the other two rather than as its own
  ## geom_vline, so it appears in the legend and cannot be styled out of step
  ## with the panels.
  for (k in LINE_KEYS) {
    if (is.finite(m[[k]]))
      p <- p + geom_vline(xintercept = m[[k]], colour = LINE_COL[[k]],
                          linetype = LINE_LTY[[k]], linewidth = LINE_LWD[[k]])
  }
  p
}

################################################################################
## TRUTH SETS
##
## The five non-Jin studies use exactly the sets `real data results.R` scores
## against, read from the record it writes. Jin's POOLED cOTU run is in that
## record too, and is the run the manuscript reports; the STRATIFIED subsets are
## not, and their measured-scale gold reference is empty at both locations, so
## their null set is the project's declared truth -- COTU-CM-2002, the one cOTU
## the pooled reference calls -- carried down to the ASV and OTU representations
## through the mapping. The record wins wherever it has an entry, so the pooled
## panel and the results tables score against the same set. Stated here rather
## than buried: with one feature removed from ~190 this barely moves the null
## mode, and the Jin panels should be read as showing the full mode against
## zero.
################################################################################
## truth_sets and truth_for() come from task_common.R, which reads the sets
## `real data results.R` scores against.
JIN_DECLARED_TRUTH <- "COTU-CM-2002"

## Jin: the declared cOTU truth, mapped down to whichever representation the
## subset uses.
jin_truth_for <- function(level, loc) {
  if (identical(level, "cOTU")) return(JIN_DECLARED_TRUTH)
  f <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates",
                 paste0("jin_results_all_levels_", DATESTAMP, ".RData"))
  if (!file.exists(f)) return(character(0))
  tl <- readRDS(f)[[level]][[loc]]$tax_lookup
  if (is.null(tl$cOTU)) return(character(0))
  as.character(tl$category[tl$cOTU %in% JIN_DECLARED_TRUTH])
}

################################################################################
## THE MEASURED theta^perp, PLACED ON THE LFC AXIS
##
## Read from theta_perp_check.csv, the same record the manuscript's scale-recovery
## table reports, so the figure and that table cannot disagree.
##
## The axis here is the COMPOSITIONAL LFC, and a feature unchanged in absolute
## terms satisfies theta = theta^|| + theta^perp = 0, so theta^|| = -theta^perp.
## A measured load shift of +0.86 therefore predicts the null features pile up at
## -0.86, not at +0.86. The line is drawn at -obs for that reason -- the same
## convention the CLR line uses, where CLR asserts theta^perp = -(its zero).
##
## Only the three studies with an external measurement have one; Jin has it at
## cOTU level only, so its ASV and OTU panels get no line.
################################################################################
F_PERP <- file.path(ROOT, "theta_perp_check", "theta_perp_check.csv")
perp_rec <- if (file.exists(F_PERP)) read.csv(F_PERP, stringsAsFactors = FALSE) else NULL

## Jin also reports the replicate spread of each scale measurement, so its
## measured theta^perp gets a band as well as a line. Same helper the
## scale-recovery table uses, so the band and the interval printed there are the
## same numbers. Assay noise only -- not between-mouse variation.
source(file.path(ROOT, "jin_scale_uncertainty.R"))
perp_ci <- jin_scale_uncertainty(file.path(ROOT, "2022_Jin_natureComm_technicalReplicates"))

## On the LFC axis, and therefore negated and with the ends swapped.
observed_perp_band <- function(key, pair) {
  if (is.null(perp_ci)) return(c(NA_real_, NA_real_))
  r <- perp_ci[perp_ci$dataset == key & !is.na(pair) & perp_ci$pair == pair, ]
  if (!nrow(r)) return(c(NA_real_, NA_real_))
  c(-r$upper[1], -r$lower[1])
}

observed_perp <- function(key, pair) {
  if (is.null(perp_rec)) return(NA_real_)
  r <- perp_rec[perp_rec$dataset == key &
                ((is.na(pair) & is.na(perp_rec$pair)) |
                 (!is.na(pair) & !is.na(perp_rec$pair) & perp_rec$pair == pair)), ]
  if (!nrow(r) || all(is.na(r$obs))) return(NA_real_)
  r$obs[!is.na(r$obs)][1]
}

################################################################################
## THE ANALYSED FIT -- the single source for both what the method estimated and
## which features it estimated it on
##
## Everything below is read off the study's saved sparse_sens_*.RData. sparse_perp()
## and attr(., "arm_out") read attributes only and fit nothing, so this is cheap.
##
## sens[[1]] is the first bandwidth variant, which is the fit the results tables
## report; the other elements are the bandwidth sensitivity arms.
##
## WHY THE FEATURE SET COMES FROM HERE. Each prep below used to re-implement its
## study's filter, and each one drifted: Zemb was two filter generations behind
## (411 features against an analysed 256, and it never dropped the extreme-dose
## tube), Stammler was still on `mean > 5` against a run using 0.05% relative
## abundance (27 against 25), and Prasse was filtering on raw prevalence for a
## panel pointed at a run no table reports (175 against 85). Three studies, three
## different rules, none of them the rule its DA script now uses.
##
## Re-implementing three filters a fourth time only sets up the fourth drift. The
## fit already records exactly which features it was computed on, as the rownames
## of arm_out, so the preps below now build the UNFILTERED table and the build
## loop subsets it to those names. D then matches by construction rather than by
## the two files happening to agree.
##
## The preps must therefore return a SUPERSET of the analysed features, in the
## same id space -- which is what agglomerating (tax_glom) but not filtering
## gives. apply_analysed_features() checks that and stops if it is not true,
## because a prep that has lost features is exactly the case a subset would
## paper over.
################################################################################
PERP_FIT <- list(
  McMurrough2014 = list(dir = DIR_SELEX),
  Zemb2020       = list(dir = DIR_ZEMB),
  Stammler2016   = list(dir = DIR_STAMMLER),
  Vandeputte2017 = list(dir = DIR_VANDEPUTTE),
  Prasse2015     = list(dir = DIR_PRASSE),
  Jin2022        = list(dir = DIR_JIN_))

## Where each study's fit files put the subset tag. Jin appends it to the stem
## (sparse_sens_ASV_prox_<stamp>), Prasse appends it after the stamp
## (sparse_sens_<stamp>_mean5_QSeq); everything else has no subset.
##
## Jin's POOLED cOTU run is the exception within Jin: sparse_DA_09012026.R writes
## it under the BARE name (sparse_sens_<stamp>.RData), keeping the location out
## of the filename, and only the stratified runs carry a tag. So pair = "both"
## resolves to the same pattern a study with no subset uses -- it is a `pair` for
## the records' sake (that is the key in scale_model_record.rds and in the
## theta^perp record), not a filename tag.
fit_pattern <- function(key, pair) {
  if (is.na(pair)) return("sparse_sens_%s.RData")
  switch(key,
    Jin2022    = if (identical(pair, "both")) "sparse_sens_%s.RData"
                 else paste0("sparse_sens_", pair, "_%s.RData"),
    Prasse2015 = paste0("sparse_sens_%s_", pair, ".RData"),
    "sparse_sens_%s.RData")
}

.fit_cache <- new.env(parent = emptyenv())

## The fit itself, cached: three callers want it and it is the largest read here.
analysed_fit <- function(key, pair) {
  spec <- PERP_FIT[[key]]
  if (is.null(spec)) return(NULL)
  ck <- paste0(key, "|", if (is.na(pair)) "NA" else pair)
  if (!is.null(.fit_cache[[ck]])) return(.fit_cache[[ck]]$v)
  f <- tryCatch(stamped(spec$dir, fit_pattern(key, pair),
                        what = paste0("fit_", ck), required = FALSE),
                error = function(e) NA_character_)
  v <- if (!is.na(f) && file.exists(f))
         tryCatch(readRDS(f)[[1]], error = function(e) NULL) else NULL
  assign(ck, list(v = v), envir = .fit_cache)
  v
}

## The feature ids the run was computed on, in the order the fit holds them.
analysed_features <- function(key, pair) {
  o <- analysed_fit(key, pair)
  if (is.null(o)) return(NULL)
  ao <- attr(o, "arm_out")
  if (is.null(ao) || is.null(ao$laplace)) return(NULL)
  rownames(ao$laplace)
}

## Subset a prep's table to the analysed features. Stops rather than warns: a
## panel built on the wrong features is the one failure this whole file exists
## to detect, and it is not something to emit a figure about.
apply_analysed_features <- function(Y, key, pair, tag) {
  want <- analysed_features(key, pair)
  if (is.null(want)) {
    warning(tag, ": no saved fit found, so the panel keeps its own ", nrow(Y),
            " features and may not match the analysed run", call. = FALSE,
            immediate. = TRUE)
    return(Y)
  }
  miss <- setdiff(want, rownames(Y))
  if (length(miss))
    stop(tag, ": the analysed run used ", length(want), " features but this ",
         "panel's preparation produced only ", nrow(Y), ", missing ",
         length(miss), " of them (", paste(utils::head(miss, 4), collapse = ", "),
         if (length(miss) > 4) ", ..." else "", "). The preparation must build a ",
         "SUPERSET of the analysed features -- agglomerate, but do not filter.",
         call. = FALSE)
  Y[want, , drop = FALSE]
}

## The estimated theta^perp: the HYBRID row, the arm the selection rule actually
## delivered, so it is the scale the reported DA results were computed against.
## It always duplicates whichever arm was chosen, so this line is never a fourth
## independent estimate.
estimated_perp <- function(key, pair) {
  o <- analysed_fit(key, pair)
  if (is.null(o)) return(NA_real_)
  tab <- tryCatch(suppressWarnings(sparse_perp(o)), error = function(e) NULL)
  if (is.null(tab) || !"arm" %in% names(tab) || !any(tab$arm == "hybrid"))
    return(NA_real_)
  as.numeric(tab$mean[tab$arm == "hybrid"][1])
}

################################################################################
## EXPECTED DIMENSIONS -- the run each panel must match
##
## scale_records comes from task_common.R. expected_D() returns NA rather than
## stopping where dims_for() would: a subset with no record still gets its panel,
## unchecked, because a missing record is a gap in the bookkeeping and not a
## reason to withhold the figure.
################################################################################
## Both dimensions, from the record for this run. N matters as much as D and was
## not being checked: the feature subsetting guarantees D, so a panel can match
## the analysed feature set exactly while computing every LFC over the wrong set
## of SAMPLES. That is not hypothetical -- Zemb's prep carried a sample drop the
## DA script had removed, and nothing in this file would have said so.
expected_dims <- function(key, pair) {
  if (is.null(scale_records)) return(c(D = NA_integer_, N = NA_integer_))
  r <- scale_records[scale_records$dataset == key &
                     ((is.na(pair) & is.na(scale_records$pair)) |
                      (!is.na(pair) & !is.na(scale_records$pair) &
                       scale_records$pair == pair)), ]
  if (!nrow(r)) c(D = NA_integer_, N = NA_integer_)
  else c(D = r$D[1], N = r$N[1])
}

expected_D <- function(key, pair) unname(expected_dims(key, pair)[["D"]])

################################################################################
## PER-STUDY PREPARATION
##
## Each returns list(Y = counts D x N, X = 0/1 group vector). The filtering in
## each is the filtering of the run the tables report -- that is what the D check
## below verifies.
################################################################################

prep_mcmurrough <- function() {
  data("selex", package = "ALDEx2", envir = environment())
  Y <- get("selex", envir = environment())
  list(Y = Y, X = as.numeric(c(rep("NS", 7), rep("S", 7)) == "S"),
       g1 = "selected", g0 = "not selected")
}

prep_zemb <- function() {
  source(file.path(ROOT, "zemb_parse.R"), local = TRUE)
  z <- parse_2020_zemb_microbiologyOpen_spike()

  ## SAMPLE SET -- empty since 2026-09-01, matching sparse_DA_09012026.R:246,
  ## "THIS script drops nothing further: the 9 samples the parser returns are the
  ## 9 that every method is fitted on." The 100 ul extreme-dose tube (oz1802-13)
  ## was excluded until then; the DA script keeps the reasoning for that, and for
  ## its reversal, at its SAMPLE EXCLUSION block. zemb_parse.R still drops the
  ## no-feces control (oz1802-19) at the source, which is where 10 becomes 9.
  ##
  ## Kept as an explicit vector rather than deleted, because this is a SAMPLE-set
  ## setting and the build loop's feature subsetting canNOT repair it: N is not
  ## something the analysed feature names can fix, so getting it wrong means
  ## every LFC in the histogram is computed over a different set of tubes. The
  ## dimension check in the build loop now tests N as well as D for that reason.
  ##
  ## Counts, metadata and qPCR are subset TOGETHER -- the scale below is
  ## positional, so dropping a sample from one table only is exactly how the
  ## scale and the counts come apart without saying so.
  ZEMB_DROP_SAMPLES <- character(0)
  stopifnot(all(ZEMB_DROP_SAMPLES %in% colnames(z$counts)))
  keep_s <- !(colnames(z$counts) %in% ZEMB_DROP_SAMPLES)
  z$counts   <- z$counts[ , keep_s, drop = FALSE]
  z$metadata <- z$metadata[keep_s, , drop = FALSE]
  z$scale    <- z$scale[   keep_s, , drop = FALSE]

  md <- cbind(z$metadata, spiked = z$metadata$added_Coli > 0)
  ## UNFILTERED on purpose -- the build loop subsets to the analysed feature set.
  ## Seq135, the spiked-in internal standard, is not a community taxon and the DA
  ## script drops it; it simply is not among the analysed names, so it falls out
  ## of the subset without needing its own rule here.
  ph <- phyloseq(otu_table(t(z$counts), taxa_are_rows = FALSE), sample_data(md))
  ## Density, not copies per tube: these two arms were weighed out 3-5x apart,
  ## so the per-tube total measures the pipetting. Same choice zemb_parse.R makes.
  dens <- z$scale$total16S_per_mg
  names(dens) <- rownames(z$metadata)
  list(Y = t(otu_table(ph)), X = as.numeric(sample_data(ph)$spiked == "TRUE"),
       g1 = "spiked", g0 = "unspiked",
       scale = dens[sample_names(ph)], scale_basis = "qPCR 16S copies/mg (density)")
}

prep_stammler <- function() {
  source(file.path(ROOT, "2016_stammler_microbiome_micehuman", "parse.R"), local = TRUE)
  out <- parse_2016_stammler_microbiome_micehuman()
  i1 <- which(out$metadata$Treatment == "Pool1")
  i6 <- which(out$metadata$Treatment == "Pool6")
  ph <- phyloseq(otu_table(t(out$counts[, c(i1, i6)]), taxa_are_rows = FALSE),
                 sample_data(out$metadata[c(i1, i6), ]),
                 tax_table(as.matrix(out$tax)))
  ## Agglomerate but do NOT filter: tax_glom defines the id space the analysed
  ## feature names live in, and the build loop does the subsetting. Leaving the
  ## filter out is what keeps this panel correct across filter changes -- a
  ## `mean(v) > 5` written here once disagreed with the DA script's then-current
  ## 0.05%-relative-abundance rule, 27 features against its 25, in silence.
  red <- tax_glom(ph, taxrank = "Genus")
  ## Copies per SAMPLE -- the only form this study reports; there is no mass
  ## column to turn it into a density, so theta^perp here is a per-tube quantity.
  sc <- as.numeric(out$scale$`16S rDNA copies per sample`)
  names(sc) <- as.character(out$scale$SampleID)
  list(Y = t(otu_table(red)),
       X = as.numeric(sample_data(red)$Treatment == "Pool6"),
       g1 = "Pool6", g0 = "Pool1",
       scale = sc[as.character(sample_data(red)$SampleID)],
       scale_basis = "16S rDNA copies per sample (per-tube)")
}

prep_vandeputte <- function() {
  md <- read.csv(file.path(ROOT, "Vandeputte", "data",
                           "41586_2017_BFnature24460_MOESM10_ESM.csv"))
  OTU  <- fread(file.path(ROOT, "Vandeputte", "data", "OTU_nochim.csv"))[-c(1:40), ]
  taxa <- fread(file.path(ROOT, "Vandeputte", "data", "taxa_assignments.csv"),
                header = TRUE) %>% column_to_rownames("V1")
  tid <- paste0("sp", seq_len(nrow(taxa))); rownames(taxa) <- tid
  OTU <- as.data.frame(OTU) %>% column_to_rownames("V1")
  colnames(OTU) <- tid; rownames(OTU) <- paste0("sa", seq_len(nrow(OTU)))
  ph <- phyloseq(otu_table(OTU, taxa_are_rows = FALSE), sample_data(md),
                 tax_table(as.matrix(taxa)))
  red <- filter_taxa(tax_glom(ph, taxrank = "Genus"),
                     function(v) sum(v > 1) > .2 * length(v), prune = TRUE)
  list(Y = t(otu_table(red)),
       X = as.numeric(sample_data(red)$Health.status == "CD"),
       g1 = "Crohn's", g0 = "healthy",
       scale = as.numeric(sample_data(red)$Average.cell.count..per.gram.of.frozen.feces.),
       scale_basis = "flow cytometry, cells per gram (density)")
}

## Jin: one prep for all six subsets, the same construction theta_perp_check.R
## uses, so the feature sets match the fits the tables report.
jin_env <- local({
  loaded <- NULL
  function() {
    if (!is.null(loaded)) return(loaded)
    lp <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates")
    dat <- read_xlsx(file.path(lp, "Supplementary Data 8.xlsx"),
                     sheet = "Sequencing-determined counts", skip = 1)
    cnt <- as.matrix(dat[, 2:56]); rownames(cnt) <- dat$`cOTU-ID`
    scl <- read_xlsx(file.path(lp, "Supplementary Data 8.xlsx"),
                     sheet = "Absolute total abundance", skip = 1)
    meta <- data.frame(row.names = scl$Sample,
                       diet = substr(scl$Sample, 1, 2),
                       location = substr(scl$Sample, 4, 7),
                       sample_type = vapply(strsplit(scl$Sample, "-"),
                                            function(z) if (length(z) > 1) z[2] else "cell",
                                            character(1)),
                       stringsAsFactors = FALSE)
    cnt  <- cnt[, colnames(cnt) %in% rownames(meta)]
    meta <- meta[rownames(meta) %in% colnames(cnt), ]
    as_mat <- function(d, idcol) {
      m <- as.matrix(d[, setdiff(colnames(d), idcol), drop = FALSE]); mode(m) <- "numeric"
      m[is.na(m)] <- 0; rownames(m) <- d[[idcol]]; m
    }
    dens <- scl$`Absolute total abudance (copies/mg)`
    names(dens) <- scl$Sample
    loaded <<- list(
      dens = dens,
      cOTU = t(cnt),
      ASV  = t(as_mat(read_xlsx(file.path(lp, "Supplementary Data 5.xlsx"),
                                sheet = "VA group", skip = 1), "ASV-ID")),
      OTU  = t(as_mat(read_xlsx(file.path(lp, "Supplementary Data 7.xlsx"),
                                sheet = "VA group", skip = 1), "OTU-ID")),
      meta = meta)
    loaded
  }
})

## `loc` is "prox", "dist", or "both" -- the POOLED run, all sixteen cell samples
## with the location clause dropped. Same three subsets sparse_DA_09012026.R
## builds (its `subsets` vector), and the same mask: the pooled run is the same
## eight mice sampled twice, not a third location.
prep_jin <- function(level, loc) {
  e <- jin_env()
  COUNTS <- e[[level]]
  md <- e$meta[match(rownames(COUNTS), rownames(e$meta)), , drop = FALSE]
  keep <- rownames(COUNTS)[!is.na(md$diet) & md$diet %in% c("VS", "VD") &
                           md$sample_type == "cell" &
                           (if (identical(loc, "both")) TRUE
                            else md$location == loc)]
  ph <- phyloseq(otu_table(COUNTS[keep, , drop = FALSE], taxa_are_rows = FALSE),
                 sample_data(e$meta[keep, , drop = FALSE]))
  red <- filter_taxa(ph, function(v) sum(v >= 5) >= 3, prune = TRUE)
  ## cOTU only: BarBIQ quantifies absolute abundance for cOTUs, and the ASV and
  ## OTU tables are gene-copy measurements with no measured scale.
  sc <- if (identical(level, "cOTU")) e$dens[sample_names(red)] else NULL
  list(Y = t(otu_table(red)), X = as.numeric(sample_data(red)$diet == "VD"),
       g1 = "VA-deficient", g0 = "VA-sufficient",
       scale = sc,
       scale_basis = if (is.null(sc)) NA_character_ else "BarBIQ copies/mg (density)")
}

## Prasse, under the rowMeans(raw counts) >= 5 filter -- the run the results
## tables report. Two subtleties carried over from sparse_DA_..._mean5.R: the
## counts are the QSeq-transformed table (qPCR-rescaled), but the filter is
## applied to the RAW genus table, because after QSeq every genus mean is between
## 1.6e3 and 7.3e7 and a mean-5 rule there would keep all 561. The two glommed
## tables are joined on genus NAME, since tax_glom picks its own archetype taxon
## on each call.
##
## No threshold constant is needed here: this prep deliberately does not filter,
## and the build loop subsets to the analysed feature set read off the saved fit.
## (A PREV_CUT <- 0.10 sat here unused after the filter moved out of this file.)
prep_prasse <- function() {
  suppressMessages(library(QSeq))
  dp  <- file.path(ROOT, "2015_prasse_appliedEnvMicrobiology_wetlands")
  TFW <- readRDS(file.path(dp, "Data", "TFW_comdat.RDS"))
  ## fileEncoding is REQUIRED: the file begins with a UTF-8 BOM, so a plain
  ## read.csv() names the first column "X...ID" and meta$ID comes back NULL. See
  ## the sample-alignment note below for what that silently did to this panel.
  meta <- read.csv(file.path(dp, "Data", "TFW_metadata.csv"),
                   fileEncoding = "UTF-8-BOM")
  rownames(meta) <- meta$ID

  ## All 57 samples, matched by name. This block used to intersect down to 54 and
  ## report that the two files disagreed about three samples; there was no such
  ## disagreement. TFW_metadata.csv begins with a UTF-8 BOM, so the plain
  ## read.csv() above named its first column "X...ID", meta$ID came back NULL,
  ## `rownames(meta) <- NULL` silently reset the rownames to the POSITIONS 1..57,
  ## and those matched TFW's names (1-21, 25-60) for only 54 of them. The three
  ## "missing" samples were an artefact of reading the file wrong. The DA script
  ## documents the same trap at sparse_DA_09012026_rawVSqseq_mean5.R and
  ## reads with fileEncoding for exactly this reason -- see the read.csv above.
  stopifnot(!is.null(meta$ID))
  stopifnot(setequal(sample_names(TFW), rownames(meta)))
  sample_data(TFW) <- sample_data(meta[sample_names(TFW), , drop = FALSE])

  Q <- QSeq(TFW, "QPCR_16S")
  ## Orient explicitly rather than assuming: the original code transposed
  ## unconditionally, which silently puts taxon names where sample names belong
  ## whenever QSeq returns a samples-as-rows table.
  as_samples_by_taxa <- function(ps) {
    m <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) t(m) else m
  }
  ## sp ids assigned once off the SHARED otu-level taxonomy, as the DA script
  ## does: QSeq rescales counts without adding, dropping or renaming taxa.
  stopifnot(identical(taxa_names(TFW), taxa_names(Q)))
  taxa <- as.matrix(tax_table(Q))
  tid  <- paste0("sp", seq_len(nrow(taxa)))
  rownames(taxa) <- tid
  build <- function(ps) {
    ot <- as_samples_by_taxa(ps); colnames(ot) <- tid
    phyloseq(otu_table(ot, taxa_are_rows = FALSE), sample_data(Q), tax_table(taxa))
  }
  genus_qseq <- tax_glom(build(Q),   taxrank = "Genus")
  genus_raw  <- tax_glom(build(TFW), taxrank = "Genus")

  ## RE-KEY QSeq ONTO THE RAW ARCHETYPE IDS, as the DA script does at
  ## sparse_DA_09012026_rawVSqseq_mean5.R. tax_glom names each genus after
  ## its most abundant constituent OTU, and because QSeq rescales samples by
  ## different factors it picks a different archetype for 543 of the 561 genera.
  ## The analysed feature names are therefore RAW archetype ids. Without this the
  ## panel's ids are the QSeq archetypes and 106 of the analysed 181 do not exist
  ## in it -- which is precisely what apply_analysed_features() refused to plot.
  G_qseq <- as.character(tax_table(genus_qseq)[, "Genus"])
  G_raw  <- as.character(tax_table(genus_raw)[,  "Genus"])
  stopifnot(setequal(G_qseq, G_raw), !anyDuplicated(G_qseq), !anyDuplicated(G_raw))
  M <- as.matrix(otu_table(genus_qseq))         # samples x taxa
  colnames(M) <- G_qseq
  M <- M[sample_names(genus_raw), G_raw, drop = FALSE]
  colnames(M) <- taxa_names(genus_raw)
  phylo_genus <- phyloseq(otu_table(M, taxa_are_rows = FALSE),
                          sample_data(genus_raw), tax_table(genus_raw))

  ## No filter: the build loop subsets to the analysed feature set. The rule that
  ## used to be here was a prevalence cut on the RAW genus table; the panel now
  ## points at mean5_QSeq -- the run T1 and T2 report -- and takes its features
  ## from that fit, so it cannot drift when the filter rule changes again.
  ## CONTRAST DIRECTION, copied from the DA script rather than chosen here.
  ## sparse_DA_09012026_rawVSqseq_mean5.R builds its contrast as
  ##     X_0 <- as.numeric(X == tail(X, 1))
  ## so group 1 is whichever TYPE the LAST sample happens to carry -- NAT in this
  ## phyloseq ordering. The analysed run, its LFCs and its theta^perp are all on
  ## NAT - REST. This panel used to hard-code TYPE == "REST", which mirrored the
  ## histogram against the fit whose lines are drawn on it, and put the observed
  ## theta^perp line on the wrong side of zero (+0.439 in the record against
  ## -0.439 recomputed here -- same magnitude, opposite sign).
  ##
  ## Note this is the direction the run HAS, not necessarily the one it wants:
  ## the DA script names indsC = NAT the control and indsT = REST the treatment,
  ## which reads as REST - NAT. tail(X, 1) is what actually decides it, and it
  ## disagrees. Flagged rather than flipped -- flipping here alone would put this
  ## figure at odds with every table built from the same fit.
  Xc <- as.character(sample_data(phylo_genus)$TYPE)
  grp1 <- tail(Xc, 1)
  stopifnot(identical(grp1, "NAT"))
  list(Y = t(otu_table(phylo_genus)),
       X = as.numeric(Xc == grp1),
       g1 = "natural", g0 = "restored",
       scale = as.numeric(sample_data(phylo_genus)$QPCR_16S),
       scale_basis = "qPCR 16S copies")
}

################################################################################
## THE PANEL LIST
################################################################################
PANELS <- list(
  list(key = "McMurrough2014", pair = NA_character_, label = "McMurrough et al (2014)",
       sub = "selex", prep = prep_mcmurrough),
  list(key = "Zemb2020", pair = NA_character_, label = "Zemb et al (2020)",
       sub = "pig faeces", prep = prep_zemb),
  ## ASCII: the default pdf() device font has no umlaut and silently drops it.
  ## The LaTeX tables spell it correctly; only the figure label is transliterated.
  list(key = "Stammler2016", pair = "Pool6", label = "Stammler et al (2016)",
       sub = "mock community", prep = prep_stammler),
  ## mean5_QSeq, not mean5. The rawVSqseq QSeq pass is what T1 and T2 report;
  ## the single-pass `mean5` run is the one audit_da_runs.R lists under
  ## "recorded but reported by no table", so pointing here at `mean5` would have
  ## the panel diagnosing a fit that appears nowhere else in the manuscript.
  ## (Under the August prevalence filter the two ran at D = 181 and D = 85; the
  ## mean filter will give different numbers, but the same argument holds.)
  list(key = "Prasse2015", pair = "mean5_QSeq", label = "Prasse et al (2015)",
       sub = "soil", prep = prep_prasse),
  list(key = "Vandeputte2017", pair = NA_character_, label = "Vandeputte et al (2017)",
       sub = "human faeces", prep = prep_vandeputte)
)
## JIN: THE POOLED cOTU RUN IS THE PANEL.
##
## The manuscript reports Jin on its POOLED cOTU run -- both cecal locations
## analysed together -- and that is the row DATASET_META carries (pair = "both")
## and the only Jin run with a truth set to score against. So it is the run this
## figure diagnoses, and its subtitle names no location: pooling means the panel
## is not about one.
##
## The six stratified subsets -- three feature representations x two locations --
## are still buildable, because they are the subject of the Jin subset tables and
## the argument there is precisely that the representations disagree. They are OFF
## by default so that the combined sheet is one panel per DATASET, matching T1 and
## T2. Set this TRUE to get them back; each writes its own
## lfc_diagnose_Jin2022_<subset>.pdf, and the combined sheet grows to twelve.
##
## The six they used to write were deleted from manuscript_figures/ when this
## went FALSE, rather than left to go stale: nothing here removes a panel it has
## stopped building, so a leftover file would sit there looking current.
JIN_SUBSET_PANELS <- FALSE

PANELS[[length(PANELS) + 1]] <- list(
  key = "Jin2022", pair = "both", label = "Jin et al (2022)",
  sub = "cOTU",
  prep = function() prep_jin("cOTU", "both"),
  jin = c(l = "cOTU", o = "both"))

if (JIN_SUBSET_PANELS) {
  for (lv in c("cOTU", "ASV", "OTU")) for (lo in c("prox", "dist")) {
    PANELS[[length(PANELS) + 1]] <- list(
      key = "Jin2022", pair = if (lv == "cOTU") lo else paste0(lv, "_", lo),
      label = "Jin et al (2022)",
      sub = paste0(lv, ", ", if (lo == "prox") "proximal" else "distal", " cecum"),
      prep = local({ l <- lv; o <- lo; function() prep_jin(l, o) }),
      jin = c(l = lv, o = lo))
  }
}

################################################################################
## BUILD
################################################################################
built <- list(); skipped <- character(0)
for (p in PANELS) {
  tag <- paste0(p$key, if (is.na(p$pair)) "" else paste0("_", p$pair))
  if (is.null(p$prep)) { skipped <- c(skipped, paste0(tag, " (no preparation implemented)")); next }
  r <- tryCatch(p$prep(), error = function(e) e)
  if (inherits(r, "error")) {
    skipped <- c(skipped, paste0(tag, " (", conditionMessage(r), ")")); next
  }
  Y <- r$Y; X <- r$X
  ## Subset to the features the run was actually computed on. The preps build an
  ## unfiltered, agglomerated table; this is where the analysed feature set is
  ## imposed, so the mode this panel draws is the mode the method saw.
  Y <- tryCatch(apply_analysed_features(Y, p$key, p$pair, tag),
                error = function(e) e)
  if (inherits(Y, "error")) {
    skipped <- c(skipped, paste0(tag, " (", conditionMessage(Y), ")")); next
  }
  ## Kept as a backstop on the record rather than on the fit: the subset above
  ## guarantees D matches the FIT, and this says whether the fit matches the
  ## scale-model RECORD. A disagreement in D means those two have come apart,
  ## which is audit_da_runs.R's business, not something to fix by re-filtering.
  ##
  ## N is the one the subsetting cannot reach, so it is the one this check is
  ## really for -- and it is an ERROR, not a warning. A panel on the wrong sample
  ## set is not a panel worth emitting: every LFC in it is a different quantity,
  ## and the fix is always in the prep, never downstream.
  ed <- expected_dims(p$key, p$pair)
  if (!is.na(ed[["N"]]) && ncol(Y) != ed[["N"]]) {
    skipped <- c(skipped, sprintf(
      "%s (built N = %d samples but the analysed run has %d -- this panel's prep is on a different sample set; fix ZEMB_DROP_SAMPLES or the equivalent in its prep)",
      tag, ncol(Y), ed[["N"]]))
    next
  }
  if (!is.na(ed[["D"]]) && nrow(Y) != ed[["D"]]) {
    warning(tag, ": subset to ", nrow(Y), " analysed features but the scale-model ",
            "record says ", ed[["D"]], " -- the saved fit and the record disagree; run ",
            "audit_da_runs.R", call. = FALSE, immediate. = TRUE)
  }
  ## The contrast is read off the group coding rather than typed into the panel
  ## list. Several subtitles had it backwards -- the LFC is always group 1 minus
  ## group 0, and naming the groups the other way inverts how the figure reads.
  contrast <- if (!is.null(r$g1) && !is.null(r$g0))
                paste0(r$g1, " vs ", r$g0) else NA_character_
  sub <- if (is.na(contrast)) p$sub else paste0(p$sub, ", ", contrast)

  lfc <- lfc_vec(Y, X)
  cz  <- clr_zero(Y, X)
  ## The record first, for every study including Jin's pooled run -- so a panel
  ## scores against what `real data results.R` scored against. Only the Jin
  ## subsets, which the record does not carry, fall back to the declared truth.
  tr  <- truth_for(p$key, p$pair)
  if (!length(tr) && identical(p$key, "Jin2022"))
    tr <- jin_truth_for(p$jin[["l"]], p$jin[["o"]])
  ## The measured theta^perp. Where the theta^perp record carries one it wins,
  ## so this figure and the scale-recovery table report the same number; the
  ## studies the record does not cover -- Stammler and Prasse -- are computed
  ## here from their own scale measurement, on whatever basis that study reports.
  ## selex has no total-abundance measurement at all and gets no line.
  ob_rec <- observed_perp(p$key, p$pair)
  ob_own <- if (is.null(r$scale) || anyNA(r$scale) || any(r$scale <= 0)) NA_real_
            else mean(log2(r$scale[X == 1])) - mean(log2(r$scale[X == 0]))
  if (is.finite(ob_rec) && is.finite(ob_own) && abs(ob_rec - ob_own) > 0.01) {
    warning(tag, ": measured theta^perp is ", sprintf("%+.3f", ob_rec),
            " in the theta^perp record but ", sprintf("%+.3f", ob_own),
            " recomputed from this panel's samples -- the two are not describing ",
            "the same sample set or the same scale basis", call. = FALSE,
            immediate. = TRUE)
  }
  ob  <- if (is.finite(ob_rec)) ob_rec else ob_own
  bd  <- observed_perp_band(p$key, p$pair)
  ## The method's own theta^perp, from the saved fit for this exact subset.
  es  <- estimated_perp(p$key, p$pair)
  m   <- lfc_modes(lfc, rownames(Y), tr, clr = cz, obs = ob, est = es)
  ## The drawn quantities lead; the three no longer on the figure follow in
  ## brackets, because they are still the numbers that say whether the panel and
  ## the fit agree even though the figure has stopped asserting it.
  message(sprintf("  %-22s D=%4d N=%3d truth=%2d  mode(null)=%s mode(all)=%+.3f gap=%s  [clr0=%+.3f  est theta^perp=%s  obs theta^perp=%s]",
                  tag, nrow(Y), ncol(Y), length(tr),
                  if (is.finite(m[["null"]])) sprintf("%+.3f", m[["null"]]) else "  n/a",
                  m[["all"]],
                  if (is.finite(m[["null"]])) sprintf("%+.3f", m[["all"]] - m[["null"]]) else "  n/a",
                  m[["clr"]],
                  if (is.finite(es)) sprintf("%+.3f", es) else "n/a",
                  if (is.finite(ob)) sprintf("%+.3f%s", ob,
                        if (is.finite(ob_rec)) "" else " [from study scale]") else "n/a"))
  d <- data.frame(value = as.numeric(lfc), feature = rownames(Y), stringsAsFactors = FALSE)
  pl <- lfc_panel(d, m, p$label, sub, band = bd)
  ggsave(file.path(DIR_FIG, paste0("lfc_diagnose_", tag, ".pdf")), pl,
         width = 6, height = 3.2)
  built[[tag]] <- list(data = d, modes = m, label = p$label, sub = sub,
                       key = p$key, pair = p$pair, band = bd, plot = pl,
                       obs = ob, from_record = is.finite(ob_rec),
                       basis = if (is.null(r$scale_basis)) NA_character_ else r$scale_basis)
}

if (length(skipped)) {
  message("\nNot built:"); for (s in skipped) message("  ", s)
}

## What each measured line actually measures. The studies do not report total
## load the same way -- most give a DENSITY (per mg, per gram), Stammler gives
## copies per tube because it reports no mass to divide by. A per-tube total
## confounds load with how much material went into the tube, which is precisely
## why Zemb is read as a density here. Printed rather than put on the panels,
## which have no room, but it belongs in the figure caption.
bas <- Filter(function(b) is.finite(b$obs), built)
if (length(bas)) {
  message("\nMeasured theta^perp, and what it is measured on:")
  for (k in names(bas)) {
    message(sprintf("  %-22s %+.3f  %-42s %s", k, bas[[k]]$obs,
                    if (is.na(bas[[k]]$basis)) "(basis not recorded)" else bas[[k]]$basis,
                    if (bas[[k]]$from_record) "[theta^perp record]" else "[study scale]"))
  }
  if (any(grepl("per-tube", vapply(bas, function(b) paste(b$basis), character(1))))) {
    message("  NOTE: not every line is on the same basis -- see the per-tube entries above.")
  }
}

## The combined sheet. Faceted rather than assembled from the individual plots,
## so the panels share a theme and the free x scale is the only thing that
## varies -- the LFC ranges differ by an order of magnitude across studies, so a
## common scale would flatten every panel but McMurrough's.
if (length(built)) {
  panel_lab <- vapply(built, function(b) paste0(b$label, "\n", b$sub), character(1))
  dd <- do.call(rbind, lapply(names(built), function(k)
    transform(built[[k]]$data[, "value", drop = FALSE], panel = panel_lab[[k]])))
  vl <- do.call(rbind, lapply(names(built), function(k) {
    m <- built[[k]]$modes
    data.frame(panel = panel_lab[[k]],
               x     = unname(m[LINE_KEYS]),
               which = unname(LINE_LAB[LINE_KEYS]),
               stringsAsFactors = FALSE)
  }))
  vl <- vl[is.finite(vl$x), ]
  dd$panel <- factor(dd$panel, levels = unname(panel_lab))
  vl$panel <- factor(vl$panel, levels = unname(panel_lab))

  ## The measurement band went with the observed theta^perp line it sat behind.
  combined <- ggplot(dd, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), colour = "grey30", fill = "white",
                   bins = 60, linewidth = 0.15) +
    geom_density(alpha = 0.2, bw = 0.5, fill = "#FF6666", colour = "grey40") +
    ## Zero comes through `vl` with the other two, so it is in the legend and
    ## cannot drift out of step with the single panels.
    geom_vline(data = vl, aes(xintercept = x, colour = which, linetype = which),
               linewidth = 0.7, show.legend = TRUE) +
    scale_colour_manual(values = setNames(unname(LINE_COL[LINE_KEYS]),
                                          unname(LINE_LAB[LINE_KEYS])),
                        breaks = unname(LINE_LAB[LINE_KEYS]), name = NULL) +
    scale_linetype_manual(values = setNames(unname(LINE_LTY[LINE_KEYS]),
                                            unname(LINE_LAB[LINE_KEYS])),
                          breaks = unname(LINE_LAB[LINE_KEYS]), name = NULL) +
    facet_wrap(~ panel, scales = "free", ncol = 3) +
    labs(x = expression("LFC (" * hat(theta)^"||" * ")"), y = "Density") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", strip.text = element_text(size = 7.5))
  ## Height follows the row count rather than being fixed at the twelve-panel
  ## sheet's 8in: with JIN_SUBSET_PANELS off there are six panels in two rows, and
  ## a fixed height stretches them.
  n_rows <- ceiling(length(built) / 3)
  ggsave(file.path(DIR_FIG, "lfc_diagnose_all.pdf"), combined,
         width = 10, height = 1.0 + 2.7 * n_rows)
  message("Wrote combined sheet: ", file.path(DIR_FIG, "lfc_diagnose_all.pdf"))
}

message("\nWrote ", length(built), " panels to ", DIR_FIG)
