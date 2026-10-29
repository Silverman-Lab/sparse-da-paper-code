################################################################################
## FILTER RULE: UNCHANGED -- this study is one of the two exceptions.
##                                                     (this run: 09012026)
##
## The project reverted to rowMeans(raw counts) >= 5 on 2026-09-01 for Prasse,
## Stammler, Zemb and selex. Jin keeps its own ">= 5 raw counts in >= 3 samples", recomputed per subset,
## because that is the source paper's own rule and reproducing its feature set is
## part of what this script is for. Nothing about the filter changed here.
##
## The datestamp moved so this study's outputs belong to the same round of fits
## as the other four: `real data results.R` and task_common.R resolve each study
## to the newest stamp that exists, and warn when a table mixes stamps. Re-running
## this file under 09012026 keeps that check quiet for the right reason.
##
## The per-method filter arguments (ANCOM-BC2 prv_cut, LinDA prev.cut, and both
## library cutoffs) are 0, so all twelve methods see the feature set this script
## decided. See count_input_policy.R.
##
## The previous version of this file is archived unchanged at
##   archive/preprocessing_2026-08/real_data/
################################################################################

################################################################################
## PATHS -- edit ROOT here; nothing below this block hardcodes an absolute path.
##
## ROOT is the project directory that used to be spelled out in full on every
## line. Every path is now built from it with file.path(), so moving between the
## cluster, a laptop, or a copy of the project is a one-line change.
################################################################################
#ROOT <- "/storage/work/wpg5129/sparse_project"   # HPC
ROOT <- local({  # manuscript_code: the folder holding .manuscript_code_root (README.md)
  p <- normalizePath(getwd(), mustWork = FALSE)
  while (!file.exists(file.path(p, ".manuscript_code_root"))) {
    if (dirname(p) == p) stop("run from inside the manuscript_code folder (see README.md)")
    p <- dirname(p)
  }
  p
})

## Expand "~" ONCE, here. R's own file functions expand it themselves, but the
## MiMIX paths get handed to the shell through shQuote(), which single-quotes
## them -- and a quoted '~/...' is never expanded by the shell, so julia is asked
## for a literal "~/..." and dies. Absolute cluster ROOTs are unaffected.
ROOT <- path.expand(ROOT)

## R library locations -- a DIFFERENT root from ROOT (home, not work).
## .libPaths() REPLACES the search path rather than appending, so handing it
## HPC-only directories off-cluster would silently discard the local user library
## and packages would appear to vanish. Non-existent entries are dropped first and
## the rest prepended: HPC lookup priority is unchanged, and off-cluster this is a
## no-op instead of a foot-gun.
# R_LIBS <- c("/storage/home/wpg5129/work/R/x86_64-pc-linux-gnu-library/4.4")
# R_LIBS <- R_LIBS[dir.exists(R_LIBS)]
# if (length(R_LIBS)) .libPaths(c(R_LIBS, .libPaths()))

## MiMIX Julia source. NOT just ROOT + a subdir: the cluster keeps this under
## Simulation_SparseDOSSA/ while the git repo calls the same directory simulation/,
## so it needs its own knob rather than riding on ROOT.
## A wrong path here does NOT stop the script -- system() just returns non-zero
## and the MiMIX block goes on to read whatever beta.tsv was left over from an
## earlier run. Pick whichever layout is present and fail loudly if neither is.
MIMIX_SRC <- file.path(ROOT, "Simulation_SparseDOSSA/src/mimix")   # layout on the cluster
if (!dir.exists(MIMIX_SRC)) MIMIX_SRC <- file.path(ROOT, "simulation/src/mimix")   # layout in the git repo
stopifnot(file.exists(file.path(MIMIX_SRC, "scripts/fit-mcmc.jl")))

## Julia binary. system() runs the command through sh, which does NOT read your
## ~/.zshrc -- so a julia that works fine in the Terminal is invisible to R, and
## you get "sh: julia: command not found". This bites hardest when R/RStudio is
## launched from Finder rather than from a shell. Resolve an absolute path once
## here instead of trusting PATH.
JULIA <- Sys.which("julia")
if (!nzchar(JULIA)) {
  cand <- path.expand(c("~/.juliaup/bin/julia",
                        "/usr/local/bin/julia",
                        "/opt/homebrew/bin/julia"))
  cand <- c(cand, Sys.glob("/Applications/Julia-*.app/Contents/Resources/julia/bin/julia"))
  cand <- cand[file.exists(cand)]
  if (length(cand)) JULIA <- cand[1]
}
if (!nzchar(JULIA)) {
  stop("julia not found on PATH or in the usual install locations. ",
       "Set JULIA <- \"/full/path/to/julia\" by hand.")
}


## Loading libraries
library(phyloseq)
library(driver)
library(plyr)
library(dplyr)
library(tidyverse)
library(ALDEx2)
library(gghighlight)
library(cowplot)
library(ggplot2)
library(gghighlight)
library(magrittr)
library(stringi)
library(phyloseq)
library(data.table)
library(DESeq2)
library(limma)
library(edgeR)
library(baySeq) # not yet
library(latex2exp)
library(ggrepel)
library(ANCOMBC)
library(LinDA)
library(numDeriv)
library(dirmult)
source(file.path(ROOT, "sparse_SSRV.R"))
## record_scale_model(): writes which scale model was selected and how many
## taxa EACH arm would have called, to <dir>/scale_model_record.rds.
source(file.path(ROOT, "scale_model_record.R"))
## Which matrix each method is fed (raw counts by default) and the DESeq2
## size-factor fallback that raw counts can need. One file, so the six study
## scripts cannot drift apart on it.
source(file.path(ROOT, "count_input_policy.R"))
## run_ckf(): CKF (the MRCKF knockoff filter) fitted HERE, on the same
## count matrix and contrast sparse_auto() gets, and merged into res_tpfp
## before it is written -- so the summary script reads a "CKF" row like any
## other method instead of loading a separate file and rbind-ing it on.
source(file.path(ROOT, "ckf_record.R"))

################################################################################
## CKF replication settings. See ckf_record.R for what each one actually buys.
##
##   CKF_R          replicates of the knockoff filter. 100 is what the six
##                  standalone scripts used, and is kept as the default so these
##                  numbers stay comparable with the ones already reported.
##   CKF_SE_TARGET  stop early once the Monte Carlo SE of BOTH mean TP and mean
##                  FP is under this, and at least CKF_R_MIN fits have run. NULL
##                  runs all CKF_R. The SE falls as 1/sqrt(R), so a study that has
##                  settled by 30 fits pays for 30. ckf_stability() prints, from
##                  the fits already done, what R would have been enough.
##   CKF_Q          the filter's FDR target. This is the real lever on the spread:
##                  knockoff.threshold() needs at least 1/q taxa above the
##                  threshold before ANY can be selected, so at q = 0.05 a study
##                  whose signal is one or two taxa can only return the empty set
##                  or a set of twenty-plus. Watch `empty_rate` in the run log --
##                  when it is high, raising q is what stabilises the answer, and
##                  more replicates will not.
##   CKF_OFFSET     0 = knockoff (what MRCKF() reports), 1 = knockoff+, which has
##                  finite-sample FDR control, selects less and varies less.
##   CKF_DERAND     derandomisation threshold: taxa selected in at least this
##                  fraction of the replicates form ONE stable selection set,
##                  reported alongside the average. Costs nothing extra.
################################################################################
CKF_R         <- 100
CKF_Q         <- 0.05
CKF_OFFSET    <- 0
CKF_SE_TARGET <- NULL
CKF_R_MIN     <- 20
CKF_DERAND    <- 0.5

set.seed(1234)

## Helper function--------------------------------------------------------------
aldex.lfc <- function(clr){
  # Use clr conditions slot instead of input
  if (is.vector(clr@conds)) {
    conditions <- clr@conds
  } else if (is.factor(clr@conds)) {
    if (length(levels(clr@conds) == 2)) {
      conditions <- clr@conds
    }
  } else if (is.matrix(clr@conds)){
    stop("currently does not support > 2 conditions.")
  } else {
    stop("please check that the conditions parameter for aldex.clr is correct.")
  }
  
  nr <- numFeatures(clr) # number of features
  rn <- getFeatureNames(clr) # feature names
  mc.s <- numMCInstances(clr)
  p <- length(conditions)
  # ---------------------------------------------------------------------
  
  # sanity check to ensure only two conditons passed to this function
  conditions <- as.factor( conditions )
  levels     <- levels( conditions )
  
  sets <- levels
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])
  
  
  ## we need to fill an array and calculate these statistics over the array
  W.est <- array(NA, dim = c(nr, p, mc.s))
  
  for(i in 1:p){
    W.est[,i,] <- getMonteCarloReplicate(clr, i)
  }
  print(paste0("Condition A is ", sets[1], " and condition B is ", sets[2], "."))
  lfc <- rep(NA, mc.s)
  lfc <- apply(W.est, MARGIN = 3, FUN = function(mat, setA, setB){rowMeans(mat[,setA]) - rowMeans(mat[,setB])}, setA =setA, setB = setB)
  return(data.frame("lfc" = rowMeans(lfc), "sd" = apply(lfc,1,sd), "p2.5" = apply(lfc,1,quantile, probs = c(0.025)), "p97.5" = apply(lfc,1,quantile, probs = c(.975))))
}

summarize_tpfp <- function(tax_gold, list_tax_methods, list_tax_methods_abs = NULL){
  
  #list_tax_methods_abs : results using absolute counts. 
  if(!is.null(list_tax_methods_abs)){
    aldex_gold <- matrix(NA, nrow = length(list_tax_methods), ncol = 2)
    rownames(aldex_gold) <- names(list_tax_methods)
    colnames(aldex_gold) <- c("TP", "FP")
    for(i in 1:length(list_tax_methods)){
      tax_method <- list_tax_methods[[i]]
      aldex_gold[i,] <- c(sum(tax_method[,"category"] %in% tax_gold$category),
                          sum(!(tax_method[,"category"] %in% tax_gold$category)))
    }
    
    # methods_gold ; take its own results based on absolute counts as gold standard 
    # and compute TP, FP using its result based on relative counts
    methods_gold <- matrix(NA, nrow = length(list_tax_methods), ncol = 2)
    rownames(methods_gold) <- names(list_tax_methods)
    colnames(methods_gold) <- c("TP", "FP")
    for(i in 1:length(list_tax_methods)){
      tax_method <- list_tax_methods[[i]]
      tax_method_abs <- list_tax_methods_abs[[i]]
      methods_gold[i,] <- c(sum(tax_method[,"category"] %in% tax_method_abs$category),
                            sum(!(tax_method[,"category"] %in% tax_method_abs$category)))
    }
    return(list(aldex_gold = aldex_gold, methods_gold = methods_gold))
    
  }else{
    
    aldex_gold <- matrix(NA, nrow = length(list_tax_methods), ncol = 2)
    rownames(aldex_gold) <- names(list_tax_methods)
    colnames(aldex_gold) <- c("TP", "FP")
    for(i in 1:length(list_tax_methods)){
      tax_method <- list_tax_methods[[i]]
      aldex_gold[i,] <- c(sum(tax_method[,"category"] %in% tax_gold$category),
                          sum(!(tax_method[,"category"] %in% tax_gold$category)))
    }
    return(list(aldex_gold = aldex_gold))
  }
}

## Feature-level bookkeeping ----------------------------------------------------
## lineage_of(), clean_taxonomy(), lowest_assigned(), sig_features_long() and
## sig_features_matrix() are shared with the Prasse analysis, so they live in one
## file rather than being copy-pasted into each script and drifting apart.
source(file.path(ROOT, "feature_record.R"))

## Which column carries each method's own test statistic. Dataset-specific, so it
## stays here rather than in the shared file: the method names and the contrast
## column name (ANCOM-BC's "dietVS") are particular to this comparison.
## An explicit NA means the method reports no effect size at all -- ANCOM gives a
## W count, and baySeq a direction string plus a log posterior, neither of which
## is a signed effect. These numbers are NOT sign-comparable across methods:
## note ANCOM-BC here reports VS-vs-VD while limma/edgeR below use VD - VS.
## See feature_record.R; use lfc_prop / lfc_abs for direction.
EFFECT_COL <- c(
  ALDEx2_original = "effect",       # ALDEx2 standardised effect, NOT a log fold change
  ALDEx2_def      = "effect",
  sparse_SSRV     = "mean",         # posterior mean of theta
  ## The same posterior mean, from the laplace-only and bootstrap-only scale
  ## models. These three ARE directly comparable to one another -- they differ
  ## only in the scale model -- unlike the cross-method statistics.
  sparse_SSRV_L   = "mean",
  sparse_SSRV_B   = "mean",
  DESeq2          = "log2FoldChange",
  edgeR           = "logFC",
  limma           = "logFC",
  baySeq          = NA_character_,
  MiMIX           = "mean",         # posterior mean of beta
  ANCOMBC         = "dietVS",
  ANCOM           = NA_character_,
  ANCOMBC2        = "lfc_dietVS",
  LinDA           = "log2FoldChange"
)

parse_2024_jin_natureComm_technicalReplicates <- function(paths = NULL) {
  local <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates/")
  
  dat <- readxl::read_xlsx(paste0(local, "Supplementary Data 8.xlsx"), sheet = "Sequencing-determined counts", skip = 1)
  counts <- as.matrix(dat[,2:56])
  rownames(counts) <- dat$`cOTU-ID`
  tax <- dat[,57:62]
  rownames(tax) <- dat$`cOTU-ID`
  scale <- readxl::read_xlsx(paste0(local, "Supplementary Data 8.xlsx"), sheet = "Absolute total abundance", skip = 1)
  abs_counts <- readxl::read_xlsx(paste0(local, "Supplementary Data 8.xlsx"), sheet = "Absolute abundance", skip = 1)
  
  #extract metadata from sample id
  metadata <- as.data.frame(matrix(NA, nrow = nrow(scale), ncol = 5))
  rownames(metadata) <- scale$Sample
  colnames(metadata) <- c("diet", "subject", "location", "sample_type", "technical_replicate")
  
  #first two letter indicate diet type (VD, VS, CE)
  metadata$diet <- substr(scale$Sample, start = 1, stop = 2)
  
  #third letter indicate subject id within group. Note; subject in different group are different mice (i.e., VDa != VSa)
  metadata$subject <- substr(scale$Sample, start = 3, stop = 3)
  
  #sample collected in two different location (distal, proximal)
  metadata$location <- substr(scale$Sample, start = 4, stop = 7)
  
  #sample measured using either filter-residue (cell-sample) and flow-through (ecDNA-sample)
  sampleid <- rownames(metadata)
  metadata$sample_type <- unlist(lapply(strsplit(sampleid, split = "-"), function(x) x[2]))
  metadata$sample_type[is.na(metadata$sample_type)] <- "cell"
  
  #technical replicate if applicable (not every diet group have technical replicates)
  metadata$technical_replicate <- as.numeric(substr(scale$Sample, start = 8, stop = 8))
  
  #remove samples that doesn't match
  counts <- counts[,colnames(counts) %in% rownames(metadata)]
  metadata <- metadata[rownames(metadata) %in% colnames(counts),]
  scale <- scale[scale$Sample %in% colnames(counts),]
  abs_counts <- abs_counts[,colnames(abs_counts) %in% colnames(counts)]
  
  all.equal(colnames(counts), rownames(metadata))

  ## ---------------------------------------------------------------------------
  ## The conventional-sequencing side of the paper: ASV and OTU tables.
  ##
  ## Jin et al. ran the SAME differential abundance analysis three ways -- on
  ## cOTU cell abundances (BarBIQ, above) and on ASV and OTU 16S rRNA GENE
  ## abundances -- and the disagreement between them is the paper's main claim.
  ## Both live in their own supplementary files, sheet "VA group", already
  ## restricted to the 16 VS/VD cell-samples and already named with the same
  ## sample ids as the cOTU table, so no metadata parsing is needed for them.
  ## ---------------------------------------------------------------------------
  asv <- readxl::read_xlsx(paste0(local, "Supplementary Data 5.xlsx"), sheet = "VA group", skip = 1)
  otu <- readxl::read_xlsx(paste0(local, "Supplementary Data 7.xlsx"), sheet = "VA group", skip = 1)

  as_count_mat <- function(d, idcol) {
    m <- as.matrix(d[, setdiff(colnames(d), idcol), drop = FALSE])
    mode(m) <- "numeric"
    m[is.na(m)] <- 0
    rownames(m) <- d[[idcol]]
    m
  }
  asv_counts <- as_count_mat(asv, "ASV-ID")
  otu_counts <- as_count_mat(otu, "OTU-ID")

  ## Neither table carries a taxonomy: the GreenGenes / RDP / Silva columns in
  ## Supplementary Data 3 and 4 are percent-identity scores, not lineages. What
  ## those files DO carry is a crosswalk -- every ASV and every OTU names the
  ## cOTU whose Bar sequence it matched -- so each feature can inherit the
  ## lineage of its cOTU, and, more usefully, be reported next to it. That
  ## crosswalk is the paper's own evidence for the splitting it describes:
  ## COTU-CM-2002 resolves to ASV-CM-28 + ASV-CM-31 and to OTU-CM-029.
  asv_ref <- readxl::read_xlsx(paste0(local, "Supplementary Data 3.xlsx"), sheet = "VA group", skip = 1)
  otu_ref <- readxl::read_xlsx(paste0(local, "Supplementary Data 4.xlsx"), sheet = "VA group", skip = 1)
  asv_map <- data.frame(feature = asv_ref$`ASV-ID`, cOTU = asv_ref$`cOTU-ID`, stringsAsFactors = FALSE)
  otu_map <- data.frame(feature = otu_ref$`OTU-ID`, cOTU = otu_ref$`cOTU-ID`, stringsAsFactors = FALSE)
  stopifnot(setequal(asv_map$feature, rownames(asv_counts)),
            setequal(otu_map$feature, rownames(otu_counts)))

  return(list(scale=scale, metadata=metadata, counts=counts, abs_counts = abs_counts, tax=tax,
              asv_counts = asv_counts, otu_counts = otu_counts,
              asv_map = asv_map, otu_map = otu_map))
}

#-------------------------------------------------------------------------------

################################################################################
## ANALYSIS DESIGN -- VS vs VD, run SEPARATELY for each cecal location
##
## Until 2026-08-19 this script ran ONE contrast: VS vs VD over the cecal
## cell-fraction samples, pooling the proximal and distal samples of each mouse
## into a single group of 8 vs 8. That pooling buys sample size and spends the
## regional structure the Jin design was built to expose (Fig 1b: the cecum was
## sampled at two locations precisely because its community is not homogeneous).
##
## The design runs the SAME contrast three ways: once pooled over both
## locations, and once per location.
##
##   both : all sixteen cell samples          -- 8 VD vs 8 VS   (PRIMARY)
##   prox : VDaprox VDbprox VDcprox VDdprox   vs   VSaprox VSbprox VScprox VSdprox
##   dist : VDadist VDbdist VDcdist VDddist   vs   VSadist VSbdist VScdist VSddist
##
## Which of the three is the headline result, as of 2026-08-28:
##
##   PRIMARY       the pooled run. It is the subset reported in the main DA
##                 results table (`real data results.R`), because it is the only
##                 one whose measured-scale gold reference is non-empty -- at
##                 4 vs 4 both stratified references are empty at padj<0.05, and
##                 TP/FP against an empty reference is not a method comparison.
##   SUPPLEMENTARY the two stratified runs. They are what the manuscript's Jin
##                 tables are built from (`manuscript_tables.R`, via
##                 `jin_DA_biological_interpretation.R`, which enumerates
##                 LOC_KEYS = prox/dist and never sees the pooled run).
##
## The pooled run is NOT a third location and must not be scored as one. It is
## the same eight mice a second time, so the cross-location section at the bottom
## of this file excludes it -- see the note there.
##
## Each subset is analysed END TO END on its own: cOTUs, filtered by the paper's
## own prevalence rule within that subset (so D differs between subsets -- 150
## prox, 153 dist, 189 pooled), its own gold-standard scale model built from
## the measured total abundance of ITS OWN samples, and therefore its own
## reference set. Nothing is shared across the runs except the parsed input.
##
## The filter is the SAME RULE everywhere -- >=MIN_COUNT raw counts in
## >=MIN_SAMPLES samples of the subset -- recomputed per subset. It is deliberately
## not a shared feature set: the pooled run tests more cOTUs than either stratified
## run because "3 samples" is a weaker demand on 16 samples than on 8. A cOTU can
## therefore be tested in the primary analysis and absent from a supplementary one;
## `tested_in` in jin_sig_by_location_*.csv is what tells that apart from a
## biological absence.
##
## KNOWN CAVEAT, measured rather than assumed. The rule is usually read as
## ">=MIN_SAMPLES biologically independent samples", and for prox and dist it is
## exactly that -- each holds one sample per mouse. Pooled does not: its 16 samples
## are 8 mice twice over, so "3 of 16" can be met by as few as 2 mice. On this data
## 21 of the 189 pooled cOTUs are backed by fewer than three mice, and one of the
## two cOTUs the pooled sparse SSRV fit calls (COTU-CM-0624) is not among them, so
## no current call rests on it. Two stricter forms were measured and NOT adopted:
## requiring the rule in each location gives 138 cOTUs, in either location 165
## (the union of the stratified sets); both leave the gold reference, ALDEx2_def
## and sparse SSRV's four calls unchanged. Kept as the plain rule for consistency
## with the stratified runs and with the source paper.
##
## Two things to keep in mind when writing this up:
##
##   1. The two runs are NOT independent replicates. The same four VD mice and
##      the same four VS mice contribute to both, once at each location, so a
##      taxon called at one location and not the other is a statement about
##      WITHIN-mouse regional structure, not about independent confirmation.
##      (The VS/VD arms carry no technical replicates -- those exist only for the
##      CE mice, which this contrast drops -- so each mouse-location cell is a
##      single sample.)
##   2. Each run is 4 vs 4. Power is lower than the pooled 8 vs 8, so a taxon
##      that drops out at one location has not been shown to be absent there;
##      read the lfc_abs column alongside the call, not the call alone.
##
## File naming. The stratified runs are written with a "_<loc>_" tag; the pooled
## cOTU run keeps the BARE names (tax_scale_09012026.RData,
## relative_results_09012026.RData, res_tpfp_09012026.RData, sparse_sens_09012026.RData)
## that this script wrote before the split, which is what `real data results.R`
## reads with JIN_TAG <- "". So the primary analysis is addressed by the same
## file names it always was, and the supplementary runs sit beside it under
## their own tags.
##
## Records (scale_model_record.rds, sparse_arm_sig.rds, theta_perp_record.rds,
## ckf_record.rds) key the pooled cOTU run as pair = "both". A pair = NA Jin row
## in any of those is a leftover from the pre-split script and should be dropped;
## `real data results.R` warns about exactly that case.
################################################################################

DIR_JIN   <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates")
DATESTAMP <- "09012026"

## Stratified subsets -- the SUPPLEMENTARY analysis. "prox" / "dist" keep only
## that location. These are the subsets jin_DA_biological_interpretation.R
## enumerates (its LOC_KEYS), so every manuscript Jin table is built from them.
LOCATIONS <- c("prox", "dist")

## Levels that ALSO get the pooled ("both") run -- the PRIMARY analysis.
##
## cOTU only, and that is not an oversight. The pooled run exists to feed the
## main DA results table, which scores methods against the measured-scale gold
## reference; that reference needs the BarBIQ total abundance, which was measured
## at cOTU level only (Supplementary Data 8). A pooled ASV or OTU run would carry
## no gold reference, would not enter the main table, and would not enter the
## supplementary tables either -- it would just be two more multi-hour fits
## (MiMIX, baySeq, ANCOM, sparse_SSRV at S draws) writing files nothing reads.
##
## Set to LEVELS to run the pooled contrast at every representation anyway.
POOLED_LEVELS <- "cOTU"

## Feature types to analyse. The paper ran the same VS-vs-VD contrast three
## ways and its headline claim is that they disagree:
##
##   cOTU : cell abundance, BarBIQ -- the number of CELLS per taxon
##   ASV  : 16S rRNA GENE abundance, conventional amplicon sequencing
##   OTU  : 16S rRNA GENE abundance, conventional amplicon sequencing
##
## Their per-location feature counts after the >=5-in->=3 filter were 150/153
## cOTUs, 236/223 ASVs and 138/131 OTUs (prox/dist); this script reproduces all
## six, which is the check that the filter is being read the way they read it.
##
## Only cOTU has a measured scale. Supplementary Data 8 carries the absolute
## total abundance per sample, Data 5 and Data 7 carry counts alone -- so the
## scale-aware gold standard, lfc_abs and the TP/FP scoring exist for cOTU only
## and are skipped for ASV and OTU. For those two the question is not "who
## matches the reference" but "who calls what", which is what the sig_matrix
## and the cross-level table at the bottom report.
LEVELS <- c("cOTU", "ASV", "OTU")

## Significance rule for the summary tables -------------------------------------
## Jin et al. do not call a feature significant on FDR alone. From the Fig. 7c
## legend: "FDR < 0.05 and a fold change > 2 (gray dotted lines) were considered
## 'significant differences' in other analyses", and the Fig. 7a/b legends repeat
## it ("Different, FDR < 0.05 and fold change > 2"). Both of their headline
## counts are under this rule -- 8 cOTUs proximal / 1 distal, 8 and 6 ASVs, 9
## and 6 OTUs.
##
## "fold change > 2" is a magnitude in either direction, so the test is
## |log2 FC| > 1. Their fold change is DESeq2's log2FoldChange computed on
## RELATIVE abundances (Methods: counts normalized by median-of-ratios), which
## is why the default basis below is lfc_prop and not lfc_abs -- lfc_abs adds the
## measured total load back and is a different quantity, larger by 0.86 log2
## proximally.
SIG_FDR <- 0.05
SIG_FC  <- 2          # fold change > 2  <=>  |log2 fold change| > 1

## Which log2 fold change the >2 test is applied to:
##   "lfc_prop" (default) the shared compositional log2 FC this script computes
##              once from the counts under a single VD-minus-VS convention. Same
##              quantity for every method, which is what makes the per-method
##              counts comparable, and the closest match to the paper.
##   "lfc_abs"  the same shifted onto the absolute scale. cOTU only -- NA at
##              ASV/OTU level, where no total load was measured.
##   "method"   each method's OWN reported effect size, falling back to lfc_prop
##              where it reports none. Faithful to the paper for DESeq2, but NOT
##              comparable across methods: ALDEx2's "effect" is a standardised
##              effect and not a log fold change at all, MiMIX's beta is on the
##              multinomial-logit scale, ANCOM reports only a W count and baySeq
##              no effect size. Offered for checking DESeq2 against the paper,
##              not for the cross-method table.
FC_BASIS <- "lfc_prop"

## Abundance filter -- the original authors' own rule ---------------------------
## Jin et al., Methods, "Differential abundance analysis":
##
##   "For the differential abundance analyses based on relative abundances, the
##    cOTUs, ASVs, or OTUs with raw counts >= 5 in at least three samples in each
##    pair of groups were used"
##
## This replaces the mean-count threshold the script used to apply. The two are
## not variations on one idea: a mean is dominated by a single large sample, so
## `mean(count) > 5` kept genera carried by one mouse and dropped genera present
## at a steady count of 4 in every sample. The paper's rule is a PREVALENCE rule
## at a count floor, which is what a 4-vs-4 contrast actually needs.
MIN_COUNT   <- 5   # "raw counts >= 5"
MIN_SAMPLES <- 3   # "in at least three samples"

## What "in each pair of groups" quantifies over. The sentence names the PAIR --
## i.e. the sample set of the comparison being run -- not the two groups within
## it, and the paper ran several such pairs (CE/VS, VS/VD, ...), which is what
## the phrase is there to distinguish. So the default counts the three samples
## anywhere among the eight in the contrast.
##
## Setting this TRUE demands three qualifying samples in EACH diet arm
## separately. That is a defensible stricter reading of the same sentence, but
## on this data it is a much heavier filter -- at 4 vs 4 it requires a genus in
## 3 of 4 mice on BOTH diets, which removes roughly half the features (prox
## 135 -> 69, dist 140 -> 70) and, by construction, most of the diet-specific
## ones. Left FALSE to follow the paper as written.
FILTER_PER_GROUP <- FALSE

## Feature level. The paper ran DESeq2 on cOTUs directly -- it never agglomerated
## -- so FALSE (the default) reproduces its unit of analysis and its feature
## count: 150 prox / 153 dist of the 860 cOTUs, against 135 / 140 when the same
## filter is applied to genus sums.
##
## Dropping the agglomeration is not only about matching the paper. This
## taxonomy carries the classifier's bootstrap confidence inside every string
## ("Barnesiella(0.52)", "Barnesiella(0.71)"), and phyloseq merges taxa whose
## taxonomy strings match exactly. The 30-odd Barnesiella cOTUs here carry 25
## different confidence values, so they were never being pooled into one genus
## in the first place: tax_glom took 860 cOTUs to 645 groups, which is exactly
## the number of distinct full lineage strings. What it merged was
## confidence-value duplicates, not genera. Agglomerating on those strings gave
## neither cOTU resolution nor genus resolution.
##
## Set TRUE to restore the genus-level agglomeration the other real-data scripts
## in this project use. Note the two are scored differently -- lineage_of() and
## the sig_matrix are keyed on whatever ids survive here -- so results across
## the two settings are not directly comparable.
AGGLOMERATE_TO_GENUS <- FALSE

## Reading in data --------------------------------------------------------------
## Parsed ONCE, outside the loop: the spreadsheet is the same for both subsets
## and re-reading it per location only invites the two runs to disagree about
## their input.

out <- parse_2024_jin_natureComm_technicalReplicates()

saveRDS(out, file.path(DIR_JIN, "output_2024_jin_natureComm_technicalReplicates.RData"))

OTU_all <- t(out$counts)
taxa    <- out$tax
## Metadata reordered to the count-table sample order, once. Every mask below is
## built from THIS frame, so masks and count columns always line up.
metadata_all <- out$metadata[match(rownames(OTU_all), rownames(out$metadata)), ]
stopifnot(identical(rownames(OTU_all), rownames(metadata_all)))

## Look the totals up BY SAMPLE NAME rather than by position. out$scale is in the
## order the spreadsheet used, which interleaves the -ecDNA rows, so applying a
## count-table-order mask to spreadsheet-order rows hands each sample a DIFFERENT
## sample's total abundance -- the bug that used to flatten the group difference
## in mean log2 total load from 0.555 to 0.019. Resolving the lookup once here,
## for ALL samples, means the per-location subsets below are just row selections
## of an already-correct vector and cannot reintroduce it.
scale_ind_all <- match(rownames(metadata_all), out$scale$Sample)
stopifnot(!anyNA(scale_ind_all))
scale_all     <- out$scale$`Absolute total abudance (copies/mg)`[scale_ind_all]
scale_std_all <- out$scale$`Standard deviation (N  = 4 or 5)`[scale_ind_all]
stopifnot(length(scale_all) == nrow(metadata_all),
          length(scale_std_all) == nrow(metadata_all))

## Overview of the measured total load, before any subsetting -----------------
ind.vs <- metadata_all$diet == "VS"
ind.vd <- metadata_all$diet == "VD"
ind.ce <- metadata_all$diet == "CE"

inds_typ1 <- metadata_all$sample_type == "cell"
inds_typ2 <- metadata_all$sample_type == "ecDNA"
inds_loc  <- metadata_all$location == "prox"

scl_rng <- c(min(scale_all, na.rm = TRUE), max(scale_all, na.rm = TRUE))

hist(scale_all[ind.vs], breaks = 10, xlim = scl_rng, main = "total load: VS (all)")
hist(scale_all[ind.vd], breaks = 10, add = TRUE, col = "red", xlim = scl_rng)

hist(scale_all[inds_loc],  breaks = 10, main = "total load: proximal")
hist(scale_all[!inds_loc], breaks = 10, main = "total load: distal")

hist(scale_all[inds_loc & ind.vs],  breaks = 10, main = "total load: proximal")
hist(scale_all[inds_loc & ind.vd],  breaks = 10, add = TRUE, col = "red", xlim = scl_rng)
hist(scale_all[!inds_loc & ind.vs], breaks = 10, main = "total load: distal")
hist(scale_all[!inds_loc & ind.vd],  breaks = 10, add = TRUE, col = "red", xlim = scl_rng)

hist(scale_all[ind.vs & inds_typ1], breaks = 10, xlim = scl_rng, main = "total load: cell fraction")
hist(scale_all[ind.vd & inds_typ1], breaks = 10, add = TRUE, col = "red", xlim = scl_rng)

hist(scale_all[ind.vs & inds_typ2], breaks = 10, xlim = scl_rng, main = "total load: ecDNA fraction")
hist(scale_all[ind.vd & inds_typ2], breaks = 10, add = TRUE, col = "red", xlim = scl_rng)

## The group difference in mean log2 total load, per location. This is the shift
## between the compositional and the absolute scale, and it is what the gold
## standard adds back -- reported here so the two locations can be compared on
## it directly before any model is fitted.
cat("\n== group difference in mean log2 total load (VD - VS), cell fraction ==\n")
for (loc_ in c("prox", "dist")) {
  k <- metadata_all$diet %in% c("VS", "VD") & inds_typ1 & metadata_all$location == loc_
  cat(sprintf("  %-4s : %+0.3f  (VD n=%d, VS n=%d)\n", loc_,
              mean(log2(scale_all[k & ind.vd])) - mean(log2(scale_all[k & ind.vs])),
              sum(k & ind.vd), sum(k & ind.vs)))
}
k <- metadata_all$diet %in% c("VS", "VD") & inds_typ1
cat(sprintf("  %-4s : %+0.3f  (VD n=%d, VS n=%d)\n", "both",
            mean(log2(scale_all[k & ind.vd])) - mean(log2(scale_all[k & ind.vs])),
            sum(k & ind.vd), sum(k & ind.vs)))

#-------------------------------------------------------------------------------

## Per-location results, filled in by the loop below. Keeping them keyed by
## location is what makes the cross-location section at the bottom possible
## without re-running anything.
## Results, keyed res_all[[level]][[location]]. Keeping every run in memory is
## what lets the summary sections below line them up without refitting.
res_all <- list()

for (level in LEVELS) {

## The measured total abundance exists for cOTUs only, so the gold standard,
## lfc_abs and the TP/FP scoring are cOTU-only. Everything else -- the twelve
## methods, the feature record, the per-location split -- runs identically at
## all three levels.
has_scale <- identical(level, "cOTU")

## Count table and taxonomy for this level. ASV and OTU inherit their lineage
## from the cOTU each one matched (see the crosswalk in the parser); the extra
## `cOTU` column is kept so a significant ASV can be read next to the cell-level
## unit it came from, which is the comparison the paper is making.
if (has_scale) {
  COUNTS_ALL <- OTU_all                       # samples x cOTUs
  TAX_ALL    <- as.matrix(taxa)
  XREF       <- NULL
} else {
  cnt_mat    <- if (identical(level, "ASV")) out$asv_counts else out$otu_counts
  XREF       <- if (identical(level, "ASV")) out$asv_map    else out$otu_map
  COUNTS_ALL <- t(cnt_mat)                    # samples x features
  TAX_ALL    <- as.matrix(taxa)[match(XREF$cOTU, rownames(as.matrix(taxa))), , drop = FALSE]
  rownames(TAX_ALL) <- XREF$feature
  TAX_ALL    <- TAX_ALL[rownames(cnt_mat), , drop = FALSE]
}

res_by_loc <- list()

## Subsets for THIS level: the pooled run first where it applies, then the two
## stratified ones. Ordering the pooled run first is deliberate -- it is the
## primary result, and it is what the per-level console summary leads with.
subsets <- c(if (level %in% POOLED_LEVELS) "both", LOCATIONS)

for (loc in subsets) {

  ## File-name tag. Every artefact this iteration writes carries it, so no two
  ## runs can overwrite each other and none can be mistaken for the pooled
  ## result. cOTU keeps the bare "_<loc>_" tag it has always written, so the
  ## files already on disk stay the ones the script maintains.
  key <- if (has_scale) loc else paste0(level, "_", loc)
  tag <- if (identical(loc, "both") && has_scale) paste0("_", DATESTAMP) else paste0("_", key, "_", DATESTAMP)

  message("\n==================================================================")
  message("Jin 2022 -- VS vs VD | level: ", level, " | location: ", loc)
  message("==================================================================")

  ## Subset: VS and VD only (CE dropped), cell fraction only (ecDNA dropped),
  ## and this location only. Built as ONE mask over metadata_all so the count
  ## rows, the metadata rows and the total-load entries are selected together.
  ## The ASV and OTU tables only ever contained the 16 VS/VD cell-samples, so
  ## for them the diet and sample_type clauses are already satisfied -- the mask
  ## is applied through the sample NAMES either way, never by position.
  keep_names <- rownames(metadata_all)[
    metadata_all$diet %in% c("VS", "VD") &
    metadata_all$sample_type == "cell" &
    (if (identical(loc, "both")) TRUE else metadata_all$location == loc)]
  keep_names <- intersect(keep_names, rownames(COUNTS_ALL))
  stopifnot(length(keep_names) > 0)

  OTU      <- COUNTS_ALL[keep_names, , drop = FALSE]
  metadata <- metadata_all[keep_names, , drop = FALSE]
  if (has_scale) {
    scale     <- scale_all[match(keep_names, rownames(metadata_all))]
    scale_std <- scale_std_all[match(keep_names, rownames(metadata_all))]
    ## Named so that perp_matched_obs() below can verify the measurement against
    ## the count table by sample NAME. The match() above is already by name; the
    ## names carry that guarantee forward to the one place that re-pairs the two.
    names(scale) <- keep_names
    stopifnot(length(scale) == nrow(metadata), length(scale_std) == nrow(metadata))
  } else {
    scale <- scale_std <- NULL
  }
  stopifnot(identical(rownames(OTU), rownames(metadata)))

  message("  samples (", nrow(metadata), "): ",
          paste(rownames(metadata), collapse = ", "))

  ##Creating a phyloseq object
  otu1 <- otu_table(OTU, taxa_are_rows = FALSE)
  sam1 <- sample_data(metadata)
  tax1 <- tax_table(TAX_ALL)

  phylo <- phyloseq(otu1, sam1, tax1)

  ##Feature level -- cOTUs by default, matching the paper; see the
  ##AGGLOMERATE_TO_GENUS block above LOCATIONS for why the agglomeration was
  ##doing less than its name suggests.
  phylo_feat <- if (AGGLOMERATE_TO_GENUS) tax_glom(phylo, taxrank = "Genus") else phylo

  ##Filtering OTU table in phyloseq
  ## The paper's rule -- raw counts >= MIN_COUNT in at least MIN_SAMPLES samples;
  ## see the block above LOCATIONS for the quote and for what "each pair of
  ## groups" is taken to mean.
  ##
  ## Applied WITHIN the subset, so each location is filtered on its own samples.
  ## That is the point of running them separately -- a cOTU that clears the
  ## threshold distally but not proximally should not be carried into the
  ## proximal analysis on the strength of the distal samples -- but it does mean
  ## the two runs test slightly different feature sets, which the cross-location
  ## section at the bottom accounts for explicitly. It also matches the paper's
  ## own practice: their filter is defined per pairwise comparison, so the
  ## feature set is a property of the contrast, not of the dataset.
  ##
  ## The counts go in RAW. filter_taxa() sees otu_table(phylo_feat), which is
  ## the un-pseudocounted matrix -- Y <- Y_raw + 1 happens further down, after
  ## pruning -- so ">= 5" really is ">= 5 raw counts" and not ">= 4".
  diet_grp <- as.character(sample_data(phylo_feat)$diet)
  stopifnot(length(unique(diet_grp)) == 2)

  filter.tax <- if (FILTER_PER_GROUP) {
    function(vec) isTRUE(all(tapply(vec, diet_grp,
                                    function(z) sum(z >= MIN_COUNT) >= MIN_SAMPLES)))
  } else {
    function(vec) sum(vec >= MIN_COUNT) >= MIN_SAMPLES
  }
  red_phylo <- filter_taxa(phylo_feat, filter.tax, prune = TRUE)

  message("  features: ", ntaxa(phylo), " ", level, "s",
          if (AGGLOMERATE_TO_GENUS) paste0(" -> ", ntaxa(phylo_feat), " genus groups") else "",
          " -> ", ntaxa(red_phylo),
          " after filter (>=", MIN_COUNT, " raw counts in >=", MIN_SAMPLES,
          if (FILTER_PER_GROUP) " samples of EACH diet arm)" else " of the ",
          nsamples(phylo_feat), " samples in the pair)")
  stopifnot(ntaxa(red_phylo) > 0)

  ##Estimate of theta^\perp from the data
  inds <- which(sample_data(red_phylo)$diet == "VS")

  normalized_GC <- apply(t(otu_table(red_phylo))+.5, 2, FUN = function(x){x/sum(x)})
  gm_GC <- apply(normalized_GC, 2, FUN = function(x){mean(log2(x))})
  cat(sprintf("  theta^perp from the data (%s): %+0.4f\n", loc,
              mean(-1*gm_GC[inds]) - mean(-1*gm_GC[-inds])))

  ## Estimate of theta^\perp from scale -- only definable where a scale was
  ## measured. Data 5 and 7 are count tables with no total-abundance sheet.
  if (has_scale)
    cat(sprintf("  theta^perp from the scale (%s): %+0.4f\n", loc,
                mean(-log2(scale[inds])) - mean(-log2(scale[-inds]))))

  ## Count matrices for this comparison, defined ONCE.
  ##   Y_raw : counts as they come out of phyloseq. THIS is what every model is
  ##           fitted on -- see count_input_policy.R for method-by-method reasons.
  ##   Y     : Y_raw + 1, kept only for derived quantities that want a zero-free
  ##           matrix (log-ratio summaries, prevalence filters written as Y > 1).
  ##           No model is fitted on it.
  ## Until 2026-08-19 every count-based method was fed Y. On the Jin data that
  ## alone was the difference between reproducing the source paper's DESeq2 counts
  ## (8 proximal / 1 distal) and not (9 / 5).
  Y_raw <- t(otu_table(red_phylo))
  Y     <- Y_raw + 1
  X <- sample_data(red_phylo)$diet
  print(dim(Y_raw))
  D = nrow(Y)
  N = ncol(Y)
  propY <- matrix(NA, nrow = D, ncol = N)
  ## Loop over `n`, NOT `i`: `i` is free here but the enclosing loop variable is
  ## `loc`, and keeping the inner counters distinct from anything the outer loop
  ## owns is what stops this pattern from biting later.
  for(n in 1:N){
    propY[,n] <- ((Y[,n]+1)/colSums(Y)[n]) # colSums = sequencing depth
  }
  LFC <- rowMeans(log2(propY[,X=="VD"])) - rowMeans(log2(propY[,X=="VS"]))
  hist(LFC, breaks = 100, main = paste0("Jin ", key, ": VD vs VS relative LFC"))

  #-----------------------------------------------------------------------------

  ## ALDEx2 models--------------------------------------------------------------
  ## First aldex2
  aldex_fit <- aldex(Y_raw,X,mc.samples = 1000)

  tax_clr = aldex_fit %>%
    rownames_to_column("category") %>%
    dplyr::select(category, effect, we.ep, we.eBH) %>%
    mutate(pval = we.ep) %>%
    mutate(padj = we.eBH) %>%
    dplyr::filter(padj < 0.05)


  ##Default scale model
  aldex_fit <- aldex(Y_raw,X,mc.samples = 1000, gamma = .5)

  tax_default = aldex_fit %>%
    rownames_to_column("category") %>%
    dplyr::select(category, effect, we.ep, we.eBH) %>%
    mutate(pval = we.ep) %>%
    mutate(padj = we.eBH) %>%
    dplyr::filter(padj < 0.05)

  ## Gold standard -- cOTU only. The measured absolute total abundance lives in
  ## Supplementary Data 8; the ASV and OTU tables (Data 5, Data 7) carry
  ## sequencing counts and nothing else, so there is no scale to build a
  ## reference from and nothing to score TP/FP against. Those levels report
  ## which features each method called and stop there.
  tax_scale <- NULL
  n_scale_rep <- rep(4, length(scale))
  if (has_scale) {

    ##Gold standard flow cytometry model
    ## Built from the measured per-sample total abundance and its replicate SD --
    ## of THIS location's samples only. This is the "golden reference will also be
    ## different for each subset" requirement: the reference is not transported
    ## from the pooled run or from the other location, it is refitted here.
    scale_samples <- matrix(NA, nrow = ncol(Y), ncol = 3000)
    for(i in 1:ncol(Y)){
      mu <- scale[i]
      sigma <- scale_std[i]/sqrt(n_scale_rep[i])

      sig_log <- sqrt(log(1 + (sigma^2 / mu^2)))
      mu_log <- log(mu) - (0.5 * sig_log^2)

      scale_samples[i, ] <- exp(rnorm(3000, mean = mu_log, sd = sig_log))
    }
    clr <- aldex.clr(Y_raw,X, gamma = scale_samples,mc.samples = 3000)
    res <- aldex.ttest(clr)
    res <- cbind(res, aldex.lfc(clr))

    hist(rowMeans(clr@scaleSamps), main = paste0("Jin ", key, ": mean scale samples"))
    hist(colMeans(clr@scaleSamps[inds,]) - colMeans(clr@scaleSamps[-inds,]),
         main = paste0("Jin ", key, ": VS - VD scale offset"))

    tax_scale = res %>%
      rownames_to_column("category") %>%
      mutate(pval = we.ep) %>%
      mutate(padj = we.eBH) %>%
      dplyr::filter(padj < 0.05)

    ## Persist the reference set for this location. Everything downstream is
    ## scored against it, so it has to be inspectable and re-usable on its own.
    saveRDS(tax_scale, file.path(DIR_JIN, paste0("tax_scale", tag, ".RData")))

    ## How close the reference came, not just how big it is. At 4 vs 4 the BH
    ## correction over ~130 cOTUs is doing most of the work, so the gap between
    ## "one hit" and "no hits" can be a single taxon moving across 0.05 -- and the
    ## reference set size is the one number every TP/FP figure downstream rests on.
    cat(sprintf("  gold (%s): %d cOTUs at padj<0.05 | min padj %.4f | %d at padj<0.10 | %d at raw p<0.05\n",
                loc, nrow(tax_scale), min(res$we.eBH, na.rm = TRUE),
                sum(res$we.eBH < 0.10, na.rm = TRUE), sum(res$we.ep < 0.05, na.rm = TRUE)))

    ## An EMPTY reference is not a failure of the run, but it does make the TP/FP
    ## scoring below degenerate: with no gold taxa every call is counted a false
    ## positive and every method scores 0 TP, which is a statement about the
    ## reference's power at n=4 vs 4, not about the methods. Say so loudly here so
    ## the resulting all-zero TP column is never read as a method comparison.
    if (nrow(tax_scale) == 0)
      warning("location ", key, ": the measured-scale gold reference is EMPTY at ",
              "padj<0.05 -- TP/FP against it is uninformative for this subset; ",
              "compare the two locations on lfc_abs and on the method-vs-method ",
              "tables instead.", call. = FALSE, immediate. = TRUE)

  }

  #-----------------------------------------------------------------------------

  ## Other methods -------------------------------------------------------------

  ## DeSeq2
  ## Conditions taken from red_phylo, the object every other method here is fed.
  ## phylo / phylo_feat have the same samples in the same order, but reading the
  ## design off a different object than the counts is how a subsetting change
  ## turns into a silently mislabelled contrast.
  coldata <- matrix(sample_data(red_phylo)$diet, ncol = 1)
  colnames(coldata) <- "Condition"
  ## Raw counts, with the size-factor fallback described in
  ## count_input_policy.R -- the default median-of-ratios estimator cannot
  ## cope when every feature has a zero somewhere, which is what the old
  ## +1 was masking.
  dds <- deseq_dataset_raw(Y_raw@.Data,
                           colData = coldata,
                           design  = ~Condition)
  dds <- DESeq(dds, betaPrior = TRUE)
  res <- deseq_results(dds)
  tax_deseq <- res %>%
    as.data.frame() %>%
    rownames_to_column("category") %>%
    dplyr::select(category, log2FoldChange, padj, lfcSE) %>%
    mutate(sig = ifelse(padj < 0.05, TRUE, FALSE)) %>%
    dplyr::filter(sig == TRUE)

  ##EdgeR

  X <- sample_data(red_phylo)$diet
  conds <- sample_data(red_phylo)$diet
  coldata <- ifelse(conds == "VD", 0, 1)
  coldata <- matrix(coldata, ncol = 1)

  y <- DGEList(counts=as.matrix(Y_raw@.Data),group=c(coldata[,1]))
  y <- calcNormFactors(y)
  y <- estimateDisp(y)
  et <- exactTest(y)
  ## n = nrow(Y), not a hardcoded ceiling: D is now set by the per-location
  ## filter and differs between runs.
  tax_edgeR <- topTags(et, n=nrow(Y), p.value = 1) %>%
    as.data.frame() %>%
    rownames_to_column("category") %>%
    mutate(sig = ifelse(FDR < 0.05, TRUE, FALSE)) %>%
    dplyr::filter(sig == TRUE)


  ##limma
  X <- sample_data(red_phylo)$diet
  conds <- sample_data(red_phylo)$diet
  coldata <- ifelse(conds == "VD", 0, 1)
  coldata <- matrix(coldata, ncol = 1)
  design <- model.matrix(~0+conds)

  y <- DGEList(counts=as.matrix(Y_raw@.Data),group=conds)
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y,design,plot=FALSE)
  fit <- lmFit(v,design)

  contr <- makeContrasts(condsVD - condsVS, levels = colnames(coef(fit)))
  contr.fit <- contrasts.fit(fit, contr)
  contr.fit <- eBayes(contr.fit)
  ## number = nrow(Y). This read `number=100` while D was 132, so limma was only
  ## ever asked for its top 100 features -- silently capping what it could call.
  tax_limma <- topTable(contr.fit, p.value = 1, number = nrow(Y)) %>%
    rownames_to_column("category") %>%
    mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>%
    dplyr::filter(sig == TRUE)


  ##baySeq
  conds <- sample_data(red_phylo)$diet
  coldata <- ifelse(conds == "VD", 0, 1)
  groups <- list(NDE=rep(1, length(conds)),
                 DE=as.numeric(coldata))
  ## The feature ids have to be handed to baySeq through @annotation. topCounts()
  ## returns its rows SORTED by posterior likelihood, and with an empty @annotation
  ## slot it returns them with no usable rownames at all (they come back as "NA",
  ## "NA.1", "NA.2", ...). The previous `rownames(res) <- taxa_names(red_phylo)`
  ## therefore stamped the UNSORTED taxon order onto SORTED rows, so every taxon
  ## baySeq was credited with was the wrong one, and `cbind(res, post)` -- res
  ## sorted, posteriors in original order -- misaligned in exactly the same way.
  dds <- new("countData", data = as.matrix(Y_raw@.Data), groups = groups, replicates = conds,
             annotation = data.frame(category = rownames(Y), stringsAsFactors = FALSE))
  libsizes(dds) <- getLibsizes(dds)
  dds <- getPriors.NB(dds, cl = NULL)
  dds <- getLikelihoods(dds,cl=NULL,bootStraps=3,verbose=FALSE)
  res <- topCounts(dds,group="DE", number = nrow(Y))
  stopifnot(!is.null(res$category), setequal(as.character(res$category), rownames(Y)))

  ##Merge to the data to see which taxa because the names are uninformative
  post <- as.matrix(dds@posteriors)
  colnames(post) = c("NDE", "DEEst")
  rownames(post) <- rownames(Y)          # posteriors keep the ORIGINAL feature order
  tax_baySeq <- cbind(res, post[as.character(res$category), , drop = FALSE]) %>%
    mutate(sig = ifelse(FWER.DE < 0.05, TRUE, FALSE)) %>%
    dplyr::filter(sig == TRUE)

  ### MiMIX
  ## Its own working directory per location. MiMIX communicates through files in
  ## a fixed directory (Y.csv/X.csv/Z.csv in, beta.tsv out); sharing one
  ## directory between the two runs means the second silently overwrites the
  ## first's inputs and leaves no way to re-inspect either fit.
  mimix_dir <- file.path(DIR_JIN, "MiMIX", key)
  dir.create(mimix_dir, recursive = TRUE, showWarnings = FALSE)

  x <-  as.numeric(X == "VD")
  write.csv(t(Y_raw), file.path(mimix_dir, "Y.csv"), row.names = F)
  write.csv(x,        file.path(mimix_dir, "X.csv"), row.names = F)
  write.csv(1:length(x), file.path(mimix_dir, "Z.csv"), row.names = F)

  src_path <- MIMIX_SRC

  julia_command <- paste(
    shQuote(JULIA), shQuote(file.path(src_path, "scripts/fit-mcmc.jl")),
    "--hyper", shQuote(file.path(src_path, "nutnet-analysis/configs/hyper.yml")),
    "--monitor", shQuote(file.path(src_path, "nutnet-analysis/configs/monitor-mimix-no-factors.yml")),
    "--inits", shQuote(file.path(src_path, "nutnet-analysis/configs/inits.yml")),
    "--factors", "0",
    "--iters", "5000",
    "--chains", "1",
    shQuote(mimix_dir),
    shQuote(mimix_dir)
  )

  ## Delete the previous fit first, then insist Julia actually produced a new one.
  ## Without this, any failure of the julia call (missing binary, bad MIMIX_SRC,
  ## unexpanded "~") leaves the old beta.tsv in place and the lines below report a
  ## stale result as if it were this run's.
  beta_file <- file.path(mimix_dir, "beta.tsv")
  if (file.exists(beta_file)) file.remove(beta_file)

  julia.time <- system.time(mimix_status <- system(julia_command))
  if (mimix_status != 0 || !file.exists(beta_file)) {
    stop("MiMIX/Julia failed (exit status ", mimix_status,
         ") for ", key, "; no beta.tsv written to ", mimix_dir)
  }

  mimix_beta <- readr::read_tsv(beta_file)

  mimix_out <- as.data.frame(cbind(t(apply(mimix_beta, 2, FUN = function(x){quantile(x, c(0.025, 0.975))}))))
  ## Posterior mean alongside the interval. The interval alone says "significant"
  ## but not by how much, and the feature-level record below reports each method's
  ## own effect size.
  mimix_out$mean <- as.numeric(colMeans(mimix_beta))
  Sig <- rep("FALSE", D)
  Sig[(mimix_out$`2.5%` * mimix_out$`97.5%`) > 0] <- "TRUE"

  mimix_out <- cbind(mimix_out, Sig = Sig)
  rownames(mimix_out) <- rownames(Y)
  tax_mimix <- mimix_out[mimix_out$Sig == "TRUE",]
  tax_mimix$category <- rownames(tax_mimix)


  # Add ANCOMBC
  ancombc.time <- system.time(ancbcfit <- ancombc(data = red_phylo, assay_name = "counts",
                                                  tax_level = NULL, #phyloseq = NULL,
                                                  formula = "diet",
                                                  ## prv_cut = 0 and lib_cut = 0: no method may filter on
                                                  ## top of the paper's rule, features OR samples. Both are
                                                  ## ancombc()'s own defaults; the 1000 here was ours, and
                                                  ## at 8 samples per contrast a silently dropped sample
                                                  ## costs 12.5% of the design.
                                                  p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                                                  group = "diet", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                  max_iter = 100, conserve = TRUE, alpha = 0.05, #global = TRUE,
                                                  n_cl = 1, verbose = TRUE))

  ancbc_out <- cbind(ancbcfit$res$lfc[,c("taxon","dietVS")], W = ancbcfit$res$W$dietVS, sig = ancbcfit$res$diff_abn[,"dietVS"])

  ## ancombc carries its own `taxon` column; assert it really is in input order
  ## before using rownames(Y) as the label, instead of assuming it (it holds here,
  ## but a silent reorder would mislabel every hit the way baySeq's did).
  stopifnot(identical(as.character(ancbc_out$taxon), rownames(Y)))
  rownames(ancbc_out) <- rownames(Y)
  tax_ancbc <- ancbc_out[ancbc_out$sig == "TRUE",]
  tax_ancbc$category <- rownames(tax_ancbc)

  ancom.time <- system.time(ancfit <- ancom(data = red_phylo, assay_name = "counts",
                                            tax_level = NULL, #phyloseq = NULL,
                                            ## prv_cut / lib_cut both 0 -- see the ANCOM-BC call above.
                                            p_adj_method = "BH", prv_cut = 0, lib_cut = 0, #prv_cut : omits taxa with prevalence less than given value
                                            main_var = "diet", #adj_formula = "age + region",
                                            rand_formula = NULL, lme_control = NULL, struc_zero = FALSE,#TRUE,
                                            neg_lb = TRUE, alpha = 0.05, n_cl = 2))

  anc_out <- data.frame(W = ancfit$res$W, sig = ancfit$res$detected_0.7)

  ## Label from ancom's own taxon column rather than assuming its rows come back in
  ## input order (this is what the Prasse script already does).
  rownames(anc_out) <- as.character(ancfit$res$taxon)
  tax_anc <- as.data.frame(anc_out[anc_out$sig == "TRUE",])
  tax_anc$category <- rownames(tax_anc)

  # Add ANCOMBC2
  ancombc2.time <- system.time(ancbc2fit <- ancombc2(data = red_phylo, assay_name = "counts",
                                                     fix_formula = "diet", rand_formula = NULL,
                                                     p_adj_method = "BH", pseudo_sens = TRUE,
                                                     ## prv_cut = 0, NOT the 0.10 default: this study's
                                                     ## feature set is fixed by the PAPER's rule
                                                     ## (MIN_COUNT reads in >= MIN_SAMPLES samples), and no
                                                     ## method may apply a second filter on top of it, or
                                                     ## the twelve stop testing the same features. The
                                                     ## other five studies filter at 10% upstream, which
                                                     ## makes this argument a no-op there; here it has to
                                                     ## be switched off explicitly instead.
                                                     prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                                     group = "diet", struc_zero = FALSE, neg_lb = FALSE,
                                                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                     global = FALSE, pairwise = TRUE,
                                                     dunnet = FALSE, trend = FALSE,
                                                     iter_control = list(tol = 1e-5, max_iter = 20,
                                                                         verbose = FALSE),
                                                     em_control = list(tol = 1e-5, max_iter = 100),
                                                     lme_control = NULL,
                                                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                                     trend_control = NULL))

  ## Keep the log fold change too. Reducing ancombc2 to a detection flag threw
  ## away the only quantity that says which diet group a taxon is enriched in.
  ancbc2_lfc_col <- grep("^lfc_diet", colnames(ancbc2fit$res), value = TRUE)[1]
  ancbc2_out <- data.frame(
    taxon     = as.character(ancbc2fit$res$taxon),
    lfc_dietVS = if (is.na(ancbc2_lfc_col)) NA_real_ else as.numeric(ancbc2fit$res[[ancbc2_lfc_col]]),
    sig       = ancbc2fit$res$diff_dietVS,
    stringsAsFactors = FALSE)

  rownames(ancbc2_out) <- ancbc2fit$res$taxon
  tax_ancbc2 <- ancbc2_out[ancbc2_out$sig == "TRUE",]
  tax_ancbc2$category <- rownames(tax_ancbc2)


  ## LinDA
  meta <- as.matrix(red_phylo@sam_data$diet)
  colnames(meta) <- "conds"
  rownames(meta) <- rownames((red_phylo@otu_table))
  ## prev.cut = 0 and lib.cut = 0 for the same reason as ANCOM-BC2 above: the
  ## paper's rule defines this study's feature set and LinDA must not filter on
  ## top of it, features or samples. prev.cut = 0 is linda()'s own default; its
  ## lib.cut default is 1, and 0 is the value that keeps every sample --
  ## keep.sam is `colSums(otu.tab) >= lib.cut`. Both 0.1 and 1000 were ours.
  linda.time <- system.time(linda.obj <- linda(t(red_phylo@otu_table), meta, formula = '~conds', alpha = 0.05,
                                               prev.cut = 0, lib.cut = 0, winsor.quan = 0.97))

  linda_out <- linda.obj$output$conds
  tax_linda <- linda_out[linda_out$reject == "TRUE",]
  tax_linda$category <- rownames(tax_linda)


  ## sparse SSRV

  ## sparse SSRV takes the RAW counts (Y_raw), NOT the pseudocounted Y that the

  ## other methods use. Y is Y_raw + 1, so passing it here would make the Dirichlet

  ## concentration Y_raw + 1.5 instead of Y_raw + 0.5. Corrected 2026-08-19.

  Y_0 <- Y_raw
  X_0 <- as.numeric(X == "VD")

  S <- 50000
  D <- nrow(Y_0)
  N <- ncol(Y_0)

  resam.array <- array(NA, dim = c(S,D,N))
  sparse_time1 <- system.time(
    for(n in 1:N){
      ## Dirichlet concentration is Y_raw + 0.5, matching the other five real-data
      ## scripts. This read +1 until 2026-08-19: Y_0 is already the unpadded count
      ## matrix, so the +1 was an extra pseudocount applied only to this dataset.
      dirichlet_samples <- rdirichlet(S, as.vector(Y_0[,n]+0.5))
      resam.array[, , n] <- (log2(dirichlet_samples))
    }
  )

  sparse_time <- system.time(
    sparse_out <- sparse_auto(dat.array = resam.array, Y_0, X_0, D, N, S,
                              filtered.ind = rep(FALSE, D))
  )

  ## `pair = key` keeps one row per (level, location) in scale_model_record.rds
  ## rather than the runs overwriting each other's row -- which arm the
  ## scale-model selection picked is itself a per-run result. `key` is the bare
  ## location for cOTU (unchanged) and "<level>_<location>" for ASV and OTU.
  record_scale_model(sparse_out, dataset = "Jin2022",
                     dir = DIR_JIN, pair = key,
                     S = S, N = N)

  ## Estimated theta^perp, scored where a total was measured -------------------
  ## theta^perp is the shift that turns this method's relative log fold changes
  ## into absolute ones; under the model it IS the log2 fold change of the total
  ## load between VD and VS. BarBIQ measured that total per sample, so at the
  ## cOTU level the estimate can be scored -- bias, and whether the credible
  ## interval covers the measured value -- and not just believed.
  ##
  ## `scale` is NULL for ASV and OTU: Supplementary Data 5 and 7 are count tables
  ## with no total-abundance sheet. Those runs record the estimate with obs = NA,
  ## which reads as "not measured here" rather than as a failure to cover; each
  ## location is its own row, so the cOTU runs give two independent checks.
  ##
  ## X_0 is 1 for VD, so this is treatment minus control -- the same direction as
  ## the "theta^perp from the scale" line printed during data prep above.
  jin_perp_obs <- if (has_scale) {
    stopifnot(length(scale) == N)
    mean(log2(scale[X_0 == 1])) - mean(log2(scale[X_0 == 0]))
  } else NA_real_

  ## ...and the same total projected onto the FITTED features. The BarBIQ figure
  ## is copies/mg for the whole community; theta^perp is the log2 fold change of
  ## the total load of the kept cOTUs, because sparse_auto() is handed
  ## proportions renormalised within that set. With r_n the retained fraction,
  ## theta^perp_S = Delta log2 T_n + Delta log2 r_n, exactly.
  ##
  ## Y_all is `phylo`: the sample-subset table BEFORE both the feature filter and
  ## the optional tax_glom(). Passing `Y_kept = Y_0` rather than an id list is
  ## what makes this correct under AGGLOMERATE_TO_GENUS -- glommed ids resolve
  ## against `phylo` but select single archetype cOTUs, which would compute the
  ## retained fraction over the wrong thing and say nothing about it.
  jin_perp_matched <- if (has_scale)
    perp_matched_obs(total = scale, X = X_0,
                     Y_all = as(otu_table(phylo), "matrix"),
                     taxa_are_rows = taxa_are_rows(phylo),
                     Y_kept = Y_0, kept_taxa_are_rows = TRUE) else NULL

  perp_report <- record_theta_perp(sparse_out, dataset = "Jin2022",
                                   dir = DIR_JIN, pair = key,
                                   S = S, N = N, obs = jin_perp_obs,
                                   matched = jin_perp_matched,
                                   notes = if (has_scale) "BarBIQ absolute total abundance"
                                           else NA_character_)

  ## All THREE scale-model results from this one fit ---------------------------
  ## sparse_out is the hybrid rule's answer. sparse_auto() now also returns the
  ## laplace-only and bootstrap-only fits, so "which features would the
  ## UNSELECTED scale model have called here?" is answerable per (level,
  ## location) -- no refitting, and both arms come off the same Dirichlet draws.
  ## They go into the feature record below as two more incidence columns.
  arm_report <- report_arm_features(sparse_out, "Jin2022", pair = key)
  tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")
  tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")

  sparse_sens <- check_sparse_bw_sens(dat.array = resam.array, Y_0, X_0, D, N, S, bandwidth = c("nrd","SJ-ste"),
                                      filtered.ind = rep(FALSE, D))
  saveRDS(sparse_sens, file.path(DIR_JIN, paste0("sparse_sens", tag, ".RData")))

  tax_ssrv_qval <- sparse_out[sparse_out$sig.q == "TRUE",]
  tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
  colnames(tax_ssrv_qval)[4] <- "Sig"

  if (has_scale) {
    cat(sprintf("  sparse_SSRV (%s): %d called, %d of them in the gold set, %d not\n",
                key, nrow(tax_ssrv_qval),
                sum(tax_ssrv_qval$category %in% tax_scale$category),
                sum(!tax_ssrv_qval$category %in% tax_scale$category)))
  } else {
    cat(sprintf("  sparse_SSRV (%s): %d called\n", key, nrow(tax_ssrv_qval)))
  }


  ## Results Summary ------------------------------------------------------------

  list_tax <- list(ALDEx2_original = tax_clr, ALDEx2_def = tax_default,
                   sparse_SSRV = tax_ssrv_qval,
                   DESeq2 = tax_deseq, edgeR = tax_edgeR, limma = tax_limma, baySeq = tax_baySeq, MiMIX = tax_mimix,
                   ANCOMBC = tax_ancbc, ANCOM = tax_anc, ANCOMBC2 = tax_ancbc2,
                   LinDA = tax_linda)

  saveRDS(list_tax, file.path(DIR_JIN, paste0("relative_results", tag, ".RData")))

  ## Feature-level record -------------------------------------------------------
  ## The interpretable half of the results: not how MANY cOTUs each method called
  ## but WHICH ones, with the lineage attached, so the Jin diet contrast
  ## (VD vs VS, cecal cell samples, THIS location) can be read biologically.
  ## Everything below is bookkeeping over objects that already exist -- no model
  ## is refitted -- so it is safe to re-run on its own after loading the saved
  ## results.

  ## The classifier writes its bootstrap confidence into the name
  ## ("Vampirovibrio(0.71)") and marks unassigned ranks "*(-)". Split those apart
  ## once, here, or every downstream grouping is done on strings that never match.
  ## taxon_label / taxon_rank give the finest rank actually assigned, which is what
  ## a figure caption needs when Genus is unassigned for a large share of cOTUs.
  ## TAX_ALL, not `taxa`: at ASV and OTU level the rows are keyed by ASV/OTU id
  ## and carry the lineage inherited from the matched cOTU, so a lookup against
  ## the cOTU-keyed table would miss every id.
  taxa_clean <- clean_taxonomy(TAX_ALL)
  rank_cols  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  taxa_clean <- cbind(taxa_clean, lowest_assigned(taxa_clean, rank_cols))

  ## Reference direction, computed ONCE from the counts and independently of any
  ## method, so a single convention holds across the whole table:
  ##   POSITIVE = higher in VD (group labels are taken verbatim from the source
  ##   data; VD and VS are the two mouse diet cohorts, CE having been dropped).
  ##   lfc_prop : log2 fold change in RELATIVE abundance -- what a compositional
  ##              method sees.
  ##   lfc_abs  : the same with the measured total load (copies/mg) added back,
  ##              i.e. on the ABSOLUTE scale the gold standard targets. It is
  ##              lfc_prop shifted by THIS LOCATION's group difference in mean
  ##              log2 total load, which is not the same number at the two
  ##              locations -- see the table printed above the loop.
  grpVS <- which(X == "VS")
  grpVD <- which(X == "VD")
  stopifnot(length(grpVS) > 0, length(grpVD) > 0,
            length(grpVS) + length(grpVD) == ncol(Y))

  Ycnt <- as.matrix(Y_raw@.Data)          # pseudocounted, the matrix the methods saw
  Yobs <- as.matrix(Y_raw@.Data)      # raw counts, for prevalence / mean abundance
  log_prop <- log2(sweep(Ycnt, 2, colSums(Ycnt), "/"))
  lfc_prop <- rowMeans(log_prop[, grpVD, drop = FALSE]) -
              rowMeans(log_prop[, grpVS, drop = FALSE])
  ## lfc_abs needs a measured total load, so it exists at cOTU level only. At
  ## ASV/OTU level it is NA rather than silently equal to lfc_prop -- those two
  ## are the same number only when the load does not differ between groups, and
  ## here it differs by 0.86 log2 proximally. `direction` falls back to the
  ## compositional fold change, which is all a gene-abundance table can support.
  if (has_scale) {
    scale_log2 <- log2(scale)         # `scale` is copies/mg on the linear scale
    load_shift <- mean(scale_log2[grpVD]) - mean(scale_log2[grpVS])
    lfc_abs    <- lfc_prop + load_shift
  } else {
    load_shift <- NA_real_
    lfc_abs    <- rep(NA_real_, length(lfc_prop))
  }
  dir_ref <- if (has_scale) lfc_abs else lfc_prop

  tax_lookup <- data.frame(
    lineage_of(rownames(Ycnt), taxa_clean),
    lfc_prop   = as.numeric(lfc_prop),
    lfc_abs    = as.numeric(lfc_abs),
    direction  = ifelse(dir_ref > 0, "up_in_VD", "up_in_VS"),
    mean_count = as.numeric(rowMeans(Yobs)),
    prev_VD    = as.numeric(rowMeans(Yobs[, grpVD, drop = FALSE] > 0)),
    prev_VS    = as.numeric(rowMeans(Yobs[, grpVS, drop = FALSE] > 0)),
    level      = level,
    location   = loc,
    stringsAsFactors = FALSE, check.names = FALSE)

  ## The crosswalk back to the cell-level unit. This is the column that makes an
  ## ASV or OTU result readable against the cOTU result -- the paper's example
  ## is COTU-CM-2002 arriving as ASV-CM-28 + ASV-CM-31.
  if (!is.null(XREF))
    tax_lookup$cOTU <- XREF$cOTU[match(tax_lookup$category, XREF$feature)]

  ## Gold reference model. Unlike Prasse there is no gamma grid here: the scale
  ## model is built from the MEASURED per-sample total abundance and its replicate
  ## standard deviation, so there is exactly one reference set -- per location.
  ## Empty at ASV/OTU level, where no scale was measured.
  gold_list <- if (has_scale) list(gold_measured = tax_scale) else list()

  ## The two scale-model arms, kept OUT of list_tax on purpose: list_tax defines
  ## the method rows of res_tpfp, and the arms are scored separately there
  ## (arm_tpfp() in `real data results.R`). Here they are two more incidence
  ## columns, so the record says which features each scale model called.
  arm_list    <- list(sparse_SSRV_L = tax_ssrv_L, sparse_SSRV_B = tax_ssrv_B)

  sig_long    <- sig_features_long(list_tax,  taxa_clean, EFFECT_COL)
  arm_long    <- sig_features_long(arm_list,  taxa_clean, EFFECT_COL)
  gold_long   <- if (length(gold_list)) sig_features_long(gold_list, taxa_clean) else NULL
  record_long <- rbind(sig_long, arm_long, gold_long)   # gold: fallback finds `lfc`
  stopifnot(!is.null(record_long), nrow(record_long) > 0)
  record_long$level    <- level
  record_long$location <- loc

  ## Union of everything anyone called, one column per method then the gold column.
  ## n_methods counts the METHODS only, so it stays a measure of cross-method
  ## agreement and is not inflated by the reference.
  method_names <- names(list_tax)
  ## The arms get their own columns but stay out of method_names: n_methods
  ## measures cross-method agreement, and the arms are two views of ONE method.
  arm_names  <- names(arm_list)
  sig_matrix <- sig_features_matrix(record_long, tax_lookup,
                                    method_levels = c(method_names, arm_names, names(gold_list)))
  sig_matrix$n_methods <- rowSums(sig_matrix[, method_names, drop = FALSE])
  sig_matrix <- sig_matrix[order(-sig_matrix$n_methods, sig_matrix$category), ]
  rownames(sig_matrix) <- NULL
  ## `location` rides along from tax_lookup, so sig_matrix is self-describing
  ## once it is written to CSV and read back next to the other location's file.

  feature_record <- list(tax_lookup  = tax_lookup,    # every cOTU tested
                         sig_long    = sig_long,      # method x significant cOTU
                         gold_long   = gold_long,     # gold   x significant cOTU
                         record_long = record_long,   # the two stacked
                         sig_matrix  = sig_matrix,    # union, wide incidence
                         location    = loc,
                         load_shift  = load_shift,    # log2(total load) VD - VS
                         samples     = rownames(metadata))

  saveRDS(feature_record, file.path(DIR_JIN, paste0("jin_feature_record", tag, ".RData")))
  ## CSVs alongside the RDS: these are the tables that get read by eye when writing
  ## up the biology, and .RData cannot be opened outside R.
  write.csv(tax_lookup,  file.path(DIR_JIN, paste0("jin_tax_lookup",  tag, ".csv")), row.names = FALSE)
  write.csv(record_long, file.path(DIR_JIN, paste0("jin_sig_long",    tag, ".csv")), row.names = FALSE)
  write.csv(sig_matrix,  file.path(DIR_JIN, paste0("jin_sig_matrix",  tag, ".csv")), row.names = FALSE)

  all_cols <- c(method_names, arm_names, names(gold_list))

  cat("\n== significant ", level, "s per method (Jin: VD vs VS, cell samples, ", loc, ") ==\n", sep = "")
  print(data.frame(
    Method = all_cols,
    n      = as.integer(colSums(sig_matrix[, all_cols, drop = FALSE])),
    row.names = NULL))

  ## WHICH features, not just how many -- the whole point of the record. The
  ## effect-size column shown is lfc_abs where a scale was measured and
  ## lfc_prop otherwise, since lfc_abs is NA at ASV/OTU level.
  show_cols <- c("category", if (!is.null(XREF)) "cOTU", "Phylum", "Family",
                 "taxon_label", "taxon_rank",
                 if (has_scale) "lfc_abs" else "lfc_prop", "direction", "n_methods")
  cat("\n== ", level, "s called by the most methods (", loc, ") ==\n", sep = "")
  print(head(sig_matrix[, show_cols], 20), row.names = FALSE)

  ## Head-to-head against the reference, spelled out taxon by taxon: agreement,
  ## method-only (false positive against this reference) and gold-only (miss).
  ## cOTU only -- there is no reference at the other two levels.
  res_tpfp <- NULL
  ## NULL at the ASV and OTU levels, where there is no measured reference to
  ## score against and CKF is not run. Carried on the result object either way so
  ## downstream scripts can tell "not run here" from "not recorded".
  ckf <- NULL
  if (has_scale) {
    cmp <- sig_matrix[sig_matrix$sparse_SSRV == 1 | sig_matrix$gold_measured == 1,
                      c("category", "Phylum", "taxon_label", "lfc_abs", "direction",
                        "sparse_SSRV", "gold_measured")]
    cmp$agreement <- ifelse(cmp$sparse_SSRV == 1 & cmp$gold_measured == 1, "both",
                     ifelse(cmp$sparse_SSRV == 1, "sparse_SSRV only", "gold only"))
    cat("\n== sparse_SSRV vs the gold model (", loc, ") ==\n", sep = "")
    print(cmp[order(cmp$agreement, -cmp$lfc_abs), ], row.names = FALSE)

    res_tpfp <- summarize_tpfp(tax_gold = tax_scale, list_tax)
    print(res_tpfp)

    ## ---- CKF -------------------------------------------------------------------
    ## Fitted on Y_0 / X_0, the SAME matrix and contrast handed to sparse_auto() a
    ## few lines above, so CKF is scored on exactly the taxa the other methods were
    ## scored on. The standalone script this replaces re-parsed the study and
    ## rebuilt the filter itself, which is how its taxon set drifted from this one.
    ## Inside `if (has_scale)`, so CKF runs for the cOTU runs that HAVE a measured
    ## reference (pooled, prox and dist) and not at the ASV/OTU levels, where there is no
    ## gold set to score against and res_tpfp is not built either. `pair = key`
    ## matches the keys scale_model_record.rds and theta_perp_record.rds already
    ## use for this study, so all three records line up run for run.
    ##
    ## The legacy jin_mrckf_out.RData is written by the POOLED run only. Every
    ## subset would otherwise write the same filename and whichever ran last would
    ## win, silently putting one subset's numbers behind the row the summary script
    ## labels another. Pooled is the subset `real data results.R` reports (JIN_TAG
    ## = ""), so it is the one that owns the legacy name. The keyed ckf_record.rds
    ## carries every subset under its own name.
    ##
    ## Note `real data results.R` does NOT read this file: res_tpfp already has the
    ## CKF row merged by add_ckf_row() below. It is kept for older readers.
    ckf <- run_ckf(Y_0, X_0, truth = tax_scale$category,
                   R = CKF_R, q = CKF_Q, offset = CKF_OFFSET,
                   se_target = CKF_SE_TARGET, R_min = CKF_R_MIN,
                   derand_thresh = CKF_DERAND, seed = 1234,
                   dataset = "Jin2022", pair = key)
    res_tpfp <- add_ckf_row(res_tpfp, ckf)
    record_ckf(ckf, dir = DIR_JIN,
               file = if (identical(loc, "both")) "jin_mrckf_out.RData" else NULL)
    ## What R would have been enough, from the fits just done.
    print(ckf_stability(ckf, se_target = 0.1))

    saveRDS(res_tpfp, file.path(DIR_JIN, paste0("res_tpfp", tag, ".RData")))
  }

  ## Everything this run produced, kept in memory for the summary sections.
  res_by_loc[[loc]] <- list(tag            = tag,
                            level          = level,
                            key            = key,
                            samples        = rownames(metadata),
                            n              = nrow(metadata),
                            D              = D,
                            load_shift     = load_shift,
                            tax_scale      = tax_scale,
                            list_tax       = list_tax,
                            tax_lookup     = tax_lookup,
                            sig_matrix     = sig_matrix,
                            record_long    = record_long,
                            feature_record = feature_record,
                            res_tpfp       = res_tpfp,
                            ## The CKF run for this subset. Its `reference_taxon`
                            ## is the taxon MRCKF drops as the log-ratio
                            ## reference and can therefore never select -- the
                            ## one fact about a CKF fit that is not recoverable
                            ## from its output afterwards, so it travels here as
                            ## well as in ckf_record.rds.
                            ckf            = ckf,
                            sparse_out     = sparse_out,
                            sparse_sens    = sparse_sens)
}

################################################################################
## CROSS-LOCATION SUMMARY, for this feature level
##
## The runs above are complete and self-contained. This section does no fitting
## -- it only lines them up so the biological question the split was made for
## ("is the diet effect regional?") can be answered, and so that a feature
## missing from one table can be told apart from one that was never tested there
## (the abundance filter is per location).
##
## Reminder when reading it: the same eight mice appear in both runs, so
## "prox only" means the signal is regional WITHIN those mice, not that it
## replicated or failed to replicate in a second cohort.
################################################################################

## The pooled ("both") run is EXCLUDED from everything below. It is the primary
## analysis, not a third location: it is the same eight mice a second time, over
## a feature set that is a superset of both stratified ones. Letting it in would
## corrupt this section in two specific ways --
##   * the "<loc> only" columns are setdiff(loc, union(other locs)), so a cOTU
##     called at prox AND in the pooled run would count as neither "prox only"
##     nor "dist only", and simply vanish from the split table;
##   * `tested_in` and `common_tested` would intersect against the pooled feature
##     set, which is larger than either location's, silently shrinking the set of
##     cOTUs the two locations are compared on.
## Both would be read as biology. Compare the stratified runs to each other here
## and leave the pooled run to `real data results.R`.
res_strat <- res_by_loc[setdiff(names(res_by_loc), "both")]

if (length(res_strat) > 1) {

  locs <- names(res_strat)

  ## 1. How many features each method called, side by side.
  method_names <- names(res_strat[[1]]$list_tax)
  ## The two scale-model arms are columns of sig_matrix but not of list_tax, so
  ## they are read off the matrix itself: present for any run made since
  ## sparse_auto() started returning all three arms, silently absent for an older
  ## one, which keeps this section runnable on results already on disk.
  arm_cols     <- intersect(c("sparse_SSRV_L", "sparse_SSRV_B"),
                            colnames(res_strat[[1]]$sig_matrix))
  all_cols     <- c(method_names, arm_cols, if (has_scale) "gold_measured")
  ## What n_calls counts: methods plus the reference, never the arms. The arms
  ## are two more views of sparse_SSRV, so counting them would treble its weight
  ## in a column that is supposed to measure cross-METHOD agreement.
  count_cols   <- c(method_names, if (has_scale) "gold_measured")
  n_by_loc <- sapply(res_strat, function(r)
    as.integer(colSums(r$sig_matrix[, all_cols, drop = FALSE])))
  rownames(n_by_loc) <- all_cols
  cat("\n== significant ", level, "s per method, by location ==\n", sep = "")
  print(as.data.frame(n_by_loc))

  ## 2. Design facts worth quoting in the text.
  cat("\n== per-location design (", level, ") ==\n", sep = "")
  design_tbl <- data.frame(
    level      = level,
    location   = locs,
    n_samples  = sapply(res_strat, `[[`, "n"),
    D_tested   = sapply(res_strat, `[[`, "D"),
    row.names  = NULL)
  if (has_scale) {
    design_tbl$load_shift <- round(sapply(res_strat, `[[`, "load_shift"), 3)
    design_tbl$n_gold     <- sapply(res_strat, function(r) nrow(r$tax_scale))
  }
  print(design_tbl)

  ## 3. The gold sets themselves. This is the reference against which every
  ## method is scored, and it is refitted per location, so whether the two
  ## references even agree is the first thing to check. cOTU only.
  if (has_scale) {
    gold_sets <- lapply(res_strat, function(r) r$tax_scale$category)
    cat("\n== gold (measured-scale) reference sets ==\n")
    cat("  shared by all locations : ", length(Reduce(intersect, gold_sets)), "\n", sep = "")
    for (l in locs) {
      cat("  ", l, " only : ", length(setdiff(gold_sets[[l]], unlist(gold_sets[setdiff(locs, l)]))),
          "\n", sep = "")
    }
  }

  ## 4. One long table: every (location, cOTU) that anyone called, plus whether
  ## that cOTU was even TESTED at the other location. `tested_at` is what stops a
  ## filter-driven absence from being read as a biological one.
  tested_at <- lapply(res_strat, function(r) r$tax_lookup$category)
  cross_long <- do.call(rbind, lapply(locs, function(l) {
    sm <- res_strat[[l]]$sig_matrix
    keep_cols <- intersect(c("category", "cOTU", "Phylum", "Class", "Order", "Family", "Genus",
                             "taxon_label", "taxon_rank", "lfc_prop", "lfc_abs",
                             "direction", "mean_count", "n_methods",
                             all_cols), colnames(sm))
    d <- sm[, keep_cols, drop = FALSE]
    d$level     <- level
    d$location  <- l
    d$tested_in <- sapply(d$category, function(cat_)
      paste(locs[sapply(tested_at[locs], function(tt) cat_ %in% tt)], collapse = "+"))
    d
  }))
  cross_long <- cross_long[order(cross_long$category, cross_long$location), ]
  rownames(cross_long) <- NULL
  cross_long_file <- file.path(DIR_JIN, paste0("jin_sig_by_location",
                                               if (has_scale) "" else paste0("_", level),
                                               "_", DATESTAMP, ".csv"))
  write.csv(cross_long, cross_long_file, row.names = FALSE)

  ## 5. Per method, the prox/dist split of its own calls. `all_locations` is the
  ## features a method called at BOTH locations -- the regionally consistent part
  ## of its answer -- and the "only" columns are the regional part. Restricted to
  ## features tested at both locations so the filter cannot inflate "only".
  ## This is the table that answers the paper's question at each level: does the
  ## proximal/distal asymmetry it reports for cOTUs survive at ASV and OTU level?
  common_tested <- Reduce(intersect, tested_at)
  split_tbl <- do.call(rbind, lapply(all_cols, function(mth) {
    sets <- lapply(res_strat, function(r) {
      sm <- r$sig_matrix
      intersect(sm$category[sm[[mth]] == 1], common_tested)
    })
    ## Written for any number of locations rather than hardcoding prox/dist, so
    ## adding "both" to LOCATIONS for a side-by-side does not break this table.
    row <- data.frame(method = mth,
                      all_locations = length(Reduce(intersect, sets)),
                      stringsAsFactors = FALSE)
    for (l in locs)
      row[[paste0(l, "_only")]] <-
        length(setdiff(sets[[l]], unlist(sets[setdiff(locs, l)])))
    row
  }))
  cat("\n== each method's calls split by location (", level, "s tested at both) ==\n", sep = "")
  print(split_tbl, row.names = FALSE)
  split_tbl$level <- level

  ## 6. Every feature anyone called, with each method's verdict at each location
  ## side by side. This is the table the user asked for and the one to read when
  ## writing up: one row per feature, a 0/1 per (method, location), and the
  ## effect size at each location so a call can be weighed against its size.
  ## Built at every level, scale or no scale.
  called <- unique(unlist(lapply(res_strat, function(r)
    r$sig_matrix$category[rowSums(r$sig_matrix[, all_cols, drop = FALSE]) > 0])))

  if (length(called)) {
    lk_all <- do.call(rbind, lapply(res_strat, function(r) r$tax_lookup))
    lk_all <- lk_all[!duplicated(lk_all$category), ]
    j <- match(called, lk_all$category)
    calls <- data.frame(level       = level,
                        category    = called,
                        stringsAsFactors = FALSE)
    if (!is.null(XREF)) calls$cOTU <- lk_all$cOTU[j]
    calls$Phylum      <- lk_all$Phylum[j]
    calls$Family      <- lk_all$Family[j]
    calls$taxon_label <- lk_all$taxon_label[j]

    ## Effect size per location, then the per-method 0/1 flags. `NA` in an
    ## effect-size column means the feature was not tested at that location.
    eff <- if (has_scale) "lfc_abs" else "lfc_prop"
    for (l in locs) {
      r <- res_strat[[l]]
      calls[[paste0(eff, "_", l)]] <- round(r$tax_lookup[[eff]][match(called, r$tax_lookup$category)], 3)
    }
    for (l in locs) {
      sm <- res_strat[[l]]$sig_matrix
      k  <- match(called, sm$category)
      for (mth in all_cols)
        calls[[paste0(mth, "_", l)]] <- ifelse(is.na(k), NA_integer_, sm[[mth]][k])
    }
    ## Named explicitly rather than matched by a "^(method|...)_" regex: with the
    ## arms in play that pattern also matches sparse_SSRV_L_prox through its
    ## sparse_SSRV prefix, which would put the arms back into the count.
    calls$n_calls <- rowSums(calls[, as.vector(outer(count_cols, locs, paste, sep = "_")),
                                   drop = FALSE], na.rm = TRUE)
    calls <- calls[order(-calls$n_calls, calls$category), ]
    rownames(calls) <- NULL

    calls_file <- file.path(DIR_JIN, paste0("jin_method_calls_by_location_", level, "_", DATESTAMP, ".csv"))
    write.csv(calls, calls_file, row.names = FALSE)
    cat("\n== features called by any method, with each method's verdict per location (", level, ") ==\n", sep = "")
    print(head(calls[, c("category", if (!is.null(XREF)) "cOTU", "taxon_label",
                         paste0(eff, "_", locs), "n_calls")], 30), row.names = FALSE)
  } else {
    calls <- NULL
    calls_file <- NA_character_
    cat("\n== no ", level, " was called significant by any method at any location ==\n", sep = "")
  }

  ## 7. The gold-vs-sparse_SSRV biology table -- cOTU only, since it is built on
  ## the reference set and on lfc_abs.
  if (has_scale) {

    ## 6. The biology table: the gold reference and sparse_SSRV, taxon by taxon,
    ## with the effect size at each location side by side. This is the table to
    ## read when writing the discussion -- a taxon with the same sign and similar
    ## lfc_abs at both locations is a diet effect; one that flips or collapses is a
    ## regional effect.
    focus <- unique(unlist(lapply(res_strat, function(r)
      c(r$tax_scale$category, r$sig_matrix$category[r$sig_matrix$sparse_SSRV == 1]))))
    lk <- res_strat[[1]]$tax_lookup
    bio <- data.frame(category    = focus,
                      Phylum      = lk$Phylum[match(focus, lk$category)],
                      Family      = lk$Family[match(focus, lk$category)],
                      taxon_label = lk$taxon_label[match(focus, lk$category)],
                      stringsAsFactors = FALSE)
    ## Fill lineage for cOTUs absent from the first location's lookup.
    for (l in locs[-1]) {
      lk2 <- res_strat[[l]]$tax_lookup
      miss <- is.na(bio$taxon_label)
      if (any(miss)) {
        j <- match(bio$category[miss], lk2$category)
        bio$Phylum[miss]      <- lk2$Phylum[j]
        bio$Family[miss]      <- lk2$Family[j]
        bio$taxon_label[miss] <- lk2$taxon_label[j]
      }
    }
    for (l in locs) {
      r  <- res_strat[[l]]
      j  <- match(bio$category, r$tax_lookup$category)
      sm <- r$sig_matrix
      k  <- match(bio$category, sm$category)
      bio[[paste0("lfc_abs_", l)]] <- round(r$tax_lookup$lfc_abs[j], 3)
      bio[[paste0("gold_",    l)]] <- ifelse(is.na(k), 0L, sm$gold_measured[k])
      bio[[paste0("ssrv_",    l)]] <- ifelse(is.na(k), 0L, sm$sparse_SSRV[k])
    }
    bio <- bio[order(-rowSums(bio[, grep("^gold_|^ssrv_", colnames(bio))], na.rm = TRUE),
                     bio$category), ]
    rownames(bio) <- NULL
    write.csv(bio, file.path(DIR_JIN, paste0("jin_gold_vs_ssrv_by_location_", DATESTAMP, ".csv")),
              row.names = FALSE)
    cat("\n== gold / sparse_SSRV calls and absolute effect sizes, by location ==\n")
    print(head(bio, 40), row.names = FALSE)

  }

  cat("\nWrote:\n  ", cross_long_file, "\n", sep = "")
  if (!is.na(calls_file)) cat("  ", calls_file, "\n", sep = "")
  if (has_scale) cat("  ", file.path(DIR_JIN, paste0("jin_gold_vs_ssrv_by_location_", DATESTAMP, ".csv")), "\n", sep = "")

  ## Kept for the cross-level section below.
  res_all[[level]] <- res_by_loc
  attr(res_all[[level]], "split_tbl") <- split_tbl
  attr(res_all[[level]], "calls")     <- calls
}

## Everything this level produced, keyed by location, in one object -- so the
## write-up can be built without re-running any model. cOTU keeps the original
## file name so anything already pointed at it keeps working.
saveRDS(res_by_loc, file.path(DIR_JIN, paste0("jin_results_by_location",
                                              if (has_scale) "" else paste0("_", level),
                                              "_", DATESTAMP, ".RData")))

}   ## end: for (level in LEVELS)

################################################################################
## CROSS-LEVEL SUMMARY -- cell abundance vs 16S rRNA gene abundance
##
## The paper's central claim, in one table. Jin et al. report that the
## proximal/distal asymmetry they see in cell abundance (cOTUs) is NOT
## reproduced by conventional gene-abundance analyses:
##
##            proximal   distal
##   cOTU        8          1        (their DESeq2 result)
##   ASV         8          6
##   OTU         9          6
##
## Their numbers come from DESeq2 alone. What follows is the same comparison run
## across all twelve methods, which is the question this project is actually
## asking: is that asymmetry a property of the measurement (cells vs genes), or
## of the estimator? A method that shows the cOTU pattern at ASV level too is
## saying something different from one that does not.
##
## There is no reference set here and no TP/FP -- only who called what.
################################################################################

if (length(res_all) > 1) {

  levs <- names(res_all)

  ## Significant features per (method, level, location). The core table.
  n_cube <- do.call(rbind, lapply(levs, function(lv) {
    rb <- res_all[[lv]]
    ## Methods plus the two scale-model arms, so the paper's central table can be
    ## read per scale model as well as per method.
    mn <- c(names(rb[[1]]$list_tax),
            intersect(c("sparse_SSRV_L", "sparse_SSRV_B"), colnames(rb[[1]]$sig_matrix)))
    do.call(rbind, lapply(names(rb), function(l) {
      sm <- rb[[l]]$sig_matrix
      data.frame(level = lv, location = l, method = mn,
                 n = as.integer(colSums(sm[, mn, drop = FALSE])),
                 stringsAsFactors = FALSE)
    }))
  }))

  ## Wide: one row per method, one column per level x location.
  n_wide <- Reduce(function(a, b) merge(a, b, by = "method", all = TRUE),
                   lapply(levs, function(lv) {
                     d <- n_cube[n_cube$level == lv, ]
                     w <- reshape(d[, c("method", "location", "n")], idvar = "method",
                                  timevar = "location", direction = "wide")
                     colnames(w)[-1] <- paste0(lv, "_", sub("^n\\.", "", colnames(w)[-1]))
                     w
                   }))
  ## Preserve the method order the analysis uses rather than merge()'s sort.
  mn0 <- c(names(res_all[[1]][[1]]$list_tax),
           intersect(c("sparse_SSRV_L", "sparse_SSRV_B"),
                     colnames(res_all[[1]][[1]]$sig_matrix)))
  n_wide <- n_wide[match(mn0, n_wide$method), ]
  rownames(n_wide) <- NULL
  cat("\n== significant features per method, by level and location ==\n")
  cat("   (cell abundance = cOTU; 16S rRNA gene abundance = ASV, OTU)\n")
  print(n_wide, row.names = FALSE)

  write.csv(n_cube, file.path(DIR_JIN, paste0("jin_n_sig_by_level_location_", DATESTAMP, ".csv")),
            row.names = FALSE)
  write.csv(n_wide, file.path(DIR_JIN, paste0("jin_n_sig_wide_", DATESTAMP, ".csv")),
            row.names = FALSE)

  ## Feature counts tested at each level, against the paper's own numbers.
  cat("\n== features tested (after the >=", MIN_COUNT, "-in->=", MIN_SAMPLES,
      " filter), against the paper ==\n", sep = "")
  paper_D <- c(cOTU_prox = 150, cOTU_dist = 153,
               ASV_prox  = 236, ASV_dist  = 223,
               OTU_prox  = 138, OTU_dist  = 131)
  dtab <- do.call(rbind, lapply(levs, function(lv)
    do.call(rbind, lapply(names(res_all[[lv]]), function(l)
      data.frame(level = lv, location = l, D = res_all[[lv]][[l]]$D,
                 paper = unname(paper_D[paste0(lv, "_", l)]), stringsAsFactors = FALSE)))))
  dtab$match <- ifelse(is.na(dtab$paper), "", ifelse(dtab$D == dtab$paper, "yes", "NO"))
  print(dtab, row.names = FALSE)

  ## The prox/dist asymmetry per method, per level -- the direct analogue of the
  ## paper's 8-vs-1 / 8-vs-6 / 9-vs-6 comparison. A positive `asymmetry` means a
  ## method finds more proximally, which is the cOTU-level pattern the paper
  ## reports and says the gene-abundance levels lack.
  if (all(c("prox", "dist") %in% n_cube$location)) {
    asym <- do.call(rbind, lapply(levs, function(lv) {
      d <- n_cube[n_cube$level == lv, ]
      pr <- d$n[d$location == "prox"][match(mn0, d$method[d$location == "prox"])]
      di <- d$n[d$location == "dist"][match(mn0, d$method[d$location == "dist"])]
      data.frame(level = lv, method = mn0, prox = pr, dist = di,
                 asymmetry = pr - di, stringsAsFactors = FALSE)
    }))
    cat("\n== proximal minus distal, per method and level ==\n")
    print(reshape(asym[, c("method", "level", "asymmetry")], idvar = "method",
                  timevar = "level", direction = "wide"), row.names = FALSE)
  }

  ## Where the same cell-level unit was called at more than one level. The
  ## crosswalk makes this possible: an ASV or OTU hit is mapped back to its cOTU,
  ## so a cOTU that BarBIQ flags and the conventional methods miss (or the
  ## reverse) can be named rather than just counted. COTU-CM-2002 is the paper's
  ## worked example -- it arrives as ASV-CM-28 + ASV-CM-31 and as OTU-CM-029.
  to_cotu <- function(lv) {
    rb <- res_all[[lv]]
    ## Methods only, on purpose: this table asks which LEVELS called a cell-level
    ## unit, and folding the two scale-model arms in would count sparse_SSRV
    ## three times in n_methods. The arms' own calls are in sig_matrix and in
    ## sparse_arm_features.csv.
    mn <- names(rb[[1]]$list_tax)
    do.call(rbind, lapply(names(rb), function(l) {
      sm  <- rb[[l]]$sig_matrix
      hit <- sm$category[rowSums(sm[, mn, drop = FALSE]) > 0]
      if (!length(hit)) return(NULL)
      cot <- if (!is.null(sm$cOTU)) sm$cOTU[match(hit, sm$category)] else hit
      data.frame(level = lv, location = l, feature = hit, cOTU = cot,
                 n_methods = sm$n_methods[match(hit, sm$category)],
                 stringsAsFactors = FALSE)
    }))
  }
  hits <- do.call(rbind, lapply(levs, to_cotu))
  if (!is.null(hits) && nrow(hits)) {
    hits <- hits[!is.na(hits$cOTU), ]
    agree <- aggregate(level ~ cOTU + location, data = hits,
                       FUN = function(x) paste(sort(unique(x)), collapse = "+"))
    colnames(agree)[3] <- "levels_calling"
    agree <- agree[order(-nchar(agree$levels_calling), agree$cOTU), ]
    rownames(agree) <- NULL
    write.csv(hits,  file.path(DIR_JIN, paste0("jin_hits_all_levels_", DATESTAMP, ".csv")), row.names = FALSE)
    write.csv(agree, file.path(DIR_JIN, paste0("jin_cotu_agreement_across_levels_", DATESTAMP, ".csv")), row.names = FALSE)
    cat("\n== cell-level units called, and which levels called them ==\n")
    print(head(agree, 40), row.names = FALSE)
  }

  saveRDS(res_all, file.path(DIR_JIN, paste0("jin_results_all_levels_", DATESTAMP, ".RData")))
  cat("\nWrote:\n  ", file.path(DIR_JIN, paste0("jin_n_sig_by_level_location_", DATESTAMP, ".csv")),
      "\n  ", file.path(DIR_JIN, paste0("jin_results_all_levels_", DATESTAMP, ".RData")), "\n", sep = "")
}

################################################################################
## SIGNIFICANT-FEATURE SUMMARY, under the paper's own rule
##
## Nested exactly as requested:
##
##   sig_summary[[level]][[location]][[method]]  ->  one data.frame
##       level    : "cOTU" | "ASV" | "OTU"
##       location : "prox" | "dist"
##       method   : the twelve methods, plus "gold_measured" at cOTU level
##
## Each data.frame holds the features that method called significant under BOTH
## halves of the paper's rule -- the method's own significance call (FDR < 0.05
## where it reports an FDR) AND |log2 fold change| > log2(SIG_FC) -- with:
##
##   the taxonomic lineage, Kingdom through Genus  (the "tree" information; note
##     there is no phylogeny in these files, only the RDP classifier's lineage,
##     so this is taxonomy and not a tree)
##   taxon_label / taxon_rank : finest rank actually assigned
##   direction                : up_in_VD / up_in_VS
##   lfc_prop, lfc_abs        : the shared log2 fold changes, one convention
##   method_lfc               : the method's OWN reported effect size
##   log2fc_used / fold_change: what the >2 test was actually applied to
##
## A note on the FDR half. Only some of these methods report an FDR: DESeq2,
## edgeR, limma, LinDA, ALDEx2 and ANCOM-BC do, while ANCOM uses a W-statistic
## cutoff, MiMIX a credible interval excluding zero, sparse_SSRV its own q-value
## rule and baySeq an FWER. Membership of `list_tax` already IS each method's
## significance call, so that is what is used; where a padj exists it is also
## checked against SIG_FDR, so tightening SIG_FDR tightens those methods and
## leaves the others at their own rule. That is stated rather than hidden
## because "FDR < 0.05" is not a thing every one of these methods computes.
##
## No model is refitted here -- it reads res_all -- so it can be re-run on its
## own after loading jin_results_all_levels_<date>.RData.
################################################################################

summarize_significant <- function(res_all, fdr = SIG_FDR, fc = SIG_FC, basis = FC_BASIS) {

  stopifnot(basis %in% c("lfc_prop", "lfc_abs", "method"))
  lfc_cut   <- log2(fc)
  rank_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

  ## Column skeleton, so a method that called nothing still returns a 0-row
  ## data.frame with the right shape instead of NULL. Downstream code that
  ## rbinds these must not have to special-case the empty ones.
  empty_row <- function(cols) {
    d <- as.data.frame(setNames(rep(list(logical(0)), length(cols)), cols),
                       stringsAsFactors = FALSE)
    d
  }

  lapply(res_all, function(rb) {                      # level
    lapply(rb, function(r) {                          # location

      lk  <- r$tax_lookup
      rl  <- r$record_long
      lvl <- r$level
      loc <- if (!is.null(r$feature_record$location)) r$feature_record$location else NA_character_
      has_xref <- !is.null(lk$cOTU)

      cols <- c("level", "location", "method", "category",
                if (has_xref) "cOTU", rank_cols, "taxon_label", "taxon_rank",
                "direction", "lfc_prop", "lfc_abs", "method_lfc", "method_lfc_col",
                "log2fc_used", "fold_change", "fc_basis", "padj",
                "mean_count", "prev_VD", "prev_VS")

      ## The CANONICAL method set, taken from list_tax and not from whatever
      ## happens to appear in record_long. A method that called nothing must
      ## still be present as a 0-row data.frame: otherwise the list has a
      ## different shape in every cell, sig_summary[[lv]][[loc]][["MiMIX"]]
      ## returns NULL in some cells and a table in others, and the count table
      ## shows NA where the honest answer is 0.
      ## The two scale-model arms are part of that canonical set: record_long
      ## carries their rows, so leaving them out here would drop them from the
      ## fold-change / FDR summary while every other table shows them.
      methods <- c(names(r$list_tax),
                   intersect(c("sparse_SSRV_L", "sparse_SSRV_B"), unique(rl$Method)),
                   if (!is.null(r$tax_scale)) "gold_measured")

      out <- lapply(methods, function(mth) {

        rows <- rl[rl$Method == mth, , drop = FALSE]
        if (!nrow(rows)) return(empty_row(cols))

        j <- match(rows$category, lk$category)

        ## The quantity the fold-change test is applied to. Under "method" a
        ## method that reports no effect size falls back to lfc_prop rather than
        ## being silently dropped -- ANCOM and baySeq would otherwise vanish
        ## from every table.
        used <- switch(basis,
                       lfc_prop = lk$lfc_prop[j],
                       lfc_abs  = lk$lfc_abs[j],
                       method   = ifelse(is.na(rows$stat), lk$lfc_prop[j], rows$stat))

        ## SIGNIFICANCE IS padj ALONE. `padj` is NA for the methods that do not
        ## report one; those keep their own significance call, which membership
        ## of list_tax already encodes.
        ##
        ## The |log2 FC| > log2(SIG_FC) half of Jin et al.'s figure-legend rule is
        ## NOT applied here any more, for two reasons.
        ##
        ##   1. Consistency. res_tpfp -- the object the main DA results table is
        ##      built from, for THIS study and all five others -- scores each
        ##      method's padj-based calls against the gold reference and never
        ##      looks at a fold change. Applying an extra magnitude filter only in
        ##      this summary made the two disagree about what a method called.
        ##   2. The filter was dropping the strongest calls, not the weakest. It
        ##      was guarded by `!is.na(used)`, and lfc_prop is NaN exactly when
        ##      BOTH arms contain a zero -- a taxon absent from one diet arm and
        ##      missing in one mouse of the other. is.na(NaN) is TRUE in R, so
        ##      those were silently dropped, while +-Inf (one arm all-zero, the
        ##      other complete) was kept by abs(Inf) > cut. Whether a structural
        ##      -zero taxon survived came down to whether one mouse in the present
        ##      arm happened to have a zero count. In the pooled cOTU run that
        ##      removed 18 of the 40 called features, including 2 of sparse SSRV's
        ##      4 calls -- which is what made jin_sig_summary report 2 where
        ##      res_tpfp reported 4.
        ##
        ## lfc_prop / lfc_abs / fold_change are still REPORTED in the output, and
        ## the paper's rule can still be recovered from them downstream. Nothing
        ## is filtered on them here.
        keep <- is.na(rows$padj) | rows$padj < fdr

        if (!any(keep)) return(empty_row(cols))
        rows <- rows[keep, , drop = FALSE]
        j    <- j[keep]
        used <- used[keep]

        d <- data.frame(level = lvl, location = loc, method = mth,
                        category = rows$category, stringsAsFactors = FALSE)
        if (has_xref) d$cOTU <- lk$cOTU[j]
        for (rk in rank_cols) d[[rk]] <- lk[[rk]][j]
        d$taxon_label    <- lk$taxon_label[j]
        d$taxon_rank     <- lk$taxon_rank[j]
        d$direction      <- lk$direction[j]
        d$lfc_prop       <- round(lk$lfc_prop[j], 4)
        d$lfc_abs        <- round(lk$lfc_abs[j], 4)
        d$method_lfc     <- round(rows$stat, 4)
        d$method_lfc_col <- rows$stat_col
        d$log2fc_used    <- round(used, 4)
        ## Signed fold change on the linear scale, the units the paper quotes
        ## ("fold change: 320-1100"). Negative means the ratio runs the other way.
        d$fold_change    <- round(sign(used) * 2^abs(used), 3)
        d$fc_basis       <- basis
        d$padj           <- rows$padj
        d$mean_count     <- round(lk$mean_count[j], 2)
        d$prev_VD        <- lk$prev_VD[j]
        d$prev_VS        <- lk$prev_VS[j]

        d <- d[order(-abs(d$log2fc_used)), cols, drop = FALSE]
        rownames(d) <- NULL
        d
      })
      names(out) <- methods
      out
    })
  })
}

sig_summary <- summarize_significant(res_all)

saveRDS(sig_summary, file.path(DIR_JIN, paste0("jin_sig_summary_", DATESTAMP, ".RData")))

## Flat version of the same thing -- every (level, location, method, feature) on
## one row. The nested list is what was asked for and is the right shape for
## programmatic access; this is the one to open in a spreadsheet.
sig_parts <- unlist(lapply(sig_summary, function(rb) lapply(rb, function(x) x)),
                    recursive = FALSE)
sig_parts <- unlist(sig_parts, recursive = FALSE)
sig_parts <- sig_parts[sapply(sig_parts, nrow) > 0]

## `cOTU` exists at ASV and OTU level only, so the three levels do not share a
## column set. Align on the union -- rbind would otherwise fail outright on the
## mismatch, and dropping to the intersection would throw away the crosswalk
## that makes the ASV/OTU rows readable.
sig_flat <- NULL
if (length(sig_parts)) {
  all_cn <- unique(unlist(lapply(sig_parts, colnames)))
  sig_parts <- lapply(sig_parts, function(d) {
    for (cn in setdiff(all_cn, colnames(d))) d[[cn]] <- NA
    d[, all_cn, drop = FALSE]
  })
  sig_flat <- do.call(rbind, sig_parts)
}
if (!is.null(sig_flat)) {
  rownames(sig_flat) <- NULL
  write.csv(sig_flat, file.path(DIR_JIN, paste0("jin_sig_summary_flat_", DATESTAMP, ".csv")),
            row.names = FALSE)
}

## What this study CALLS significant (padj alone -- see summarize_significant),
## beside what Jin et al.'s figure-legend rule WOULD have called on the same fits.
##
## The second column is reported, never applied. It is here so the paper's own
## numbers (cOTU 8/1, ASV 8/6, OTU 9/6) stay comparable without that rule being
## imposed on this project's calls, which are padj-based for every study.
##
## Read the gap between the columns with the zero-inflation caveat in mind: a
## taxon absent from one diet arm has an infinite compositional fold change, and
## one that is ALSO missing in a single mouse of the other arm has an undefined
## (NaN) one. Neither is a small effect, but only the first passes an
## "abs(log2 FC) > cut" test -- so n_padj_and_fc understates the strong calls.
cat("\n== significant features: padj < ", SIG_FDR,
    "  (reference column: also fold change > ", SIG_FC,
    ", basis ", FC_BASIS, ") ==\n", sep = "")
lfc_cut_report <- log2(SIG_FC)
sig_counts <- do.call(rbind, lapply(names(sig_summary), function(lv)
  do.call(rbind, lapply(names(sig_summary[[lv]]), function(l) {
    parts <- sig_summary[[lv]][[l]]
    mn    <- names(parts)
    data.frame(level = lv, location = l, method = mn,
               n_padj = as.integer(sapply(parts, nrow)),
               n_padj_and_fc = as.integer(sapply(parts, function(x)
                 sum(!is.na(x$log2fc_used) & abs(x$log2fc_used) > lfc_cut_report))),
               stringsAsFactors = FALSE)
  }))))
rownames(sig_counts) <- NULL
sig_counts$cell <- paste0(sig_counts$level, "_", sig_counts$location)
sig_wide <- reshape(sig_counts[, c("method", "cell", "n_padj")],
                    idvar = "method", timevar = "cell", direction = "wide")
colnames(sig_wide) <- sub("^n_padj\\.", "", colnames(sig_wide))
ord <- unique(sig_counts$method)
sig_wide <- sig_wide[match(ord, sig_wide$method), ]
rownames(sig_wide) <- NULL
print(sig_wide, row.names = FALSE)

write.csv(sig_counts, file.path(DIR_JIN, paste0("jin_sig_counts_fdr_fc_", DATESTAMP, ".csv")),
          row.names = FALSE)

cat("\n== how to read the nested list ==\n")
cat("  sig_summary[[\"cOTU\"]][[\"prox\"]][[\"DESeq2\"]]   # a data.frame\n")
cat("  names(sig_summary)                     ->", paste(names(sig_summary), collapse = ", "), "\n")
cat("  names(sig_summary[[1]])                ->", paste(names(sig_summary[[1]]), collapse = ", "), "\n")
cat("  names(sig_summary[[1]][[1]])           ->", paste(names(sig_summary[[1]][[1]]), collapse = ", "), "\n")

res_by_loc_out <- readRDS(file.path(DIR_JIN, paste0("jin_results_by_location_", DATESTAMP, ".RData")))
View(cbind(res_by_loc_out$prox$sparse_out, res_by_loc_out$dist$sparse_out[rownames(res_by_loc_out$prox$sparse_out),]))


summary_flat <- read.csv(file.path(DIR_JIN, "jin_sig_summary_flat_09012026.csv"))

## these six report VS relative to VD; put them on VD relative to VS like the rest
REVERSED <- c("DESeq2", "edgeR", "LinDA", "ANCOMBC", "ALDEx2_original", "ALDEx2_def")
flip <- summary_flat$method %in% REVERSED
summary_flat$method_lfc[flip] <- -summary_flat$method_lfc[flip]

## each method's own direction; NA for ANCOM and baySeq, which report no effect size
summary_flat$method_direction <- ifelse(summary_flat$method_lfc > 0, "up_in_VD", "up_in_VS")

row.keep <- summary_flat$location == "both" & summary_flat$level == "cOTU"
row.discard <- summary_flat$method %in% c("sparse_SSRV_L", "sparse_SSRV_B", "ALDEx2_original")
col.keep <- c("method", "category", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "method_lfc", "method_direction")
summary_report <- summary_flat[row.keep & (!row.discard), col.keep]
View(summary_report)

write.csv(summary_report, file.path(DIR_JIN, "jin_summary_report.csv"))
