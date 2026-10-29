################################################################################
## FILTER RULE: rowMeans(raw counts) >= 5.   (this run: 09012026)
##
## Reverted 2026-09-01 from the prevalence / relative-abundance rules the August
## scripts used, back to the mean-abundance rule the project used before them.
## ONE rule now decides the feature set in every study except two:
##
##   Prasse, Stammler, Zemb, selex   rowMeans(counts) >= 5   <- this file's family
##   Vandeputte                      20% prevalence, its paper's own rule
##   Jin                             >=5 counts in >=3 samples, its paper's own rule
##
## Vandeputte and Jin are deliberately excluded: both reproduce a published
## analysis whose feature set is fixed by the source paper, so their filters are
## part of what is being reproduced and are not ours to change.
##
## CONSEQUENCE FOR THE PER-METHOD FILTER ARGUMENTS -- the reason this is not a
## one-line change. rowMeans(counts) >= 5 is an ABUNDANCE rule and implies NO
## prevalence floor at all: a genus with every read in a single sample of 12 sits
## at prevalence 0.083 and can still clear a mean of 5. Under the August
## prevalence rules ANCOM-BC2's prv_cut = 0.10 and LinDA's prev.cut = 0.1 were
## provable no-ops; under THIS rule they are not, and would hand two of the
## twelve methods a smaller feature set than the other ten, silently. They are
## therefore 0 here, as are lib_cut / lib.cut. See count_input_policy.R.
##
## The previous version of this file is archived unchanged at
##   archive/preprocessing_2026-08/real_data/
################################################################################

################################################################################
## VARIANT: RAW COUNTS vs QSeq-RESCALED COUNTS, 10% PREVALENCE FILTER
##
## Runs the SAME NAT-vs-REST analysis twice, changing exactly one thing: the
## count matrix every method is fitted on.
##
##   INPUT = "raw"   the observed 16S read counts from TFW.
##   INPUT = "QSeq"  those counts after QSeq(TFW, "QPCR_16S"), i.e. what
##                   sparse_DA_09012026_mean5.R actually feeds the methods.
##
## WHY THIS VARIANT EXISTS
##
## In sparse_DA_09012026_mean5.R, `phylo` is built from Q.TFW (line 282), so
## `Y_raw` at line 382 is the QSeq-rescaled table, NOT raw reads -- the "raw" in
## that name distinguishes it from `Y <- Y_raw + 1` on the next line, i.e. it
## means un-pseudocounted. QSeq multiplies each sample by a factor between
## 1.1e3 and 3.5e6, which has two consequences worth measuring rather than
## arguing about:
##
##   1. aldex.clr() draws from Dirichlet(Y + 0.5). At counts of ~1e9 that
##      posterior collapses to a point for any PRESENT taxon while absent taxa
##      stay pinned at the 0.5 prior, so the Monte Carlo uncertainty ALDEx2
##      exists to propagate is destroyed for exactly the features that carry
##      signal.
##   2. qPCR then enters twice -- once through QSeq, again through the gamma
##      scale model built from log2(QPCR_16S).
##
## Neither is settled by the taxon filter, which is why this script varies the
## input instead. Everything else -- the 10% prevalence cut, the twelve methods,
## the gold reference, CKF, the feature record -- is unchanged.
##
## WHAT IS HELD FIXED ACROSS THE TWO PASSES
##
## QSeq is a per-sample rescaling: after matching genera by name, the zero
## pattern is identical between the two tables and the ratio is constant within
## a sample to 3e-5 (rounding). So the two passes can be made to test literally
## the same features, and they are:
##   * the prevalence filter is computed ONCE, on the RAW genus counts, and the
##     resulting keep-set is applied to both;
##   * tax_glom() picks its archetype OTU by abundance and therefore chooses
##     differently before and after rescaling (543 of 561 genera differ), so the
##     QSeq genus table is re-keyed onto the RAW table's archetype ids below.
##     Without that the two passes would carry different row names and nothing
##     downstream would line up.
##   * scale_mean = log2(QPCR_16S) is the same in both, so the gold reference is
##     targeting the same absolute scale either way.
##
## OUTPUT
##
##   res_tpfp   nested list, res_tpfp[[INPUT]][[gamma_index]], where INPUT is
##              "raw" or "QSeq" and each element is what summarize_tpfp() +
##              add_ckf_row() returned for that gamma -- the same shape
##              sparse_DA_09012026_mean5.R writes, one level deeper.
##
## Every file this script writes carries a _raw or _QSeq suffix, so it cannot
## overwrite anything from sparse_DA_09012026_mean5.R.
################################################################################
################################################################################
## PATHS -- edit ROOT here; nothing below this block hardcodes an absolute path.
##
## ROOT is the project directory that used to be spelled out in full on every
## line. Every path is now built from it with file.path(), so moving between the
## cluster, a laptop, or a copy of the project is a one-line change.
################################################################################
#ROOT <- "/storage/work/wpg5129/sparse_project"   # HPC   (was left active and silently overridden by the line below)
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
R_LIBS <- c("/storage/home/wpg5129/work/R/x86_64-pc-linux-gnu-library/4.4")
R_LIBS <- R_LIBS[dir.exists(R_LIBS)]
if (length(R_LIBS)) .libPaths(c(R_LIBS, .libPaths()))

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
#library(ggpattern) # not working yet
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
library(dirmult)
library(vegan)
library(QSeq)
library(viridis)
library(LinDA)
library(numDeriv)
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

#source(file.path(ROOT, "sparse_SSRV_03122025.R"))
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
## sig_features_matrix() are shared with the Jin analysis, so they live in one
## file rather than being copy-pasted into each script and drifting apart.
source(file.path(ROOT, "feature_record.R"))

## Which column carries each method's own test statistic. Dataset-specific, so it
## stays here rather than in the shared file: the method names and the contrast
## column name (ANCOM-BC's "TYPEREST") are particular to this comparison.
## An explicit NA means the method reports no effect size at all -- ANCOM gives
## detection only, and baySeq gives a direction string plus a log posterior.
## These numbers are NOT sign-comparable across methods; see feature_record.R.
EFFECT_COL <- c(
  ALDEx2_clr  = "effect",           # ALDEx2 standardised effect, NOT a log fold change
  ALDEx2_def  = "effect",
  sparse_SSRV = "mean",             # posterior mean of theta
  ## The same posterior mean, from the laplace-only and bootstrap-only scale
  ## models. They are on the SAME scale as sparse_SSRV and as each other -- the
  ## three differ only in the scale model, so unlike the cross-method statistics
  ## these three ARE directly comparable.
  sparse_SSRV_L = "mean",
  sparse_SSRV_B = "mean",
  DESeq2      = "log2FoldChange",
  edgeR       = "logFC",
  limma       = "logFC",
  baySeq      = NA_character_,
  MiMIX       = "mean",             # posterior mean of beta
  `ANCOM-BC`  = "TYPEREST",
  ANCOM       = NA_character_,
  `ANCOM-BC2` = "lfc_TYPEREST",
  LinDA       = "log2FoldChange"
)


#-------------------------------------------------------------------------------

## Reading in data------------------------------------------------

setwd(file.path(ROOT, "2015_prasse_appliedEnvMicrobiology_wetlands/"))

# import TFW, prepare meta data
TFW<-readRDS("Data/TFW_comdat.RDS")
## fileEncoding is required: TFW_metadata.csv begins with a UTF-8 BOM, so a plain
## read.csv() names the first column "X...ID" instead of "ID". meta$ID is then
## NULL, `rownames(meta) <- NULL` silently resets rownames to the positions
## 1..57, and those do not match TFW's sample names (1-21, 25-60) -- the two
## agree as sets only once the real ID column is used.
meta<-read.csv("Data/TFW_metadata.csv", fileEncoding = "UTF-8-BOM")
stopifnot(!is.null(meta$ID))
rownames(meta)<-meta$ID
sample_data(TFW)<-sample_data(meta)
sample_names(TFW)<-paste("Sample_", sample_names(TFW), sep="")

# total abundance relative standard deviation: ####
sd(sample_data(TFW)$QPCR_16S)/mean(sample_data(TFW)$QPCR_16S) # 0.98

# Qseq transformation ####
Q.TFW<-QSeq(TFW, "QPCR_16S")

## QSeq only rescales counts; it does not add, drop or rename taxa or samples.
## That is what makes the raw/QSeq comparison below a controlled one, so assert
## it rather than assume it.
stopifnot(identical(taxa_names(TFW),   taxa_names(Q.TFW)),
          identical(sample_names(TFW), sample_names(Q.TFW)),
          taxa_are_rows(TFW), taxa_are_rows(Q.TFW))

## sp ids, assigned ONCE to the shared OTU-level taxonomy so both inputs carry
## the same feature names. (The original script assigns these off Q.TFW; since
## the two tables share taxa_names, the same ids apply to TFW unchanged.)
taxa <- Q.TFW@tax_table
taxa_id <- paste0("sp", 1:nrow(taxa))
rownames(taxa) = taxa_id
taxa <- as.matrix(taxa)

sam1 <- Q.TFW@sam_data
tax1 <- tax_table(taxa)

## Two OTU tables, samples x taxa, same column order (taxa_are_rows = FALSE, to
## match the orientation the original script's `phylo` object has).
OTU_qseq <- as.data.frame(t(Q.TFW@otu_table@.Data)); colnames(OTU_qseq) <- taxa_id
OTU_raw  <- as.data.frame(t(TFW@otu_table@.Data));   colnames(OTU_raw)  <- taxa_id

phylo_qseq <- phyloseq(otu_table(OTU_qseq, taxa_are_rows=FALSE), sam1, tax1)
phylo_raw  <- phyloseq(otu_table(OTU_raw,  taxa_are_rows=FALSE), sam1, tax1)

##Amalgamating to genus level
genus_qseq <- tax_glom(phylo_qseq, taxrank = "Genus")
genus_raw  <- tax_glom(phylo_raw,  taxrank = "Genus")

## Re-key the QSeq genus table onto the RAW table's archetype ids ---------------
## tax_glom() names each genus after its most abundant constituent OTU, and
## because QSeq rescales samples by different factors it picks a different
## archetype for 543 of the 561 genera. Matching on Genus NAME and then adopting
## the raw table's ids makes the two passes share row names, ordering and
## tax_table, so every downstream lookup (lineage_of, taxa[,"Genus"], the sig
## matrix, res_tpfp) is comparable between them.
G_qseq <- as.character(tax_table(genus_qseq)[, "Genus"])
G_raw  <- as.character(tax_table(genus_raw)[,  "Genus"])
stopifnot(setequal(G_qseq, G_raw), !anyDuplicated(G_qseq), !anyDuplicated(G_raw))

M_qseq <- as.matrix(otu_table(genus_qseq))            # samples x taxa
colnames(M_qseq) <- G_qseq
M_qseq <- M_qseq[sample_names(genus_raw), G_raw, drop = FALSE]
colnames(M_qseq) <- taxa_names(genus_raw)
genus_qseq <- phyloseq(otu_table(M_qseq, taxa_are_rows = FALSE),
                       sample_data(genus_raw), tax_table(genus_raw))

## The two tables should now differ ONLY by a per-sample scalar. Check it: the
## zero patterns must agree exactly, and within a sample the ratio must be
## constant up to QSeq's rounding.
M_raw_chk <- as.matrix(otu_table(genus_raw))
stopifnot(identical(dimnames(M_qseq), dimnames(M_raw_chk)),
          identical(M_qseq == 0, M_raw_chk == 0))
.ratio <- M_qseq / M_raw_chk; .ratio[!is.finite(.ratio)] <- NA
stopifnot(max(apply(.ratio, 1, function(v) sd(v, na.rm=TRUE)/mean(v, na.rm=TRUE)),
              na.rm = TRUE) < 1e-3)
message(sprintf("QSeq per-sample rescaling factor: %.0fx to %.0fx",
                min(rowSums(M_qseq)/rowSums(M_raw_chk)),
                max(rowSums(M_qseq)/rowSums(M_raw_chk))))

## MEAN-ABUNDANCE filter -------------------------------------------------------
## rowMeans(raw counts) >= MIN_MEAN, the project-wide rule (see the header
## block). It replaces the 10% prevalence rule the August version used.
##
## Computed ONCE, on the RAW genus counts, and applied to BOTH inputs. That
## placement matters MORE under this rule than under the old one. Prevalence was
## a property of the observed data and identical between the passes anyway (the
## zero patterns agree exactly); a MEAN is not -- QSeq multiplies each sample by
## its own factor, and those factors span three orders of magnitude here, so
## "mean >= 5" evaluated on the QSeq table would keep all 561 genera and filter
## nothing. Computing it here, once, on genus_raw, is what makes the threshold
## mean what it says AND holds the feature set identical across the two passes,
## which is the entire point of the contrast.
##
## M_raw_chk is samples x taxa (taxa_are_rows = FALSE on these objects), so the
## per-taxon mean is colMeans(), not rowMeans(). Same rule, transposed layout.
##
## PREVALENCE CONSEQUENCE: rowMeans >= 5 implies no prevalence floor, so
## ANCOM-BC2's prv_cut and LinDA's prev.cut are 0 in this file rather than the
## 0.10 / 0.1 that were exact no-ops under the prevalence rule; lib_cut /
## lib.cut go from 1000 to 0 for the same reason.
MIN_MEAN <- 5
keep_ids <- taxa_names(genus_raw)[colMeans(M_raw_chk) >= MIN_MEAN]
message(sprintf("mean raw count >= %g: kept %d of %d genera (min prevalence among kept %.2f)",
                MIN_MEAN, length(keep_ids), ntaxa(genus_raw),
                min(colMeans(M_raw_chk[, keep_ids, drop = FALSE] > 0))))

## Set to TRUE to also run the theta^perp / Dirichlet diagnostic block that the
## original script has between the filter and the gold standard. It is purely
## descriptive, costs an extra aldex.clr(mc.samples = 1000) plus a 1000-iteration
## loop per input, and nothing downstream depends on it.
RUN_THETA_PERP_DIAG <- FALSE

################################################################################
## MAIN LOOP -- one pass per input matrix.
##
## The body below is sparse_DA_09012026_mean5.R from its `indsC <- ...` line
## through res_tpfp, unchanged except that (a) phylo_genus is picked by INPUT at
## the top, and (b) every output path and every `pair =` key carries the INPUT
## suffix so the two passes cannot overwrite each other.
################################################################################
INPUT_TYPES <- c("raw", "QSeq")

res_tpfp        <- setNames(vector("list", length(INPUT_TYPES)), INPUT_TYPES)
tax_scale_all   <- setNames(vector("list", length(INPUT_TYPES)), INPUT_TYPES)
list_tax_all    <- setNames(vector("list", length(INPUT_TYPES)), INPUT_TYPES)
sparse_out_all  <- setNames(vector("list", length(INPUT_TYPES)), INPUT_TYPES)
ckf_by_gamma_all<- setNames(vector("list", length(INPUT_TYPES)), INPUT_TYPES)
feature_rec_all <- setNames(vector("list", length(INPUT_TYPES)), INPUT_TYPES)

for (INPUT in INPUT_TYPES) {

set.seed(1234)   # re-seeded per pass so the two are not offset by each other

SUF     <- paste0("_09012026_mean5_", INPUT)   # file-name suffix
PAIR    <- paste0("mean5_", INPUT)             # key for the *_record.rds files
OUT_DIR <- file.path(ROOT, "2015_prasse_appliedEnvMicrobiology_wetlands")

cat("\n\n################################################################\n")
cat("##  INPUT =", INPUT, "\n")
cat("################################################################\n")

phylo_genus <- prune_taxa(keep_ids, if (INPUT == "raw") genus_raw else genus_qseq)

indsC <- which(sample_data(phylo_genus)$TYPE == "NAT")
indsT <- which(sample_data(phylo_genus)$TYPE == "REST")

Y <- t(otu_table(phylo_genus))
X <- as.character(sample_data(phylo_genus)$TYPE)

# keep.ind <- rowMeans(Y[,indsC]>0) > 0 & rowMeans(Y[,indsT]>0)
#
# Y <- Y[keep.ind,]
red_phylo <- phyloseq(otu_table(t(Y)), sam_data = sample_data(phylo_genus))
## Count matrices for this comparison, defined ONCE.
##   Y_raw : counts as they come out of phyloseq -- for INPUT = "raw" these are
##           the observed reads, for INPUT = "QSeq" the rescaled table. THIS is
##           what every model is fitted on; see count_input_policy.R.
##   Y     : Y_raw + 1, kept only for derived quantities that want a zero-free
##           matrix (log-ratio summaries, prevalence filters written as Y > 1).
##           No model is fitted on it.
Y_raw <- t(otu_table(red_phylo))
Y     <- Y_raw + 1

D_feat <- nrow(Y_raw)
cat(sprintf("  D = %d genera, N = %d samples | counts range %.0f .. %.0f | %.1f%% zeros\n",
            D_feat, ncol(Y_raw), min(Y_raw), max(Y_raw), 100*mean(Y_raw == 0)))

propY <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
for(i in 1:ncol(Y)){
  propY[,i] <- ((Y[,i]+1)/colSums(Y)[i]) # colSums = sequencing depth
}
LFC <- rowMeans(log2(propY[,indsT]),na.rm = T) - rowMeans(log2(propY[,indsC]),na.rm = T)
hist(LFC, breaks = 100, main = paste("Prasse et al., 2015 --", INPUT))

## Estimate of theta^\perp from the data -- descriptive only (see the flag above)
if (RUN_THETA_PERP_DIAG) {
  normalized_GC <- apply(t(otu_table(phylo_genus))+.5, 2, FUN = function(x){x/sum(x)})
  gm_GC <- apply(normalized_GC, 2, FUN = function(x){mean(log2(x))})
  print(mean(-1*gm_GC[indsC]) - mean(-1*gm_GC[indsT]))

  clr.samps <- aldex.clr(t(otu_table(red_phylo)), as.character(X), mc.samples = 1000, denom = "all")
  dir.sams <- clr.samps@dirichletData
  gm_imp <- rep(NA,1000)
  for(i in 1:1000){
    tmp <- matrix(NA,nrow = dim(red_phylo@otu_table)[2], ncol = dim(red_phylo@otu_table)[1])
    for(k in 1:dim(red_phylo@otu_table)[1]){
      tmp[,k] <- dir.sams[[k]][,i]
    }
    tmp_gm <- apply(tmp, 2, FUN = function(x){mean(log2(x))})
    gm_imp[i] <- mean(-1*mean(tmp_gm[indsC])) - mean(-1*mean(tmp_gm[-indsC]))
  }
  print(mean(gm_imp))
}

## Estimate of theta^\perp due to the qPCR measurements. IDENTICAL in both
## passes -- this is the external scale the gold reference targets, and it is
## deliberately not affected by which count matrix the methods were given.
scale_mean <- log2(sample_data(red_phylo)$QPCR_16S)

# Gold Standard (using given total abundance)
## Choosing scale based off of the estimated technical variation from the other methods.
## Doing a sensitivity analysis as well
sen_res <- list()
gamma <- c(1e-3,.1,.25,.5,0.7,1,1.5,2,3,4,5)
tax_scale <- list()
for(j in 1:length(gamma)){
  scale_var <- rep(gamma[j], ncol(Y))

  scale_samples <- matrix(NA, nrow = ncol(Y), ncol = 1000)
  for(i in 1:ncol(Y)){
    scale_samples[i,] <- 2^rnorm(1000, scale_mean[i], scale_var[i])
  }
  clr <- aldex.clr(Y_raw,X, gamma = scale_samples,mc.samples = 1000)
  res <- aldex.ttest(clr) # paired.test = FALSE
  res <- cbind(res, aldex.lfc(clr))
  sen_res[[j]] <- res
  tax_scale[[j]] <- sen_res[[j]] %>%
    rownames_to_column("category") %>%
    mutate(pval = we.ep) %>%
    mutate(padj = we.eBH) %>%
    dplyr::filter(padj < 0.05)
  cat(sprintf("  gold gamma = %-6g -> %d significant\n", gamma[j], nrow(tax_scale[[j]])))
}

saveRDS(tax_scale, file.path(OUT_DIR, paste0("tax_scale", SUF, ".RData")))

# Ordinary ALDEx2
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

## DeSeq2
coldata <- matrix(as.factor(X), ncol = 1)
colnames(coldata) <- "Condition"
rownames(coldata) <- colnames(as.matrix(Y))

## Raw counts, with the size-factor fallback described in count_input_policy.R.
dds <- deseq_dataset_raw(Y_raw@.Data,
                         colData = coldata,
                         design  = ~ Condition,
                         label   = INPUT)

dds <- DESeq(dds)
res <- deseq_results(dds)
tax_deseq <- res %>%
  as.data.frame() %>%
  rownames_to_column("category") %>%
  dplyr::select(category, log2FoldChange, padj, lfcSE) %>%
  mutate(sig = ifelse(padj < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)

##EdgeR
conds <- as.character(sample_data(red_phylo)$TYPE)
coldata <- ifelse(conds == as.character(head(X,1)), 0, 1)
coldata <- matrix(coldata, ncol = 1)

y <- DGEList(counts=as.matrix(Y_raw@.Data),group=c(coldata[,1]))
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tax_edgeR <- topTags(et, n=1600, p.value = 1) %>%
  as.data.frame() %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(FDR < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)

##limma
conds <- as.character(sample_data(red_phylo)$TYPE)
design <- model.matrix(~0+conds)

y <- DGEList(counts=as.matrix(Y_raw@.Data),group=conds)
y <- calcNormFactors(y, method = "TMM")
v <- voom(y,design,plot=FALSE)
fit <- lmFit(v,design)

contr <- makeContrasts(condsREST - condsNAT, levels = colnames(coef(fit)))
contr.fit <- contrasts.fit(fit, contr)
contr.fit <- eBayes(contr.fit)
tax_limma <- topTable(contr.fit, p.value = 1,number=100) %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)

##baySeq
conds <- as.character(sample_data(red_phylo)$TYPE)
coldata <- ifelse(conds == as.character(head(X,1)), 0, 1)
groups <- list(NDE=rep(1, length(conds)),
               DE=as.numeric(coldata))
## The feature ids have to be handed to baySeq through @annotation -- topCounts()
## returns its rows SORTED by posterior likelihood, so labelling them afterwards
## with the UNSORTED taxon order credits every taxon to the wrong feature.
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
## Per-INPUT data directory. The original writes Y/X/Z.csv into Data/ and lets
## julia write beta.tsv back into the same place; running two passes through one
## directory would make the second silently inherit the first's fit.
x <-  as.numeric(X == unique(X)[1])
mimix_dir <- file.path(ROOT, "2015_prasse_appliedEnvMicrobiology_wetlands",
                       "Data", paste0("mimix_", INPUT))
dir.create(mimix_dir, recursive = TRUE, showWarnings = FALSE)
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
## Without this, any failure of the julia call leaves the old beta.tsv in place
## and the lines below report a stale result as if it were this run's.
beta_file <- file.path(mimix_dir, "beta.tsv")
if (file.exists(beta_file)) file.remove(beta_file)

julia.time <- system.time(mimix_status <- system(julia_command))
if (mimix_status != 0 || !file.exists(beta_file)) {
  stop("MiMIX/Julia failed (exit status ", mimix_status,
       ") for INPUT = ", INPUT, "; no beta.tsv written to ", mimix_dir)
}

mimix_beta <- readr::read_tsv(beta_file)
mimix_out <- as.data.frame(cbind(t(apply(mimix_beta, 2, FUN = function(x){quantile(x, c(0.025, 0.975))}))))
mimix_out$mean <- as.numeric(colMeans(mimix_beta))
Sig <- rep("FALSE", nrow(mimix_out))
Sig[(mimix_out$`2.5%` * mimix_out$`97.5%`) > 0] <- "TRUE"

mimix_out <- cbind(mimix_out, Sig = Sig)
rownames(mimix_out) <- rownames(Y)
tax_mimix <- mimix_out[mimix_out$Sig == "TRUE",]
tax_mimix$category <- rownames(tax_mimix)

# Add ANCOMBC
main.factor <- "TYPE"
ancombc.time <- system.time(ancbcfit <- ancombc(data = red_phylo, assay_name = "counts",
                                                tax_level = NULL,
                                                formula = main.factor,
                                                p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                                                group = main.factor, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                max_iter = 100, conserve = TRUE, alpha = 0.05,
                                                n_cl = 1, verbose = TRUE))

ancbc_out <- cbind(ancbcfit$res$lfc[,c(1,3)], W = ancbcfit$res$W[,3], adj.p = ancbcfit$res$q_val[,3], sig = ancbcfit$res$diff_abn[,3])

rownames(ancbc_out) <- rownames(Y)
ancbc_out <- cbind(ancbc_out, taxa = taxa[,"Genus"][match(rownames(ancbc_out), rownames(taxa))])
tax_ancbc <- ancbc_out[ancbc_out$sig == "TRUE",]
tax_ancbc$category <- rownames(tax_ancbc)

# Add ANCOM
ancom.time <- system.time(ancfit <- ancom(data = red_phylo, assay_name = "counts",
                                          tax_level = NULL,
                                          p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                                          main_var = main.factor,
                                          rand_formula = NULL, lme_control = NULL, struc_zero = FALSE,
                                          neg_lb = TRUE, alpha = 0.05, n_cl = 2))

anc_out <- data.frame(sig = ancfit$res$detected_0.7)
rownames(anc_out) <- ancfit$res$taxon
anc_out <- cbind(anc_out, taxa = taxa[,"Genus"][match(rownames(anc_out), rownames(taxa))])
tax_anc <- anc_out[anc_out$sig == "TRUE",]
tax_anc$category <- rownames(tax_anc)

# Add ANCOMBC2
ancombc2.time <- system.time(ancbc2fit <- ancombc2(data = red_phylo, assay_name = "counts",
                                                   fix_formula = main.factor, rand_formula = NULL,
                                                   p_adj_method = "BH", pseudo_sens = TRUE,
                                                   prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                                   group = main.factor, struc_zero = FALSE, neg_lb = FALSE,
                                                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                   global = FALSE, pairwise = TRUE,
                                                   dunnet = FALSE, trend = FALSE,
                                                   iter_control = list(tol = 1e-5, max_iter = 20,
                                                                       verbose = FALSE),
                                                   em_control = list(tol = 1e-5, max_iter = 100),
                                                   lme_control = NULL,
                                                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                                   trend_control = NULL))

## Keep the log fold change too -- reducing ancombc2 to a detection flag throws
## away the only quantity that says which group a taxon is enriched in.
ancbc2_lfc_col <- grep("^lfc_TYPE", colnames(ancbc2fit$res), value = TRUE)[1]
ancbc2_out <- data.frame(
  taxon        = as.character(ancbc2fit$res$taxon),
  lfc_TYPEREST = if (is.na(ancbc2_lfc_col)) NA_real_ else as.numeric(ancbc2fit$res[[ancbc2_lfc_col]]),
  sig          = ancbc2fit$res$diff_TYPEREST,
  stringsAsFactors = FALSE)

rownames(ancbc2_out) <- ancbc2fit$res$taxon
tax_ancbc2 <- ancbc2_out[ancbc2_out$sig == "TRUE",]
tax_ancbc2$category <- rownames(tax_ancbc2)

## LinDA
meta_linda <- as.matrix(X)
colnames(meta_linda) <- main.factor
linda.time <- system.time(linda.obj <- linda(Y_raw, meta_linda, formula = paste0("~",main.factor), alpha = 0.05,
                                             prev.cut = 0, lib.cut = 0, winsor.quan = 0.97))

linda_out <- linda.obj$output[[1]]
linda_out <- cbind(linda_out, taxa = taxa[,"Genus"][match(rownames(linda_out), rownames(taxa))])
tax_linda <- linda_out[linda_out$reject == "TRUE",]
tax_linda$category <- rownames(tax_linda)

# sparse SSRV
## sparse SSRV takes Y_raw, NOT the pseudocounted Y that some derived quantities
## use: Y is Y_raw + 1, so passing it here would make the Dirichlet concentration
## Y_raw + 1.5 instead of Y_raw + 0.5.
Y_0 <- Y_raw
X_0 <- as.numeric(X == tail(X,1))

S <- 10000
D <- nrow(Y_0)
N <- ncol(Y_0)

propY <- matrix(NA, nrow = D, ncol = N)
for(i in 1:N){
  propY[,i] <- ((Y[,i]+1)/colSums(Y)[i]) # colSums = sequencing depth
}
LFC <- rowMeans(log2(propY[,X_0 == 1]),na.rm = T) - rowMeans(log2(propY[,X_0 == 0]),na.rm = T)
hist(LFC, breaks = 100, main = paste("TFW --", INPUT))

resam.array <- array(NA, dim = c(S,D,N))
sparse_time1 <- system.time(
  for(n in 1:N){
    dirichlet_samples <- rdirichlet(S, as.vector(Y_0[,n]+0.5))
    resam.array[, , n] <- (log2(dirichlet_samples))
  }
)

sparse_time <- system.time(
  sparse_out <- sparse_auto(dat.array = resam.array, Y_0, X_0, D, N, S,
                            filtered.ind = rep(FALSE, D))
)

## `pair` names THIS PASS, so the raw and QSeq rows sit next to each other in
## scale_model_record.rds instead of one replacing the other.
record_scale_model(sparse_out, dataset = "Prasse2015",
                   dir = OUT_DIR, pair = PAIR, S = S, N = N)
saveRDS(sparse_out, file.path(OUT_DIR, paste0("sparse_out", SUF, ".RData")))

## All THREE scale-model results from this one fit: the hybrid (what the
## selection rule picked) plus the laplace-only and bootstrap-only arms, so
## "which genera would the UNSELECTED scale model have called?" is answerable
## without refitting.
arm_report <- report_arm_features(sparse_out, "Prasse2015", pair = PAIR)
tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")
tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")

sparse_sens <- check_sparse_bw_sens(dat.array = resam.array, Y_0, X_0, D, N, S,
                                    bandwidth = c("nrd","SJ-ste"),
                                    filtered.ind = rep(FALSE, D))
saveRDS(sparse_sens, file.path(OUT_DIR, paste0("sparse_sens", SUF, ".RData")))
rm(resam.array); gc(verbose = FALSE)

tax_ssrv_qval <- sparse_out[sparse_out$sig.q == "TRUE",]
tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
colnames(tax_ssrv_qval)[4] <- "Sig"

cat(sprintf("  sparse_SSRV: hybrid %d (arm = %s) | laplace %d | bootstrap %d | lap_fail %d\n",
            nrow(tax_ssrv_qval), attr(sparse_out, "scale_model"),
            nrow(tax_ssrv_L), nrow(tax_ssrv_B), attr(sparse_out, "n_laplace_fail")))

## Results Summary ----------------------------------------------------------------

list_tax <- list(ALDEx2_clr = tax_clr, ALDEx2_def = tax_default,
                 sparse_SSRV = tax_ssrv_qval,
                 DESeq2 = tax_deseq,
                 edgeR = tax_edgeR, limma = tax_limma, baySeq = tax_baySeq, MiMIX = tax_mimix,
                 `ANCOM-BC` = tax_ancbc, ANCOM = tax_anc, `ANCOM-BC2` = tax_ancbc2, LinDA = tax_linda)

saveRDS(list_tax, file.path(OUT_DIR, paste0("results", SUF, ".RData")))

## Feature-level record ---------------------------------------------------------
## Not how MANY genera each method called but WHICH ones, with the lineage
## attached, so the contrast (NAT = natural, REST = restored) reads biologically.
##   POSITIVE = higher in REST (restored).
##   lfc_prop : log2 fold change in RELATIVE abundance.
##   lfc_abs  : the same with the measured 16S load added back.
grpC <- which(X == "NAT")     # control   = natural
grpT <- which(X == "REST")    # treatment = restored
stopifnot(length(grpC) > 0, length(grpT) > 0,
          length(grpC) + length(grpT) == ncol(Y))

Ycnt <- as.matrix(Y_raw@.Data)          # the matrix the methods saw, this pass
## Prevalence and mean abundance are reported off the OBSERVED reads in BOTH
## passes. In the original script Yobs was just a second copy of Ycnt, which on
## the QSeq table makes "mean_count" a rescaled quantity in the hundreds of
## thousands and not comparable to anything. Using the raw genus counts here
## keeps these two columns identical between the raw and QSeq records, which is
## what makes the two tax_lookup tables joinable.
Yobs <- t(M_raw_chk)[rownames(Ycnt), colnames(Ycnt), drop = FALSE]

log_prop <- log2(sweep(Ycnt, 2, colSums(Ycnt), "/"))
lfc_prop <- rowMeans(log_prop[, grpT, drop = FALSE]) -
            rowMeans(log_prop[, grpC, drop = FALSE])
lfc_abs  <- lfc_prop + (mean(scale_mean[grpT]) - mean(scale_mean[grpC]))

tax_lookup <- data.frame(
  lineage_of(rownames(Ycnt), taxa),
  lfc_prop   = as.numeric(lfc_prop),
  lfc_abs    = as.numeric(lfc_abs),
  direction  = ifelse(lfc_abs > 0, "up_in_REST", "up_in_NAT"),
  mean_count = as.numeric(rowMeans(Yobs)),
  prev_NAT   = as.numeric(rowMeans(Yobs[, grpC, drop = FALSE] > 0)),
  prev_REST  = as.numeric(rowMeans(Yobs[, grpT, drop = FALSE] > 0)),
  stringsAsFactors = FALSE, check.names = FALSE)

## Gold reference model: the scale-aware ALDEx2 fit, one table per gamma.
stopifnot(length(tax_scale) == length(gamma))
gold_names <- sprintf("gold_gamma_%s", gamma)
names(tax_scale) <- gold_names

sig_long    <- sig_features_long(list_tax,  taxa, EFFECT_COL)   # method x significant genus
arm_list    <- list(sparse_SSRV_L = tax_ssrv_L, sparse_SSRV_B = tax_ssrv_B)
arm_long    <- sig_features_long(arm_list,  taxa, EFFECT_COL)   # arm x significant genus
gold_long   <- sig_features_long(tax_scale, taxa)               # gamma x significant genus
record_long <- rbind(sig_long, arm_long, gold_long)
stopifnot(!is.null(record_long), nrow(record_long) > 0)

## Union of everything anyone called, one column per method then one per gamma.
## n_methods counts the METHODS only, so it stays a measure of cross-method
## agreement and is not inflated by the gold columns or by the two arms (which
## are two more views of ONE method).
method_names <- names(list_tax)
arm_names    <- names(arm_list)
sig_matrix <- sig_features_matrix(record_long, tax_lookup,
                                  method_levels = c(method_names, arm_names, gold_names))
sig_matrix$n_methods <- rowSums(sig_matrix[, method_names, drop = FALSE])
sig_matrix <- sig_matrix[order(-sig_matrix$n_methods, sig_matrix$category), ]
rownames(sig_matrix) <- NULL

feature_record <- list(tax_lookup  = tax_lookup,
                       sig_long    = sig_long,
                       gold_long   = gold_long,
                       record_long = record_long,
                       sig_matrix  = sig_matrix,
                       gamma       = gamma,
                       input       = INPUT,
                       ## the incidence columns, named explicitly, so the
                       ## summary below does not have to guess which columns of
                       ## sig_matrix are methods and which are lineage metadata
                       method_levels = c(method_names, arm_names, gold_names))

saveRDS(feature_record, file.path(OUT_DIR, paste0("prasse_feature_record", SUF, ".RData")))
write.csv(tax_lookup,  file.path(OUT_DIR, paste0("prasse_tax_lookup",  SUF, ".csv")), row.names = FALSE)
write.csv(record_long, file.path(OUT_DIR, paste0("prasse_sig_long",    SUF, ".csv")), row.names = FALSE)
write.csv(sig_matrix,  file.path(OUT_DIR, paste0("prasse_sig_matrix",  SUF, ".csv")), row.names = FALSE)

cat("\n== significant genera per method (INPUT =", INPUT, ") ==\n")
print(data.frame(
  Method = c(method_names, arm_names, gold_names),
  n      = as.integer(colSums(sig_matrix[, c(method_names, arm_names, gold_names), drop = FALSE])),
  row.names = NULL))

## ---- CKF -------------------------------------------------------------------
## Fitted ONCE on Y_0 / X_0 -- the same matrix and contrast sparse_auto() got --
## and then scored against each gamma's reference set. The knockoff filter does
## not depend on the truth, so refitting per gamma would be wasted fits.
ckf_fit <- run_ckf(Y_0, X_0, truth = tax_scale[[4]]$category,
                   R = CKF_R, q = CKF_Q, offset = CKF_OFFSET,
                   se_target = CKF_SE_TARGET, R_min = CKF_R_MIN,
                   derand_thresh = CKF_DERAND, seed = 1234,
                   dataset = "Prasse2015", pair = PAIR)
print(ckf_stability(ckf_fit, se_target = 0.1))

res_tpfp_input <- list()
ckf_by_gamma   <- vector("list", length(gamma))
for(j in 1:length(gamma)){
  res_tpfp_input[[j]] <- summarize_tpfp(tax_gold = tax_scale[[j]], list_tax)
  ## Re-scoring only; no refit.
  ckf_by_gamma[[j]] <- score_ckf(ckf_fit, tax_scale[[j]]$category,
                                 pair = paste0(PAIR, "_gamma", j))
  res_tpfp_input[[j]] <- add_ckf_row(res_tpfp_input[[j]], ckf_by_gamma[[j]])
}

saveRDS(ckf_by_gamma,     file.path(OUT_DIR, paste0("prasse_ckf_sensitivity", SUF, ".RData")))
saveRDS(ckf_fit$selected, file.path(OUT_DIR, paste0("prasse_ckf_taxa",        SUF, ".RData")))
## gamma[4] is the reported one, so that is the row the keyed record carries.
record_ckf(score_ckf(ckf_fit, tax_scale[[4]]$category, pair = PAIR), dir = OUT_DIR)

## ---- stash this pass --------------------------------------------------------
res_tpfp[[INPUT]]         <- res_tpfp_input
tax_scale_all[[INPUT]]    <- tax_scale
list_tax_all[[INPUT]]     <- list_tax
sparse_out_all[[INPUT]]   <- sparse_out
ckf_by_gamma_all[[INPUT]] <- ckf_by_gamma
feature_rec_all[[INPUT]]  <- feature_record

}   ## end for (INPUT in INPUT_TYPES)

################################################################################
## OUTPUT
##
## res_tpfp[["raw"]][[j]] and res_tpfp[["QSeq"]][[j]], j indexing gamma. Each
## element is exactly what sparse_DA_09012026_mean5.R stores in its res_tpfp[[j]]
## -- a list with an $aldex_gold matrix of TP/FP, one row per method plus the CKF
## row add_ckf_row() appended -- so anything that already reads that file works
## after one extra [[INPUT]].
################################################################################
saveRDS(res_tpfp,         file.path(OUT_DIR, "res_tpfp_09012026_mean5_rawVSqseq.RData"))
saveRDS(tax_scale_all,    file.path(OUT_DIR, "tax_scale_09012026_mean5_rawVSqseq.RData"))
saveRDS(list_tax_all,     file.path(OUT_DIR, "results_09012026_mean5_rawVSqseq.RData"))
saveRDS(sparse_out_all,   file.path(OUT_DIR, "sparse_out_09012026_mean5_rawVSqseq.RData"))
saveRDS(feature_rec_all,  file.path(OUT_DIR, "prasse_feature_record_09012026_mean5_rawVSqseq.RData"))

################################################################################
## SIDE-BY-SIDE SUMMARY
################################################################################

## 1. Significant genera per method, at the headline gamma[4] = 0.5 -------------
##    (gamma is the same vector in both passes, so gamma[4] resolves either way.)
gamma_ref <- 4L
n_by_method <- function(input) {
  fr <- feature_rec_all[[input]]
  cols <- fr$method_levels
  setNames(as.integer(colSums(fr$sig_matrix[, cols, drop = FALSE])), cols)
}
## Both passes use the same twelve methods, two arms and four gammas, so the
## column sets match; assert it rather than silently union them.
stopifnot(identical(feature_rec_all$raw$method_levels,
                    feature_rec_all$QSeq$method_levels))
summary_counts <- data.frame(
  Method = feature_rec_all$raw$method_levels,
  raw    = as.integer(n_by_method("raw")),
  QSeq   = as.integer(n_by_method("QSeq")),
  row.names = NULL)

cat("\n\n=========== significant genera per method: raw vs QSeq input ===========\n")
print(summary_counts, row.names = FALSE)

## 2. The gold reference across gamma -------------------------------------------
gold_counts <- data.frame(
  gamma = gamma,
  raw   = sapply(tax_scale_all$raw,  nrow),
  QSeq  = sapply(tax_scale_all$QSeq, nrow),
  row.names = NULL)
cat("\n=========== gold reference (scale-aware ALDEx2) by gamma ===========\n")
print(gold_counts, row.names = FALSE)

## 3. sparse SSRV: hybrid, both arms, and how wide the posterior actually is -----
ssrv_summary <- do.call(rbind, lapply(INPUT_TYPES, function(input) {
  so <- sparse_out_all[[input]]
  ns <- attr(so, "n_sig")
  data.frame(input       = input,
             arm_chosen  = attr(so, "scale_model"),
             scale_rule  = attr(so, "scale_rule"),
             hybrid      = sum(so$sig.q == "TRUE", na.rm = TRUE),
             laplace     = as.integer(ns[["laplace"]]),
             bootstrap   = as.integer(ns[["bootstrap"]]),
             n_lap_fail  = attr(so, "n_laplace_fail"),
             n_var_reject= attr(so, "n_var_reject"),
             median_CI_width = round(median(so[["97.5%"]] - so[["2.5%"]], na.rm = TRUE), 2),
             row.names = NULL)
}))
cat("\n=========== sparse SSRV: hybrid vs each scale-model arm ===========\n")
print(ssrv_summary, row.names = FALSE)

## 4. TP/FP against each pass's OWN gold, at gamma[4] ---------------------------
##    Each column is scored against the reference built on the SAME input, which
##    is the like-for-like comparison; a method's TP moving between the two
##    columns is the input changing both the method and its reference.
tpfp_ref <- do.call(rbind, lapply(INPUT_TYPES, function(input) {
  m <- as.data.frame(res_tpfp[[input]][[gamma_ref]]$aldex_gold)
  data.frame(input = input, Method = rownames(m), TP = m$TP, FP = m$FP, row.names = NULL)
}))
tpfp_wide <- merge(
  subset(tpfp_ref, input == "raw",  select = c("Method","TP","FP")),
  subset(tpfp_ref, input == "QSeq", select = c("Method","TP","FP")),
  by = "Method", all = TRUE, suffixes = c("_raw", "_QSeq"))
cat("\n=========== TP/FP vs each pass's own gold, gamma =", gamma[gamma_ref], "===========\n")
print(tpfp_wide, row.names = FALSE)

write.csv(summary_counts, file.path(OUT_DIR, "prasse_rawVSqseq_method_counts.csv"), row.names = FALSE)
write.csv(gold_counts,    file.path(OUT_DIR, "prasse_rawVSqseq_gold_by_gamma.csv"), row.names = FALSE)
write.csv(ssrv_summary,   file.path(OUT_DIR, "prasse_rawVSqseq_ssrv_arms.csv"),     row.names = FALSE)
write.csv(tpfp_wide,      file.path(OUT_DIR, "prasse_rawVSqseq_tpfp.csv"),          row.names = FALSE)

cat("\nDone. res_tpfp[[\"raw\"]] and res_tpfp[[\"QSeq\"]] are in the workspace and saved to\n",
    file.path(OUT_DIR, "res_tpfp_09012026_mean5_rawVSqseq.RData"), "\n")
