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
# R_LIBS <- c("/storage/home/wpg5129/R/x86_64-pc-linux-gnu-library/4.4",
#             "/storage/home/wpg5129/work/R/x86_64-pc-linux-gnu-library/4.4")
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
# Install the latest version from github 
#devtools::install_github("jsilve24/driver")
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
library(LinDA)
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

## THE Zemb parser -- one definition, shared. Drops the no-feces control tube
## (oz1802-19) at the source, so counts, metadata and qPCR are subset together
## and every script that parses this dataset gets the same 9 samples. See
## zemb_parse.R for the tube-by-tube reading of the source paper's design.
## THIS script drops nothing further: the 9 samples the parser returns are the
## 9 that every method is fitted on. The block below records why the 100 ul dose
## was excluded until 2026-09-01 and what that exclusion cost.
source(file.path(ROOT, "zemb_parse.R"))
#-------------------------------------------------------------------------------

## Reading in data------------------------------------------------

out <- parse_2020_zemb_microbiologyOpen_spike()

################################################################################
## SAMPLE EXCLUSION -- the 100 ul E. coli dose (oz1802-13).
##
## Experiment 2 is a dose series: 0 / 1 / 5 / 10 / 100 ul of the SAME 10^7
## cells/ul suspension, one replicate of each sequenced (tubes 1/4/7/10/13). The
## contrast every method in this script is given is BINARY -- added_Coli > 0 --
## so all four spiked tubes are summarised by a single group effect.
##
## oz1802-13 is poorly described by that effect. It received TEN TIMES the
## next-highest dose, and the deposited count table shows what that does: the
## spiked E. coli (Seq1) is 1.7%, 9.1% and 20.4% of the reads in tubes 4, 7 and
## 10, and 66.5% in tube 13. A taxon at two thirds of the library compresses
## every other taxon in the composition, so this one sample supplies most of the
## compositional dilution the pooled contrast sees. It is an intentional extreme
## dose, not a defective observation -- it is dropped because the binary
## contrast cannot represent it, not because anything went wrong in the tube.
##
## The exclusion costs signal rather than manufacturing it, which is what makes
## it defensible rather than merely convenient:
##
##   Seq1 relative log2 fold change    5.51  ->  4.69
##   theta^perp measured by qPCR       0.896 ->  0.637   (still positive)
##
## Both the single true positive and the measured load shift survive, smaller.
## What shrinks further is the WITHIN-Experiment-2 part of theta^perp: tube 1
## against tubes 4/7/10 alone is only +0.06, so with 13 gone nearly all of the
## remaining +0.637 comes from the 11-14 mg vs 38-55 mg feces-mass difference
## between the two experiments the arms are drawn from (zemb_parse.R, "STILL
## OPEN"), not from the spike. That confound is neither created nor addressed
## here.
##
## SET TO character(0) ON 2026-09-01 -- the 100 ul tube is BACK IN and this is
## now the full 9-sample analysis. Everything above is the argument for dropping
## it, kept because that is the argument to answer if it is ever dropped again;
## it is no longer the setting in force.
##
## What returns with the tube, from the measurements above: Seq1's relative log2
## fold change goes 4.69 -> 5.51 and the qPCR-measured theta^perp 0.637 -> 0.896.
## Both were shown to SURVIVE the exclusion, so restoring the tube strengthens
## signals that were already there rather than creating new ones. The standing
## objection is unchanged and now simply accepted: a binary added_Coli > 0
## contrast cannot represent a dose ten times the next-highest, and the pooled
## contrast takes most of its compositional dilution from this one sample, where
## Seq1 is 66.5% of the library against 1.7-20.4% in the other three spiked tubes.
##
## REMOVED 2026-09-01. There is no sample-exclusion step in this script any more
## -- no ZEMB_DROP_SAMPLES, no conditional subset, no message. Every tube the
## parser returns is analysed, so this script's design is exactly what
## zemb_parse.R hands it: 9 samples, the no-feces control tube (oz1802-19) having
## already been dropped at the source.
##
## Everything above is retained as the record of why the 100 ul tube was once
## excluded and what excluding it cost. It is the argument to answer if anyone
## proposes dropping it again; it is not describing this script's behaviour.
################################################################################

OTU <- t(out$counts)
dim(OTU)
rowSums(OTU)

#abs_counts <-

#taxa <- out$tax
metadata = out$metadata#[(out$metadata$location == "prox" & out$metadata$sample_type = "ecDNA"),]
metadata <- cbind(metadata, spiked =  metadata$added_Coli > 0)
## The DENSITY (copies per mg), not copies per tube -- see zemb_parse.R for why
## the per-tube total measures the pipetting on this dataset and flips the sign.
## This one assignment sets the units for everything below it: the histograms,
## the theta^perp-from-scale estimate, and the ALDEx2 gold-standard scale model.
scale = out$scale$total16S_per_mg

spiked.ind <- metadata$added_Coli > 0

hist(scale[spiked.ind], breaks = 10, xlim = c(min(scale),max(scale)))
hist(scale[!spiked.ind], breaks = 10, xlim = c(min(scale),max(scale)))

##Creating a phyloseq object
otu1 <- otu_table(OTU, taxa_are_rows=FALSE)
sam1 <- sample_data(metadata) 
#tax1 <- tax_table(as.matrix(taxa))

phylo <- phyloseq(otu1, sam1)

## The no-feces control tube is already gone -- parse_2020_zemb_microbiologyOpen_spike()
## drops it, and drops it from counts, metadata and qPCR together, so `scale` and
## `spiked.ind` above are aligned with these samples by construction rather than by
## a hand-maintained index. zemb_parse.R has the reasoning and the tube list.
stopifnot(length(scale) == nsamples(phylo))

# ##From methods: "taxa unclassified at the genus level or present in less than 20% of samples were excluded from statistical analysis
# #Amalgamating to genus level
# phylo_genus <- tax_glom(phylo)
# 
# #Filtering OTU table in phyloseq
# filter.tax <- function(vec){
#  return(sum(vec > 1) > .2*length(vec))
# }
## FEATURE FILTER: mean raw count >= MIN_MEAN across the samples in the design,
## i.e. rowMeans(counts) >= 5 -- the project-wide rule (see the header block).
##
## This replaces "at least MIN_COUNT raw counts in at least MIN_SAMPLES samples"
## (5 and 2) that the August version used. The two rules are close and they
## differ in BOTH directions, so neither is a relaxation of the other:
##
##   kept by the OLD rule, dropped by this one
##     a feature at exactly 5 counts in 2 of 9 samples and 0 elsewhere: it
##     clears 5 twice, but its mean is 1.11.
##   kept by THIS rule, dropped by the old one
##     a feature at 40 counts in ONE sample and 0 elsewhere: mean 5, but only a
##     single sample reaches 5.
##
## WHY THE SECOND CASE IS TOLERABLE HERE, and it is the reason this study needs
## the note at all. Zemb is a dose series. The old MIN_SAMPLES = 2 was chosen to
## protect dose-dependent features in what was then a 3-sample spiked arm (4 now
## that the 100 ul tube is back) from an all-of-arm requirement. rowMeans >= 5 imposes no sample-count requirement whatsoever, so
## it protects those features a fortiori: a feature responding at only the
## highest dose survives, provided the response is large. The burden moves onto
## multiplicity, where the q-values carry it.
##
## PREVALENCE CONSEQUENCE: the old rule guaranteed prevalence >= 2/9 = 0.22;
## this one guarantees nothing -- a single-sample feature sits at 1/9 = 0.11 and
## passes if its count is big enough. ANCOM-BC2's prv_cut and LinDA's prev.cut
## are therefore 0 in this file rather than 0.10 / 0.1 -- and note how close 0.11
## lands to their old 0.10 default, which is the margin that would have decided
## those features. The message prints the realised minimum so it stays measured
## rather than assumed.
##
## prune_taxa() rather than assigning into @otu_table: a raw slot assignment
## leaves @tax_table at the old size, which is what broke ANCOM-BC2 in the
## Stammler script.
MIN_MEAN <- 5
Mc <- as(otu_table(phylo), "matrix")
if (!taxa_are_rows(phylo)) Mc <- t(Mc)             # taxa x samples
keep_ids <- taxa_names(phylo)[rowMeans(Mc) >= MIN_MEAN]
## Seq135 is the internal standard -- a spiked-in construct, not part of the
## community, so it is not a taxon any method should be testing.
##
## setdiff(), NOT `keep_ids %in% "Seq135"`. The latter turns keep_ids into a
## LOGICAL vector as long as the KEPT set rather than as long as phylo, and
## prune_taxa() rejects that outright ("logical argument to taxa is wrong
## length"), so the script could not run past this line. It also inverted the
## intent -- it selected the standard instead of removing it.
keep_ids  <- setdiff(keep_ids, "Seq135")
red_phylo <- prune_taxa(keep_ids, phylo)
message(sprintf("mean count >= %g: kept %d of %d taxa (min prevalence among kept %.2f)",
                MIN_MEAN, ntaxa(red_phylo), ntaxa(phylo),
                min(rowMeans(Mc[taxa_names(red_phylo), , drop = FALSE] > 0))))

##Estimate of theta^\perp from the data
inds <- which(sample_data(red_phylo)$added_Coli > 0)

normalized_GC <- apply(t(otu_table(red_phylo))+.5, 2, FUN = function(x){x/sum(x)})
gm_GC <- apply(normalized_GC, 2, FUN = function(x){mean(log2(x))})
mean(-1*gm_GC[inds]) - mean(-1*gm_GC[-inds])

# ## Value implied by the dirichlet samples
# ## CLR estimate
# clr.samps <- aldex.clr(t(otu_table(red_phylo)), sample_data(red_phylo)$spiked, mc.samples = 1000, denom = "all")
# dir.sams <- clr.samps@dirichletData
# 
# gm_imp <- rep(NA,1000)
# 
# for(i in 1:1000){
#   tmp <- matrix(NA,nrow = dim(red_phylo@otu_table)[2], ncol = dim(red_phylo@otu_table)[1])
#   for(k in 1:dim(red_phylo@otu_table)[1]){
#     tmp[,k] <- dir.sams[[k]][,i]
#   }
#   tmp_gm <- apply(tmp, 2, FUN = function(x){mean(log2(x))})
#   gm_imp[i] <- mean(-1*mean(tmp_gm[inds])) - mean(-1*mean(tmp_gm[-inds]))
# }
# 
# mean(gm_imp)

## Estimate of theta^\perp from scale
mean(log2(scale[inds])) - mean(log2(scale[-inds]))

## Count matrices for this comparison, defined ONCE.
##   Y_raw : counts as they come out of phyloseq (filtered table, red_phylo).
##   Y     : the single +1-pseudocounted matrix. Nothing downstream adds
##           another +1 -- doing that is what gave DESeq2 counts+2 while
##           edgeR/limma/baySeq got counts+1.
Y_raw <- t(otu_table(red_phylo))
Y     <- Y_raw + 1

## edgeR, limma and baySeq used to read t(otu_table(phylo)) -- the UNFILTERED
## object (428 taxa) -- while ALDEx2 and sparse SSRV used red_phylo (187). That
## was unintended, so they now use Y like every other method and all methods in
## this study are scored on the same taxon set.

X <- sample_data(red_phylo)$added_Coli > 0

D = nrow(Y)
N = ncol(Y)
propY <- matrix(NA, nrow = D, ncol = N)
for(i in 1:N){
  propY[,i] <- ((Y[,i]+1)/colSums(Y)[i]) # colSums = sequencing depth
}
LFC <- rowMeans(log2(propY[,inds])) - rowMeans(log2(propY[,-inds]))
hist(LFC, breaks = 100, main = "zemb spiked vs not spiked")
which.max(LFC)

#-------------------------------------------------------------------------------

## ALDEx2 models----------------------------------------------------------------
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

##Gold standard flow cytometry model

## Choosing scale based off of the estimated technical variation from the other methods. See "Data" for the calculation in the excel file "41586_2017_BFnature24460_MOESM11_ESM.xlsx"
## Doing a sensitivity analysis as well
sen_res <- list()
scale_mean <- log2(scale)
gamma <- c(1e-3,.5)
for(j in 1:length(gamma)){
  scale_var <- rep(gamma[j], ncol(Y))
  
  scale_samples <- matrix(NA, nrow = ncol(Y), ncol = 1000)
  for(i in 1:ncol(Y)){
    scale_samples[i,] <- 2^rnorm(1000, scale_mean[i], scale_var[i])
  }
  clr <- aldex.clr(Y_raw,X>0, gamma = scale_samples,mc.samples = 1000)
  res <- aldex.ttest(clr)
  res <- cbind(res, aldex.lfc(clr))
  sen_res[[j]] <- res
}


## Choosing scale based off of the estimated technical variation from the other methods. See "Data" for the calculation in the excel file "41586_2017_BFnature24460_MOESM11_ESM.xlsx"
tax_scale = sen_res[[1]] %>% 
  rownames_to_column("category") %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)

tax_scale2 = sen_res[[2]] %>% 
  rownames_to_column("category") %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)


#-------------------------------------------------------------------------------

## Other methods ---------------------------------------------------------------

## DeSeq2
coldata <- matrix(sample_data(red_phylo)$spiked, ncol = 1)
colnames(coldata) <- "Condition"
## Raw counts, with the size-factor fallback described in
## count_input_policy.R -- the default median-of-ratios estimator cannot
## cope when every feature has a zero somewhere, which is what the old
## +1 was masking.
dds <- deseq_dataset_raw(Y_raw@.Data,
                         colData = coldata,
                         design  = ~Condition)
dds <- DESeq(dds)
res <- deseq_results(dds)
tax_deseq <- res %>% 
  as.data.frame() %>% 
  rownames_to_column("category") %>% 
  dplyr::select(category, log2FoldChange, padj, lfcSE) %>% 
  mutate(sig = ifelse(padj < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)


##EdgeR

X <- sample_data(red_phylo)$spiked
conds <- sample_data(red_phylo)$spiked 
coldata <- ifelse(conds == "TRUE", 0, 1)
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
X <- sample_data(red_phylo)$spiked
conds <- sample_data(red_phylo)$spiked 
coldata <- ifelse(conds == "TRUE", 0, 1)
coldata <- matrix(coldata, ncol = 1)
design <- model.matrix(~0+conds)

y <- DGEList(counts=as.matrix(Y_raw@.Data),group=conds)
y <- calcNormFactors(y, method = "TMM")
v <- voom(y,design,plot=FALSE)
fit <- lmFit(v,design)

contr <- makeContrasts(condsTRUE - condsFALSE, levels = colnames(coef(fit)))
contr.fit <- contrasts.fit(fit, contr)
contr.fit <- eBayes(contr.fit)
tax_limma <- topTable(contr.fit, p.value = 1,number=100) %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)


##baySeq
conds <- sample_data(red_phylo)$spiked 
coldata <- ifelse(conds == "TRUE", 0, 1)
groups <- list(NDE=rep(1, length(conds)), 
               DE=as.numeric(coldata))
dds <- new("countData", data = as.matrix(Y_raw@.Data), groups = groups, replicates = conds)
libsizes(dds) <- getLibsizes(dds)
dds <- getPriors.NB(dds, cl = NULL)
dds <- getLikelihoods(dds,cl=NULL,bootStraps=3,verbose=FALSE)
res <- topCounts(dds,group="DE", number = 1600)
rownames(res) <- taxa_names(red_phylo)

##Merge to the data to see which taxa because the names are uninformative
post <- as.matrix(dds@posteriors)
colnames(post) = c("NDE", "DEEst")
tax_baySeq <- cbind(res, post) %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(FWER.DE < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)

### MiMIX
x <-  as.numeric(X == "TRUE")
write.csv(t(Y_raw), paste0(file.path(ROOT, "2020_zemb_microOpen_spike/MiMIX"), "/Y.csv"), row.names = F)
write.csv(x, paste0(file.path(ROOT, "2020_zemb_microOpen_spike/MiMIX"), "/X.csv"), row.names = F)
write.csv(1:length(x), paste0(file.path(ROOT, "2020_zemb_microOpen_spike/MiMIX/"), "/Z.csv"), row.names = F)
#write.csv(metadata_df, paste0(output_dir, "/metadata_df.csv"))

src_path <- MIMIX_SRC
res_path <- file.path(ROOT, "2020_zemb_microOpen_spike")
#dir.create(file.path(src_path, "SparseDOSSA-analysis-results/08042024/", paste("results",safe_param_string,i, sep = "_")),recursive=TRUE,showWarnings=TRUE) # change date accordingly

julia_command <- paste(
  shQuote(JULIA), shQuote(file.path(src_path, "scripts/fit-mcmc.jl")),
  "--hyper", shQuote(file.path(src_path, "nutnet-analysis/configs/hyper.yml")),
  "--monitor", shQuote(file.path(src_path, "nutnet-analysis/configs/monitor-mimix-no-factors.yml")),
  "--inits", shQuote(file.path(src_path, "nutnet-analysis/configs/inits.yml")),
  "--factors", "0",
  "--iters", "5000",
  "--chains", "1",
  shQuote(file.path(res_path, "MiMIX")),
  shQuote(file.path(res_path, "MiMIX"))
)

## Delete the previous fit first, then insist Julia actually produced a new one.
## Without this, any failure of the julia call (missing binary, bad MIMIX_SRC,
## unexpanded "~") leaves the old beta.tsv in place and the lines below report a
## stale result as if it were this run's.
mimix_dir <- file.path(res_path, "MiMIX")
beta_file <- file.path(mimix_dir, "beta.tsv")
if (file.exists(beta_file)) file.remove(beta_file)

julia.time <- system.time(mimix_status <- system(julia_command))
if (mimix_status != 0 || !file.exists(beta_file)) {
  stop("MiMIX/Julia failed (exit status ", mimix_status,
       "); no beta.tsv written to ", mimix_dir)
}

mimix_beta <- readr::read_tsv(beta_file)

mimix_out <- as.data.frame(cbind(t(apply(mimix_beta, 2, FUN = function(x){quantile(x, c(0.025, 0.975))}))))
Sig <- rep("FALSE", nrow(mimix_out))
Sig[(mimix_out$`2.5%` * mimix_out$`97.5%`) > 0] <- "TRUE"

mimix_out <- cbind(mimix_out, Sig = Sig)
rownames(mimix_out) <- rownames(Y)
tax_mimix <- mimix_out[mimix_out$Sig == "TRUE",]
tax_mimix$category <- rownames(tax_mimix)
tax_mimix


## LinDA
## red_phylo, not phylo: LinDA was another method reading the UNFILTERED table.
meta <- as.matrix(red_phylo@sam_data$spiked)
colnames(meta) <- "conds"
rownames(meta) <- rownames((red_phylo@otu_table))
linda.time <- system.time(linda.obj <- linda(t(red_phylo@otu_table), meta, formula = '~conds', alpha = 0.05,
                                             prev.cut = 0, lib.cut = 0, winsor.quan = 0.97))

linda_out <- linda.obj$output$conds
tax_linda <- linda_out[linda_out$reject == "TRUE",]
tax_linda$category <- rownames(tax_linda)
tax_linda

# Add ANCOMBC
## red_phylo, not phylo. ANCOM-BC was running on the UNFILTERED table (428 taxa)
## while Y is the filtered one (187) -- which is why
## `rownames(ancbc_out) <- rownames(Y)` below failed with a length mismatch.
## Every other study calls ancombc(data = red_phylo, ...).
ancombc.time <- system.time(ancbcfit <- ancombc(data = red_phylo, assay_name = "counts",
                                                tax_level = NULL, #phyloseq = NULL,
                                                formula = "spiked",
                                                p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                                                group = "spiked", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                max_iter = 100, conserve = TRUE, alpha = 0.05, #global = TRUE,
                                                n_cl = 1, verbose = TRUE))

ancbc_out <- cbind(ancbcfit$res$lfc[,c("taxon","spikedTRUE")], W = ancbcfit$res$W$spikedTRUE, sig = ancbcfit$res$diff_abn[,"spikedTRUE"])

rownames(ancbc_out) <- rownames(Y)
tax_ancbc <- ancbc_out[ancbc_out$sig == "TRUE",]
tax_ancbc$category <- rownames(tax_ancbc)
tax_ancbc

# Add ANCOM
## red_phylo, not phylo -- same unfiltered-table problem as ANCOM-BC above.
ancom.time <- system.time(ancfit <- ancom(data = red_phylo, assay_name = "counts",
                                          tax_level = NULL, #phyloseq = NULL,
                                          p_adj_method = "BH", prv_cut = 0, lib_cut = 0, #prv_cut : omits taxa with prevalence less than given value
                                          main_var = "spiked", #adj_formula = "age + region",
                                          rand_formula = NULL, lme_control = NULL, struc_zero = FALSE,#TRUE,
                                          neg_lb = TRUE, alpha = 0.05, n_cl = 2))

anc_out <- data.frame(W = ancfit$res$W, sig = ancfit$res$detected_0.7)

rownames(anc_out) <- rownames(Y)
tax_anc <- as.data.frame(anc_out[anc_out$sig == "TRUE",])
tax_anc$category <- rownames(tax_anc)
tax_anc


# Add ANCOMBC2
## The hardcoded 13-taxon exclusion that used to sit here
##   ancbc2.inds <- !colnames(red_phylo@otu_table) %in% c("Seq263", "Seq266", ...)
##   red_phylo@otu_table <- red_phylo@otu_table[, ancbc2.inds]
## has been REMOVED (2026-08-27), along with the unguarded first fit and the
## retry loop that followed it. What was wrong with that arrangement:
##
##   * The list was a CACHED error message from an older feature set. It was
##     still fully live under the prevalence filter -- all 13 taxa survive it --
##     so ANCOM-BC2 was silently fitted on 280 features while the other eleven
##     methods used all 293. Verified 2026-08-27: ancombc2() now succeeds on all
##     293, so the exclusion buys nothing and costs comparability. (The
##     prevalence filter drops the 118 single-sample taxa that were the likely
##     source of the zero-variance errors the list was cached from.)
##   * The first fit was NOT wrapped in tryCatch, so the retry loop below it
##     could only ever run after that fit had already succeeded -- and then it
##     refitted the whole model a second time and threw the answer away, because
##     ancbc2_out reads `ancbc2fit` from the first call, not `result`.
##   * That loop also called return() at top level, which is an error in R
##     ("no function to return from"), and had no iteration cap.
##
## What replaces it: ONE fit, wrapped, with the zero-variance recovery the loop
## was reaching for -- prune_taxa() rather than a raw slot assignment (so
## @tax_table cannot fall out of step with @otu_table, which is what broke the
## equivalent line in the Stammler script), a cap on the retries, and the
## dropped taxa recorded in `ancbc2_dropped` so a shrunken feature set is
## visible instead of silent.
ancbc2_dropped <- character(0)
ancombc2.time <- system.time({
  ancbc2fit <- NULL
  for (attempt in 1:10) {
    ancbc2fit <- tryCatch(ancombc2(data = red_phylo, assay_name = "counts",
                                                   fix_formula = "spiked", rand_formula = NULL,
                                                   p_adj_method = "BH", pseudo_sens = TRUE,
                                                   prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                                   group = "spiked", struc_zero = FALSE, neg_lb = FALSE,
                                                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                   global = FALSE, pairwise = TRUE, 
                                                   dunnet = FALSE, trend = FALSE,
                                                   iter_control = list(tol = 1e-5, max_iter = 20, 
                                                                       verbose = FALSE),
                                                   em_control = list(tol = 1e-5, max_iter = 100),
                                                   lme_control = NULL, 
                                                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                                                   trend_control = NULL),
      error = function(e) {
        ## The ONE error this recovers from: ancombc2() refusing to fit taxa with
        ## zero variance. Anything else is re-raised -- a silent catch-all here
        ## would turn a real modelling failure into a quietly shrinking feature
        ## set. Parsed out of the message because the package offers no other way
        ## to learn which taxa it objected to.
        msg <- conditionMessage(e)
        hit <- stringr::str_match(
          msg, "Zero variances have been detected for the following taxa:\\s*(.*)\\n")
        if (is.na(hit[2])) stop("ancombc2() failed for a reason this loop does not handle: ", msg)
        drop <- unlist(strsplit(hit[2], ",\\s*"))
        message("  ancombc2: dropping zero-variance taxa: ", paste(drop, collapse = ", "))
        ancbc2_dropped <<- c(ancbc2_dropped, drop)
        ## prune_taxa(), NOT a raw @otu_table assignment -- the latter leaves
        ## @tax_table at the old size, which is exactly what broke the equivalent
        ## line in the Stammler script.
        red_phylo <<- prune_taxa(setdiff(taxa_names(red_phylo), drop), red_phylo)
        NULL
      })
    if (!is.null(ancbc2fit)) break
  }
  if (is.null(ancbc2fit))
    stop("ancombc2() still failing after 10 zero-variance removals; investigate rather than raising the cap")
})

if (length(ancbc2_dropped)) {
  message(sprintf("ancombc2 was fitted on %d of the %d filtered taxa (%d dropped for zero variance: %s)",
                  ntaxa(red_phylo), ntaxa(red_phylo) + length(ancbc2_dropped),
                  length(ancbc2_dropped), paste(ancbc2_dropped, collapse = ", ")))
} else {
  message(sprintf("ancombc2 was fitted on all %d filtered taxa", ntaxa(red_phylo)))
}

ancbc2_out <- as.data.frame(cbind(taxon = ancbc2fit$res[,c("taxon")], sig = ancbc2fit$res$diff_spikedTRUE))

rownames(ancbc2_out) <- ancbc2fit$res$taxon
tax_ancbc2 <- ancbc2_out[ancbc2_out$sig == "TRUE",]
tax_ancbc2$category <- rownames(tax_ancbc2)
tax_ancbc2

# Sparse SSRV

## sparse SSRV takes the RAW counts (Y_raw), NOT the pseudocounted Y that the

## other methods use. Y is Y_raw + 1, so passing it here would make the Dirichlet

Y_0 <- Y_raw
X_0 <- as.numeric(X == "TRUE")

S <- 10000
D <- nrow(Y_0)
N <- ncol(Y_0)

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

record_scale_model(sparse_out, dataset = "Zemb2020",
                   dir = file.path(ROOT, "2020_zemb_microOpen_spike"),
                   S = S, N = N,
                   ## So a saved record says which samples it was fitted on
                   ## instead of only how many. No sample is excluded any more,
                   ## and the note says so rather than going silent -- a blank
                   ## here would read the same as a record written before the
                   ## exclusion was ever tracked.
                   notes = "all 9 samples; no sample excluded")

## Estimated theta^perp, scored against the qPCR -----------------------------
## theta^perp is the shift that turns this method's relative log fold changes
## into absolute ones; under the model it IS the log2 fold change of the total
## load between spiked and unspiked samples. Zemb measured that total directly by
## qPCR, so the estimate is checkable rather than merely asserted -- which is the
## point of running this study at all.
##
## The total has to be a DENSITY -- copies per mg of feces -- and NOT copies per
## tube. This is not a presentational choice, it is what theta^perp is:
##
##   tube n holds mg_n of feces at density d_dn (copies of taxon d per mg), so it
##   holds a_dn = d_dn * mg_n copies, and the sequenced composition is
##       w_dn = d_dn * mg_n / (mg_n * sum_d d_dn) = d_dn / sum_d d_dn.
##
## The mass CANCELS. The composition carries no information about how much feces
## was loaded, so the shift that turns a relative log fold change into an absolute
## one is the log2 ratio of total DENSITY, not of total copies in the tube.
##
## It matters a great deal here. This design puts 11-14 mg in the spiked tubes and
## 38-55 mg in four of the five unspiked ones -- a -1.32 log2 mass imbalance -- so
## copies_in_tube reports theta^perp = -0.684 where the density reports +0.637
## (on the 8 samples analysed; -0.437 vs +0.896 with the 100 ul dose kept).
## The tube figure is measuring the pipetting, not the biology, and it is the
## reason this study looked like a miss. Jin (copies/mg) and Vandeputte (cells per
## gram) are already densities, so Zemb was the only one on the wrong footing.
##
## `scale` IS the density already -- it is out$scale$total16S_per_mg, divided by
## the feces mass once, in the parser. Do NOT divide by mgfeces again here; this
## block used to do exactly that back when `scale` was the per-tube total, and
## repeating it now would apply the mass correction twice.
##
## It is in red_phylo's sample order and aligns with X_0 by construction, because
## the parser subsets counts, metadata and the qPCR table together; the stopifnot()
## is what stops a future change in the sample filtering from silently comparing
## misaligned vectors. X_0 is 1 for the spiked arm, so the difference is treatment
## minus control -- the direction theta^perp is defined in.
##
## The 0 mg control tube is already gone -- the parser drops it, before the taxa
## filter and before any method saw it -- so every remaining sample has a defined
## density and the measurement covers exactly the samples the model was fitted on.
## The is.finite() guard is a safety net rather than a live filter: it fires only
## if a future sample set reintroduces a zero mass (the parser writes NA there,
## not Inf), and says so instead of quietly averaging over nothing.
stopifnot(length(scale) == N)
zemb_ok <- is.finite(scale)
zemb_perp_obs <- mean(log2(scale[zemb_ok & X_0 == 1])) -
                 mean(log2(scale[zemb_ok & X_0 == 0]))
if (any(!zemb_ok))
  message("  [theta^perp] Zemb: ", sum(!zemb_ok), " sample(s) have no measured ",
          "density (", paste(colnames(Y_0)[!zemb_ok], collapse = ", "), ") -- in ",
          "the fit, out of the measured total")

## THE MEASUREMENT AND THE ESTIMATE ARE NOT OVER THE SAME COMMUNITY.
##
## `scale` is the qPCR total for everything in the tube. theta^perp as fitted is
## the log2 fold change of the total load of the KEPT features only, because
## sparse_auto() is handed proportions renormalised within that set. Writing
## T_S,n = r_n T_n with r_n the retained fraction gives
## theta^perp_S = Delta log2 T_n + Delta log2 r_n exactly, and
## perp_matched_obs() evaluates it by forming W_dn = p_dn * T_n over the whole
## community and summing over the kept features.
##
## Seq135 leaves the denominator as well as the numerator: it is spiked-in
## standard DNA, and the qPCR total is derived FROM it, not inclusive of it.
##
## On Zemb this is worth about -0.005 log2 -- r_n runs 0.98-0.99, because an
## abundance filter drops features that carry almost no reads. It does not move
## this study's bias or its coverage verdict, and recording it is how that stays
## a measured fact rather than an assumption. The studies where it does move
## (Prasse, Vandeputte) are the ones with a tax_glom() step.
##
## `scale` is named so the projection checks its sample order against the counts
## instead of trusting position; the parser's own stopifnot is what makes those
## names trustworthy. Non-finite entries are dropped inside perp_matched_obs(),
## from both halves at once, matching zemb_ok above.
zemb_scale_named <- scale
names(zemb_scale_named) <- colnames(Y_0)
zemb_M_all <- as(otu_table(phylo), "matrix")
if (taxa_are_rows(phylo)) zemb_M_all <- t(zemb_M_all)     # samples x taxa
zemb_M_all <- zemb_M_all[, setdiff(colnames(zemb_M_all), "Seq135"), drop = FALSE]
zemb_perp_matched <- perp_matched_obs(total = zemb_scale_named, X = X_0,
                                      Y_all = zemb_M_all, taxa_are_rows = FALSE,
                                      keep = keep_ids)

perp_report <- record_theta_perp(sparse_out, dataset = "Zemb2020",
                                 dir = file.path(ROOT, "2020_zemb_microOpen_spike"),
                                 S = S, N = N, obs = zemb_perp_obs,
                                 matched = zemb_perp_matched,
                                 notes = paste0("qPCR: average16S copies per mg feces ",
                                                "(density); all 9 samples"))

## All THREE scale-model results from this one fit ----------------------------
## sparse_out is the hybrid rule's answer. sparse_auto() now also returns the
## laplace-only and bootstrap-only fits, so the run log says WHICH taxa the
## unselected scale model would have called, not just how many. Printed here and
## written to sparse_arm_features.csv by record_scale_model() above.
arm_report <- report_arm_features(sparse_out, "Zemb2020")
tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")
tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")
## The spiked E. coli sequence -- this study's single true positive, spelled the
## same way `real data results.R` spells it when it scores the arms.
zemb_truth <- "Seq1"
cat("spike-in", zemb_truth, "called by -- hybrid:",
    zemb_truth %in% rownames(sparse_out)[sparse_out$sig.q %in% TRUE],
    " laplace:",   zemb_truth %in% tax_ssrv_L$category,
    " bootstrap:", zemb_truth %in% tax_ssrv_B$category, "\n")

sparse_sens <- check_sparse_bw_sens(dat.array = resam.array, Y_0, X_0, D, N, S, bandwidth = c("nrd","SJ-ste"),
                                    filtered.ind = rep(FALSE, D))
saveRDS(sparse_sens, file.path(ROOT, "2020_zemb_microOpen_spike/sparse_sens_09012026.RData"))


tax_ssrv <- sparse_out[sparse_out$sig == "TRUE",]
tax_ssrv$category <- rownames(tax_ssrv)
colnames(tax_ssrv)[4] <- "Sig"
tax_ssrv_qval <- sparse_out[sparse_out$sig.q == "TRUE",]
tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
colnames(tax_ssrv_qval)[4] <- "Sig"



## Results Summary ----------------------------------------------------------------

list_tax <- list(ALDEx2_original = tax_clr, ALDEx2_def = tax_default, 
                 sparse_SSRV = tax_ssrv_qval,
                 DESeq2 = tax_deseq, edgeR = tax_edgeR, limma = tax_limma, baySeq = tax_baySeq, MiMIX = tax_mimix, 
                 ANCOMBC = tax_ancbc, ANCOM = tax_anc, ANCOMBC2 = tax_ancbc2, LinDA = tax_linda)


saveRDS(list_tax, file.path(ROOT, "2020_zemb_microOpen_spike/results_09012026.RData"))
list_tax <- readRDS(file.path(ROOT, "2020_zemb_microOpen_spike/results_09012026.RData"))

## TRUTH. Seq1 is the added E. coli spike-in and is the only taxon known to
## differ between the two arms, so it is the truth set for EVERY method in this
## table -- CKF included.
##
## This used to read `tax_gold <- tax_ssrv_qval`, i.e. sparse_SSRV's own
## significant calls. Under that reference the table did not measure recovery of
## the spike-in at all, it measured agreement with sparse_SSRV, which is why
## sparse_SSRV scored 10/0 in it by construction while CKF -- scored against
## Seq1 by its own separate script -- was on a different footing entirely.
tax_gold <- data.frame(category = zemb_truth, stringsAsFactors = FALSE)
res_tpfp <- summarize_tpfp(tax_gold = tax_gold, list_tax)
res_tpfp


## ---- CKF -------------------------------------------------------------------
## Fitted on Y_0 / X_0, the SAME matrix and contrast handed to sparse_auto() a
## few lines above, so CKF is scored on exactly the taxa the other methods were
## scored on. The standalone script this replaces re-parsed the study and
## rebuilt the filter itself, which is how its taxon set drifted from this one.
## TRUTH SET. zemb_truth is "Seq1", the added E. coli spike-in -- the same
## reference summarize_tpfp() is given below, so every row of this table,
## including this one, is scored against the same thing.
ckf <- run_ckf(Y_0, X_0, truth = zemb_truth,
               R = CKF_R, q = CKF_Q, offset = CKF_OFFSET,
               se_target = CKF_SE_TARGET, R_min = CKF_R_MIN,
               derand_thresh = CKF_DERAND, seed = 1234,
               dataset = "Zemb2020", pair = NA_character_)
res_tpfp <- add_ckf_row(res_tpfp, ckf)
record_ckf(ckf, dir = file.path(ROOT, "2020_zemb_microOpen_spike"), file = "zemb_mrckf_out.RData")
## What R would have been enough, from the fits just done.
print(ckf_stability(ckf, se_target = 0.1))

saveRDS(res_tpfp, file.path(ROOT, "2020_zemb_microOpen_spike/res_tpfp_09012026.RData"))
#res_tpfp <- readRDS(file.path(ROOT, "2020_zemb_microOpen_spike/res_tpfp_09012026.RData"))
