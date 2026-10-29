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
# ROOT <- "/storage/work/wpg5129/sparse_project"   # HPC
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

## Selex analysis

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
#library(cmdstanr)
library(numDeriv)
library(dirmult)
## ANCOMBC was never attached here, although this script calls ancombc(),
## ancombc2() and ancom() -- the other five study scripts all load it. A
## standalone run died at the first ancombc() call with "could not find
## function", after MiMIX had already completed.
library(ANCOMBC)
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
## 5, not the 10 every other study uses (2026-09-01). CKF is by far the most
## expensive step in THIS script and in no other: the knockoff fit is on 1600
## features against 14 samples, where the other studies run 27-207 features, and
## a single fit takes tens of minutes. At R = 10 the stage ran past 70 minutes
## without finishing; 5 halves it.
##
## 5 is CKF_R_MIN below -- the floor this script already treats as enough fits to
## report -- so the derandomised (freq >= 0.50) selection still means something.
## The cost is precision, not validity: the Monte Carlo SE on mean TP and mean FP
## widens by about sqrt(2). selex is the only study at this setting, so its CKF
## row carries more MC noise than the others' and should be read that way.
CKF_R         <- 5
CKF_Q         <- 0.05
CKF_OFFSET    <- 0
CKF_SE_TARGET <- NULL
CKF_R_MIN     <- 5
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

#-------------------------------------------------------------------------------

## Reading in the selex data------------------------------------------------

data("selex")
head(selex)

Y <- selex

hist(colSums(Y)[1:7], col = "red", breaks = 7, xlim = c(100000, 750000))
hist(colSums(Y)[8:14], col = "blue", breaks = 7, add = T)

conds <- c(rep("NS", 7), rep("S", 7))
x <- as.numeric(conds == "S")

samples_df <- as.data.frame(cbind(sample = paste0("sample", 1:length(conds)), conds, CellCount = as.numeric(colSums(Y))))
samples_df <- samples_df %>% tibble::column_to_rownames("sample")

otu_mat <- as.matrix(Y)
colnames(otu_mat) <- rownames(samples_df)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
samples = sample_data(samples_df)
phylo <- phyloseq(OTU, sam_data = samples)
phylo_genus <- phylo


## FEATURE FILTER: mean raw count >= MIN_MEAN over all 14 samples, i.e.
## rowMeans(counts) >= 5 -- the project-wide rule (see the header block).
##
## This study previously ran UNFILTERED, on all 1600 sequences, on the argument
## that a SELEX library has no rare-taxon tail to remove. That argument is still
## true; the reason to filter anyway is uniformity. Prasse, Stammler and Zemb all
## decide their feature set with this rule, and a feature count that means
## "cleared mean 5" in three studies and "everything there was" in the fourth is
## not a number the cross-study table can carry.
##
## Evaluated on the RAW counts: phylo_genus here IS the unrescaled selex table
## (this study has no QSeq step), so there is none of the rescaled/raw ambiguity
## the Prasse scripts have to navigate.
##
## PREVALENCE CONSEQUENCE. rowMeans >= 5 implies no prevalence floor in general,
## but on THIS table it cannot bite: the least prevalent of the 1600 sequences is
## already in 7 of 14 samples, so every survivor is at prevalence >= 0.50. The
## message below prints the realised minimum so that stays checked rather than
## assumed. prv_cut / prev.cut / lib_cut / lib.cut are 0 regardless.
MIN_MEAN   <- 5
filter.tax <- function(vec) mean(vec) >= MIN_MEAN
red_phylo  <- filter_taxa(phylo_genus, filter.tax, prune = TRUE)
.Mkeep <- as(otu_table(red_phylo), "matrix")
if (!taxa_are_rows(red_phylo)) .Mkeep <- t(.Mkeep)   # taxa x samples
message(sprintf("mean count >= %g: kept %d of %d features (min prevalence among kept %.2f)",
                MIN_MEAN, ntaxa(red_phylo), ntaxa(phylo_genus),
                min(rowSums(.Mkeep > 0)) / nsamples(red_phylo)))
rm(.Mkeep)

##Estimate of theta^\perp from the data
inds <- which(sample_data(red_phylo)$conds == "S")
normalized_GC <- apply(t(otu_table(red_phylo))+.5, 2, FUN = function(x){x/sum(x)})
gm_GC <- apply(normalized_GC, 2, FUN = function(x){mean(log2(x))})
mean(-1*gm_GC[inds]) - mean(-1*gm_GC[-inds])

## Value implied by the dirichlet samples
## CLR estimate
clr.samps <- aldex.clr((otu_table(red_phylo)), sample_data(red_phylo)$conds, mc.samples = 2000, denom = "all") # removed transpose for otu_table()
dir.sams <- clr.samps@dirichletData

gm_imp <- rep(NA,2000)

for(i in 1:2000){
  tmp <- matrix(NA,nrow = dim(red_phylo@otu_table)[1], ncol = dim(red_phylo@otu_table)[2])
  for(k in 1:dim(red_phylo@otu_table)[2]){
    tmp[,k] <- dir.sams[[k]][,i]
  }
  tmp_gm <- apply(tmp, 2, FUN = function(x){mean(log2(x))}) # used to be log2, but I changed it to log
  gm_imp[i] <- mean(-1*mean(tmp_gm[inds])) - mean(-1*mean(tmp_gm[-inds]))
}

mean(gm_imp)


#-------------------------------------------------------------------------------

## ALDEx2 models----------------------------------------------------------------
## First aldex2
## Count matrices for this comparison, defined ONCE.
##   Y_raw : counts as they come out of phyloseq. THIS is what every model is
##           fitted on -- see count_input_policy.R for method-by-method reasons.
##   Y     : Y_raw + 1, kept only for derived quantities that want a zero-free
##           matrix (log-ratio summaries, prevalence filters written as Y > 1).
##           No model is fitted on it.
## Until 2026-08-19 every count-based method was fed Y. On the Jin data that
## alone was the difference between reproducing the source paper's DESeq2 counts
## (8 proximal / 1 distal) and not (9 / 5).
Y_raw <- (otu_table(red_phylo))
Y     <- Y_raw + 1
X <- as.numeric(sample_data(red_phylo)$conds == "S")
aldex_fit <- aldex(Y_raw,X,mc.samples = 1000)

tax_clr = aldex_fit %>% 
  rownames_to_column("category") %>%
  dplyr::select(category, effect, we.ep, we.eBH) %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)


##Default scale model
aldex_fit <- aldex(Y_raw,X,mc.samples = 2000, gamma = .5)

tax_default = aldex_fit %>% 
  rownames_to_column("category") %>%
  dplyr::select(category, effect, we.ep, we.eBH) %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)

selex_truth <- read.table(file.path(ROOT, "selex/selex-truth.txt"))

sum(rownames(selex_truth) %in% tax_default$category)
sum(!(tax_default$category %in% rownames(selex_truth)))

aldex_def_obj <- aldex.clr(Y_raw,X,mc.samples = 2000, gamma = .5, verbose = FALSE)
scale.samps <- getScaleSamples(aldex_def_obj)
scale.samps
dim(scale.samps)

#-------------------------------------------------------------------------------

## Other methods ---------------------------------------------------------------

## DeSeq2
coldata <- matrix(sample_data(red_phylo)$conds, ncol = 1)
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

X <- sample_data(red_phylo)$conds
conds <- matrix(sample_data(red_phylo)$conds, ncol = 1)
coldata <- ifelse(conds == "S", 0, 1)

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
X <- sample_data(red_phylo)$conds
conds <- sample_data(red_phylo)$conds
coldata <- ifelse(conds == "S", 0, 1)
coldata <- matrix(coldata, ncol = 1)
design <- model.matrix(~0+conds)

y <- DGEList(counts=as.matrix(Y_raw@.Data),group=coldata)
y <- calcNormFactors(y, method = "TMM")
v <- voom(y,design,plot=FALSE)
fit <- lmFit(v,design)

contr <- makeContrasts(condsS - condsNS, levels = colnames(coef(fit)))
contr.fit <- contrasts.fit(fit, contr)
contr.fit <- eBayes(contr.fit)
tax_limma <- topTable(contr.fit, p.value = 1,number=100) %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)


##baySeq
conds <- sample_data(red_phylo)$conds
coldata <- ifelse(conds == "S", 0, 1)
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
src_path <- MIMIX_SRC
res_path <- file.path(ROOT, "selex")

## This script never wrote its MiMIX inputs, and selex/data/ did not exist, so
## the julia call had nothing to read. Y is D x N here (taxa are rows in this
## phyloseq object) and MiMIX wants samples x taxa, hence the transpose.
mimix_in <- file.path(res_path, "data")
dir.create(mimix_in, showWarnings = FALSE, recursive = TRUE)
x <- as.numeric(X == "S")
write.csv(t(Y_raw), file.path(mimix_in, "Y.csv"), row.names = FALSE)
write.csv(x, file.path(mimix_in, "X.csv"), row.names = FALSE)
write.csv(1:length(x), file.path(mimix_in, "Z.csv"), row.names = FALSE)
#dir.create(file.path(src_path, "SparseDOSSA-analysis-results/08042024/", paste("results",safe_param_string,i, sep = "_")),recursive=TRUE,showWarnings=TRUE) # change date accordingly

julia_command <- paste(
  shQuote(JULIA), shQuote(file.path(src_path, "scripts/fit-mcmc.jl")),
  "--hyper", shQuote(file.path(src_path, "nutnet-analysis/configs/hyper.yml")),
  "--monitor", shQuote(file.path(src_path, "nutnet-analysis/configs/monitor-mimix-no-factors.yml")),
  "--inits", shQuote(file.path(src_path, "nutnet-analysis/configs/inits.yml")),
  "--factors", "0",
  "--iters", "5000",
  "--chains", "1",
  shQuote(file.path(res_path, "data")),
  shQuote(file.path(res_path, "MiMIX results"))
)

## Delete the previous fit first, then insist Julia actually produced a new one.
## Without this, any failure of the julia call (missing binary, bad MIMIX_SRC,
## unexpanded "~") leaves the old beta.tsv in place and the lines below report a
## stale result as if it were this run's.
mimix_dir <- file.path(res_path, "MiMIX results")
dir.create(mimix_dir, showWarnings = FALSE, recursive = TRUE)   # julia will not create it
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

library(LinDA)
meta <- as.matrix(red_phylo@sam_data$conds)
colnames(meta) <- "conds"
rownames(meta) <- rownames(t(red_phylo@otu_table))
linda.time <- system.time(linda.obj <- linda((red_phylo@otu_table), meta, formula = '~conds', alpha = 0.05,
                                             prev.cut = 0, lib.cut = 0, winsor.quan = 0.97))

linda_out <- linda.obj$output$conds
tax_linda <- linda_out[linda_out$reject == "TRUE",]
tax_linda$category <- rownames(tax_linda)
tax_linda

# Add ANCOM
# Add ANCOMBC
ancombc.time <- system.time(ancbcfit <- ancombc(data = red_phylo, assay_name = "counts",
                                                tax_level = NULL, #phyloseq = NULL,
                                                formula = "conds",
                                                p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                                                group = "conds", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                max_iter = 100, conserve = TRUE, alpha = 0.05, #global = TRUE,
                                                n_cl = 1, verbose = TRUE))

ancbc_out <- cbind(ancbcfit$res$lfc[,c("taxon","condsS")], W = ancbcfit$res$W$condsS, sig = ancbcfit$res$diff_abn[,"condsS"])

rownames(ancbc_out) <- rownames(Y)
tax_ancbc <- ancbc_out[ancbc_out$sig == "TRUE",]
tax_ancbc$category <- rownames(tax_ancbc)
tax_ancbc

ancom.time <- system.time(ancfit <- ancom(data = red_phylo, assay_name = "counts",
                                          tax_level = NULL, #phyloseq = NULL,
                                          p_adj_method = "BH", prv_cut = 0, lib_cut = 0, #prv_cut : omits taxa with prevalence less than given value
                                          main_var = "conds", #adj_formula = "age + region",
                                          rand_formula = NULL, lme_control = NULL, struc_zero = FALSE,#TRUE,
                                          neg_lb = TRUE, alpha = 0.05, n_cl = 2))

anc_out <- data.frame(W = ancfit$res$W, sig = ancfit$res$detected_0.7)

rownames(anc_out) <- rownames(Y)
tax_anc <- as.data.frame(anc_out[anc_out$sig == "TRUE",])
tax_anc$category <- rownames(tax_anc)
tax_anc

# Add ANCOMBC2
ancbc2.time <- system.time(ancbc2fit <- ancombc2(data = red_phylo, assay_name = "counts",
                                                 fix_formula = "conds", rand_formula = NULL,
                                                 p_adj_method = "BH", pseudo_sens = TRUE,
                                                 prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                                 group = "conds", struc_zero = FALSE, neg_lb = FALSE,
                                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                 global = FALSE, pairwise = TRUE, 
                                                 dunnet = FALSE, trend = FALSE,
                                                 iter_control = list(tol = 1e-5, max_iter = 20, 
                                                                     verbose = FALSE),
                                                 em_control = list(tol = 1e-5, max_iter = 100),
                                                 lme_control = NULL, 
                                                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                                                 trend_control = NULL))

ancbc2_out <- cbind(ancbc2fit$res[,c("taxon","diff_condsS")], sig = ancbc2fit$res$diff_condsS)

rownames(ancbc2_out) <- rownames(Y)
tax_ancbc2 <- ancbc2_out[ancbc2_out$sig == "TRUE",]
tax_ancbc2$category <- rownames(tax_ancbc2)
tax_ancbc2


#-------------------------------------------------------------------------------

## Sparse model ----------------------------------------------------------------
## Y is already the +1 matrix defined above; do not re-derive it raw here.
X <- as.numeric(sample_data(red_phylo)$conds == "S")

S <- 20000
## sparse SSRV runs on the UNPADDED counts. Y is Y_raw + 1, a pseudocount added
## for the competitor methods above (ALDEx2/DESeq2/edgeR/baySeq); feeding it here
## too would make the Dirichlet concentration Y_raw + 1.5 instead of Y_raw + 0.5,
## over-smoothing the per-taxon posteriors.
D <- nrow(Y_raw)
N <- ncol(Y_raw)

resam.array <- array(NA, dim = c(S,D,N))
sparse_time1 <- system.time(
  for(n in 1:N){
    dirichlet_samples <- rdirichlet(S, as.vector(Y_raw[,n] + 0.5))
    resam.array[, , n] <- (log2(dirichlet_samples))
  }
)
sparse_out <- sparse_auto_sens(dat.array = resam.array, Y_raw, X, D, N, S, bandwidth = "nrd0", file_path = NULL,
                               filtered.ind = rep(FALSE, D))

record_scale_model(sparse_out, dataset = "McMurrough2014",
                   dir = file.path(ROOT, "selex"),
                   S = S, N = N, notes = "sparse_auto_sens, bandwidth=nrd0")

## All THREE scale-model results from this one fit ----------------------------
## sparse_out is the hybrid rule's answer. sparse_auto() now also returns the
## laplace-only and bootstrap-only fits, so the run log says WHICH genes the
## unselected scale model would have called, not just how many. Printed here and
## written to sparse_arm_features.csv by record_scale_model() above.
arm_report <- report_arm_features(sparse_out, "McMurrough2014")
tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")
tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")
cat("TP/FP against selex-truth --",
    " laplace:",   sum(tax_ssrv_L$category %in% rownames(selex_truth)), "/",
                   sum(!tax_ssrv_L$category %in% rownames(selex_truth)),
    " bootstrap:", sum(tax_ssrv_B$category %in% rownames(selex_truth)), "/",
                   sum(!tax_ssrv_B$category %in% rownames(selex_truth)), "\n")

tax_ssrv_qval <- sparse_out[sparse_out$sig.q == "TRUE",]
tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
colnames(tax_ssrv_qval)[4] <- "Sig"
dim(tax_ssrv_qval)

sum(tax_ssrv_qval$category %in% rownames(selex_truth))
sum(!tax_ssrv_qval$category %in% rownames(selex_truth))

#sparse_out_bws <- sparse_auto_sens_multi(dat.array = resam.array, Y, X, D, N, S, filtered.ind = rep(FALSE,D))

sparse_sens <- check_sparse_bw_sens(dat.array = resam.array, Y_raw, X, D, N, S, bandwidth = c("nrd","SJ-ste"), file_path = file.path(ROOT, "selex/"),
                                    filtered.ind = rep(FALSE, D))
saveRDS(sparse_sens, file.path(ROOT, "selex/sparse_sens_09012026.RData"))

# propY <- matrix(NA, nrow = D, ncol = N)
# for(i in 1:N){
#   propY[,i] <- ((Y[,i]+1)/colSums(Y)[i]) # colSums = sequencing depth
# }
# LFC <- rowMeans(log2(propY[,X==1])) - rowMeans(log2(propY[,X==0]))
# hist(LFC, main = "SELEX", breaks = 100)
# abline(v = LFC[rownames(Y) %in% rownames(selex_truth)], col = "red")
# abline(v = LFC[rownames(Y) %in% tax_ssrv_qval$category], col = "blue")

#-------------------------------------------------------------------------------

## Results Summary ----------------------------------------------------------------

list_tax <- list(ALDEx2_original = tax_clr, ALDEx2_def = tax_default, 
                 sparse_SSRV = tax_ssrv_qval,
                 DESeq2 = tax_deseq, edgeR = tax_edgeR, limma = tax_limma, baySeq = tax_baySeq, MiMIX = tax_mimix, 
                 `ANCOM-BC` = tax_ancbc, ANCOM = tax_anc, `ANCOM-BC2` = tax_ancbc2, LinDA = tax_linda)
tax_selex_truth <- as.matrix(rownames(selex_truth))
colnames(tax_selex_truth) <- "category"
res_tpfp <- summarize_tpfp(as.data.frame(tax_selex_truth), list_tax)
res_tpfp

saveRDS(list_tax, file.path(ROOT, "selex/results_09012026.RData")) 

## ---- CKF -------------------------------------------------------------------
## Fitted on Y_raw / X, the SAME matrix and contrast handed to sparse_auto() a
## few lines above, so CKF is scored on exactly the taxa the other methods were
## scored on. The standalone script this replaces re-parsed the study and
## rebuilt the filter itself, which is how its taxon set drifted from this one.
## NOTE. The selex_mrckf_out.RData currently on disk is NOT usable: its FP
## line was missing a "!", so FP was computed with the same expression as TP
## and every selex CKF cell read x/(x) -- 0.74/(0.74) in the summary table.
## This run overwrites it with a correctly scored one.
ckf <- run_ckf(Y_raw, X, truth = rownames(selex_truth),
               R = CKF_R, q = CKF_Q, offset = CKF_OFFSET,
               se_target = CKF_SE_TARGET, R_min = CKF_R_MIN,
               derand_thresh = CKF_DERAND, seed = 1234,
               dataset = "McMurrough2014", pair = NA_character_)
res_tpfp <- add_ckf_row(res_tpfp, ckf)
record_ckf(ckf, dir = file.path(ROOT, "selex"), file = "selex_mrckf_out.RData")
## What R would have been enough, from the fits just done.
print(ckf_stability(ckf, se_target = 0.1))

saveRDS(res_tpfp, file.path(ROOT, "selex/res_tpfp_09012026.RData"))

readRDS(file.path(ROOT, "selex/res_tpfp_09012026.RData"))
