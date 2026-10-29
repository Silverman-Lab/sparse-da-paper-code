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
# R_LIBS <- character(0)   # this script set no .libPaths(); add paths here if needed
# R_LIBS <- R_LIBS[dir.exists(R_LIBS)]
# if (length(R_LIBS)) .libPaths(c(R_LIBS, .libPaths()))

## MiMIX Julia source. NOT just ROOT + a subdir: the cluster keeps this under
## Simulation_SparseDOSSA/ while the git repo calls the same directory simulation/,
## so it needs its own knob rather than riding on ROOT.
## Pick whichever layout is actually present rather than hardcoding one: a wrong
## path here does NOT stop the script -- system() just returns non-zero and the
## MiMIX block goes on to read whatever beta.tsv was left over from an earlier
## run, silently reporting stale results. The stopifnot() makes that fail loudly.
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
library(baySeq)
library(latex2exp)
library(ggrepel)
library(ANCOMBC)
library(dirmult)
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

parse_2016_stammler_microbiome_micehuman <- function(paths = NULL) {
  local <- file.path(ROOT, "2016_stammler_microbiome_micehuman/")
  
  dat <- read_table(paste0(local, "40168_2016_175_MOESM9_ESM.txt"), skip = 1)
  counts <- as.matrix(dat[,2:37])
  rownames(counts) <- dat$OTU_ID
  
  tax <- as.matrix(dat[,38:44])
  rownames(tax) <- dat$OTU_ID
  
  scale <- readxl::read_xlsx(paste0(local, "40168_2016_175_MOESM1_ESM.xlsx"), skip = 5, col_names = FALSE, sheet = "Tabelle1")[,c(1,5)]
  colnames(scale) <- c("Sample", "16S rDNA copies per sample")
  
  idmap <- read.csv(paste0(local, "40168_2016_175_MOESM11_ESM.csv"), sep = "\t")
  
  metadata <- as.data.frame(read_table(paste0(local, "40168_2016_175_MOESM7_ESM.txt"))[-c(23,24),])
  rownames(metadata) <- metadata$SampleID
  
  samID <- idmap$SampleID[match(scale$Sample, idmap$SampleNumber)]
  scale <- scale[samID %in% rownames(metadata),]
  scale$SampleID <- samID[samID %in% rownames(metadata)]
  
  counts <- counts[,match(rownames(metadata), colnames(counts))]
  
  #remove samples that doesn't match
  counts <- counts[,colnames(counts) %in% rownames(metadata)]
  metadata <- metadata[metadata$SampleID %in% colnames(counts),]
  rownames(metadata) <- metadata$SampleID
  scale <- scale[scale$SampleID %in% colnames(counts),]
  
  all.equal(colnames(counts), rownames(metadata))
  all.equal(colnames(counts), scale$SampleID)
  all.equal(rownames(metadata), scale$SampleID)
  
  return(list(scale=scale, metadata=metadata, counts=counts, tax=tax))
}

#-------------------------------------------------------------------------------

## Reading in data------------------------------------------------

out <- parse_2016_stammler_microbiome_micehuman()
saveRDS(out, file.path(ROOT, "2016_stammler_microbiome_micehuman/output_2016_stammler_microbiome_micehuman.RData"))

scale <- out$scale$`16S rDNA copies per sample`

inds1 <- which(out$metadata$Treatment == "Pool1")
id1 <- out$metadata$SampleID[inds1]
inds6 <- which(out$metadata$Treatment == "Pool6")
id6 <- out$metadata$SampleID[inds6]

# only include sample of interest (pool 1 and pool 6)
OTU <- t(out$counts[,c(inds1,inds6)])
taxa <- out$tax
metadata = out$metadata[c(inds1,inds6),]
scale = as.numeric(out$scale$`16S rDNA copies per sample`[c(inds1,inds6)])


##Creating a phyloseq object
otu1 <- otu_table(OTU, taxa_are_rows=FALSE)
sam1 <- sample_data(metadata) 
tax1 <- tax_table(as.matrix(taxa))

phylo <- phyloseq(otu1, sam1, tax1)

##From methods: taxa unclassified at the genus level were excluded, so the
##analysis is run on genus groups.
phylo_genus <- tax_glom(phylo, taxrank = "Genus")

## FEATURE FILTER: mean raw count >= MIN_MEAN over the samples in this
## comparison, i.e. rowMeans(counts) >= 5 -- the project-wide rule (see header).
##
## This replaces the relative-abundance rule (>= 0.05% of a sample's own total,
## in >= 3 samples) the August version used. That rule normalised per sample,
## which is the better-behaved choice when depth varies; this one does not, and
## is used anyway because it is the rule the other studies apply and cross-study
## feature counts have to be produced the same way to be comparable.
##
## Computed on phylo_genus, the table actually analysed, so a genus is judged on
## its summed counts and not on its tax_glom representative OTU alone. It IS a
## filter_taxa() predicate now -- a mean needs only the taxon's own vector,
## unlike the relative-abundance form, which needed the sample totals and so had
## to be written out by hand.
##
## PREVALENCE CONSEQUENCE: rowMeans >= 5 implies no prevalence floor -- a genus
## holding all of its reads in one of 12 samples sits at prevalence 0.083 and can
## still clear it. ANCOM-BC2's prv_cut and LinDA's prev.cut are therefore 0 in
## this file, not the 0.10 / 0.1 that were safe under the old rule. The message
## prints the realised minimum so the assumption stays measured.
MIN_MEAN   <- 5
filter.tax <- function(vec) mean(vec) >= MIN_MEAN
red_phylo  <- filter_taxa(phylo_genus, filter.tax, prune = TRUE)
.Mkeep <- as(otu_table(red_phylo), "matrix")
if (!taxa_are_rows(red_phylo)) .Mkeep <- t(.Mkeep)   # taxa x samples
message(sprintf("mean count >= %g: kept %d of %d genera (min prevalence among kept %.2f)",
                MIN_MEAN, ntaxa(red_phylo), ntaxa(phylo_genus),
                min(rowSums(.Mkeep > 0)) / nsamples(red_phylo)))
rm(.Mkeep)

##Estimate of theta^\perp from the data
inds <- which(sample_data(red_phylo)$Treatment == "Pool6")

normalized_GC <- apply(t(otu_table(red_phylo))+.5, 2, FUN = function(x){x/sum(x)})
gm_GC <- apply(normalized_GC, 2, FUN = function(x){mean(log2(x))})
mean(-1*gm_GC[inds]) - mean(-1*gm_GC[-inds])

## Value implied by the dirichlet samples
## CLR estimate
clr.samps <- aldex.clr(t(otu_table(red_phylo)), sample_data(red_phylo)$Treatment, mc.samples = 1000, denom = "all")
dir.sams <- clr.samps@dirichletData

gm_imp <- rep(NA,1000)

for(i in 1:1000){
  tmp <- matrix(NA,nrow = dim(red_phylo@otu_table)[2], ncol = dim(red_phylo@otu_table)[1])
  for(k in 1:dim(red_phylo@otu_table)[1]){
    tmp[,k] <- dir.sams[[k]][,i]
  }
  tmp_gm <- apply(tmp, 2, FUN = function(x){mean(log2(x))})
  gm_imp[i] <- mean(-1*mean(tmp_gm[inds])) - mean(-1*mean(tmp_gm[-inds]))
}

mean(gm_imp)

## Estimate of theta^\perp from scale
mean(-log2(scale[inds])) - mean(-log2(scale[-inds]))

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
X <- sample_data(red_phylo)$Treatment

D = nrow(Y_raw)
N = ncol(Y_raw)
propY <- matrix(NA, nrow = D, ncol = N)
for(n in 1:N){
  propY[,n] <- ((Y_raw[,n]+1)/colSums(Y_raw)[n]) # colSums = sequencing depth
}
LFC <- rowMeans(log2(propY[,X=="Pool6"])) - rowMeans(log2(propY[,X=="Pool1"]))
hist(LFC, breaks = 100, main = paste("Treatment", unique(X)[1], "vs", unique(X)[2]))

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
gamma <-  c(1e-3,.1,.25,.5,0.7,1,1.5,2, 3,4,5,6,7,8,9,10)
gamma <- c(1e-3,.5)
for(j in 1:length(gamma)){
  scale_var <- rep(gamma[j], ncol(Y))
  
  scale_samples <- matrix(NA, nrow = ncol(Y), ncol = 1000)
  for(i in 1:ncol(Y)){
    scale_samples[i,] <- 2^rnorm(1000, scale_mean[i], scale_var[i])
  }
  clr <- aldex.clr(Y_raw,X, gamma = scale_samples,mc.samples = 1000)
  res <- aldex.ttest(clr)
  res <- cbind(res, aldex.lfc(clr))
  sen_res[[j]] <- res
}
hist(rowMeans(clr@scaleSamps))
hist(scale_mean)

dim(clr@scaleSamps)

hist(colMeans(clr@scaleSamps[inds,]) -colMeans(clr@scaleSamps[-inds,]))

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

hist(LFC, breaks = 100, main = paste("Treatment", unique(X)[1], "vs", unique(X)[2]))
abline(v = LFC[rownames(Y) %in% tax_scale2$category], col = "red")

## Other methods ---------------------------------------------------------------

## DeSeq2
coldata <- matrix(sample_data(red_phylo)$Treatment, ncol = 1)
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

X <- sample_data(red_phylo)$Treatment
conds <- sample_data(red_phylo)$Treatment 
coldata <- ifelse(conds == "Pool1", 0, 1)
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

X <- sample_data(red_phylo)$Treatment
conds <- sample_data(red_phylo)$Treatment 
coldata <- ifelse(conds == "Pool1", 0, 1)
coldata <- matrix(coldata, ncol = 1)
design <- model.matrix(~0+conds)

y <- DGEList(counts=as.matrix(Y_raw@.Data),group=conds)
y <- calcNormFactors(y, method = "TMM")
v <- voom(y,design,plot=FALSE)
fit <- lmFit(v,design)

contr <- makeContrasts(condsPool1 - condsPool6, levels = colnames(coef(fit)))
contr.fit <- contrasts.fit(fit, contr)
contr.fit <- eBayes(contr.fit)
tax_limma <- topTable(contr.fit, p.value = 1,number=100) %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)


##baySeq
conds <- sample_data(red_phylo)$Treatment 
coldata <- ifelse(conds == "Pool1", 0, 1)
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
x <-  as.numeric(X == "Pool6")
write.csv(t(Y_raw), paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX"), "/Y.csv"), row.names = F)
write.csv(x, paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX"), "/X.csv"), row.names = F)
write.csv(1:length(x), paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX/"), "/Z.csv"), row.names = F)

src_path <- MIMIX_SRC
res_path <- file.path(ROOT, "2016_stammler_microbiome_micehuman")
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

julia.time <- system.time(system(julia_command))

mimix_beta <- readr::read_tsv(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX/beta.tsv"))
mimix_out <- as.data.frame(cbind(t(apply(mimix_beta, 2, FUN = function(x){quantile(x, c(0.025, 0.975))}))))
Sig <- rep("FALSE", D)
Sig[(mimix_out$`2.5%` * mimix_out$`97.5%`) > 0] <- "TRUE"

mimix_out <- cbind(mimix_out, Sig = Sig)
rownames(mimix_out) <- rownames(Y)
tax_mimix <- mimix_out[mimix_out$Sig == "TRUE",]
tax_mimix$category <- rownames(tax_mimix)
tax_mimix


# Add ANCOMBC
ancombc.time <- system.time(ancbcfit <- ancombc(data = red_phylo, assay_name = "counts",
                                                tax_level = NULL, #phyloseq = NULL,
                                                formula = "Treatment",
                                                p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                                                group = "Treatment", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                max_iter = 100, conserve = TRUE, alpha = 0.05, #global = TRUE,
                                                n_cl = 1, verbose = TRUE))

ancbc_out <- cbind(ancbcfit$res$lfc[,c("taxon","TreatmentPool6")], W = ancbcfit$res$W$TreatmentPool6, sig = ancbcfit$res$diff_abn[,"TreatmentPool6"])

rownames(ancbc_out) <- rownames(Y)
tax_ancbc <- ancbc_out[ancbc_out$sig == "TRUE",]
tax_ancbc$category <- rownames(tax_ancbc)
tax_ancbc

# Add ANCOM
ancom.time <- system.time(ancfit <- ancom(data = red_phylo, assay_name = "counts",
                                          tax_level = NULL, #phyloseq = NULL,
                                          p_adj_method = "BH", prv_cut = 0, lib_cut = 0, #prv_cut : omits taxa with prevalence less than given value
                                          main_var = "Treatment", #adj_formula = "age + region",
                                          rand_formula = NULL, lme_control = NULL, struc_zero = FALSE,#TRUE,
                                          neg_lb = TRUE, alpha = 0.05, n_cl = 2))

anc_out <- data.frame(W = ancfit$res$W, sig = ancfit$res$detected_0.7)

rownames(anc_out) <- rownames(Y)
tax_anc <- as.data.frame(anc_out[anc_out$sig == "TRUE",])
tax_anc$category <- rownames(tax_anc)
tax_anc


## linDA
meta <- as.matrix(red_phylo@sam_data$Treatment)
colnames(meta) <- "conds"
rownames(meta) <- rownames((red_phylo@otu_table))
linda.time <- system.time(linda.obj <- linda(t(red_phylo@otu_table), meta, formula = '~conds', alpha = 0.05,
                                             prev.cut = 0, lib.cut = 0, winsor.quan = 0.97))

linda_out <- linda.obj$output$conds
tax_linda <- linda_out[linda_out$reject == "TRUE",]
tax_linda$category <- rownames(tax_linda)
tax_linda


# Add ANCOMBC2
## The `ancbc2.inds <- !colnames(...) %in% c("EU454563")` exclusion that used to
## sit here has been REMOVED (2026-08-27). Two reasons:
##
##   1. It was already a no-op. EU454563 does not survive the old mean(vec) > 5
##      filter, so `ancbc2.inds` was all TRUE and the line below subset nothing.
##      Under the prevalence filter EU454563 IS retained, at which point the line
##      became live -- and it assigned straight into the @otu_table slot without
##      touching @tax_table, leaving 48 taxa in one and 49 in the other. ANCOM-BC2
##      converts the phyloseq object to a SummarizedExperiment and that mismatch
##      fails validity: "'x@assays' is not parallel to 'x'".
##   2. Even if the mechanism were fixed (prune_taxa() keeps the slots in sync),
##      excluding a taxon HERE would leave ANCOM-BC2 testing 48 features while the
##      other eleven methods test 49 -- the same silent asymmetry that setting
##      prv_cut / lib_cut to 0 was meant to remove.
##
## ancombc2() runs to completion on the full feature set; verified 2026-08-27.
## Three genera are absent from an entire pool (EF602801, AATA01000105,
## GU272326) and get NA effect sizes with a "NaNs produced" warning, which is
## ANCOM-BC2 reporting that it cannot estimate them -- not a failure.
ancombc2.time <- system.time(ancbc2fit <- ancombc2(data = red_phylo, assay_name = "counts",
                                                   fix_formula = "Treatment", rand_formula = NULL,
                                                   p_adj_method = "BH", pseudo_sens = TRUE,
                                                   prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                                   group = "Treatment", struc_zero = FALSE, neg_lb = FALSE,
                                                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                   global = FALSE, pairwise = TRUE, 
                                                   dunnet = FALSE, trend = FALSE,
                                                   iter_control = list(tol = 1e-5, max_iter = 20, 
                                                                       verbose = FALSE),
                                                   em_control = list(tol = 1e-5, max_iter = 100),
                                                   lme_control = NULL, 
                                                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                                                   trend_control = NULL))

ancbc2_out <- as.data.frame(cbind(taxon = ancbc2fit$res[,c("taxon")], sig = ancbc2fit$res$diff_TreatmentPool6))

rownames(ancbc2_out) <- ancbc2fit$res$taxon
tax_ancbc2 <- ancbc2_out[ancbc2_out$sig == "TRUE",]
tax_ancbc2$category <- rownames(tax_ancbc2)
tax_ancbc2


##--------- sparse SSRV

## sparse SSRV takes the RAW counts (Y_raw), NOT the pseudocounted Y that the

## other methods use. Y is Y_raw + 1, so passing it here would make the Dirichlet

## concentration Y_raw + 1.5 instead of Y_raw + 0.5. Corrected 2026-08-19.

Y_0 <- Y_raw
X_0 <- as.numeric(X == "Pool6")

S <- 50000
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

record_scale_model(sparse_out, dataset = "Stammler2016",
                   dir = file.path(ROOT, "2016_stammler_microbiome_micehuman"), pair = "Pool6",
                   S = S, N = N)

## All THREE scale-model results from this one fit ----------------------------
## sparse_out is the hybrid rule's answer. sparse_auto() now also returns the
## laplace-only and bootstrap-only fits, so the run log says WHICH taxa the
## unselected scale model would have called, not just how many. Printed here and
## written to sparse_arm_features.csv by record_scale_model() above.
arm_report <- report_arm_features(sparse_out, "Stammler2016", pair = "Pool6")
tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")
tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")
## The two spiked sequences -- this study's true positives.
stammler_truth <- c("AB076660XXXX", "AB247615XXXX")
cat("TP/FP against the spike-ins --",
    " laplace:",   sum(tax_ssrv_L$category %in% stammler_truth), "/",
                   sum(!tax_ssrv_L$category %in% stammler_truth),
    " bootstrap:", sum(tax_ssrv_B$category %in% stammler_truth), "/",
                   sum(!tax_ssrv_B$category %in% stammler_truth), "\n")

sparse_sens <- check_sparse_bw_sens(dat.array = resam.array, Y_0, X_0, D, N, S,
                                    filtered.ind = rep(FALSE, D))
saveRDS(sparse_sens, file.path(ROOT, "2016_stammler_microbiome_micehuman/sparse_sens_09012026.RData"))

tax_ssrv <- sparse_out[sparse_out$sig == "TRUE",]
tax_ssrv$category <- rownames(tax_ssrv)
colnames(tax_ssrv)[4] <- "Sig"
tax_ssrv_qval <- sparse_out[sparse_out$sig.q == "TRUE",]
tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
colnames(tax_ssrv_qval)[4] <- "Sig"



## Results Summary ----------------------------------------------------------------

list_tax <- list(ALDEx2_original = tax_clr, ALDEx2_def = tax_default, 
                 Sparse_SSRV = tax_ssrv_qval,
                 DESeq2 = tax_deseq, edgeR = tax_edgeR, limma = tax_limma, baySeq = tax_baySeq, MiMIX = tax_mimix, 
                 `ANCOM-BC` = tax_ancbc, ANCOM = tax_anc, `ANCOM-BC2` = tax_ancbc2, LinDA = tax_linda)

saveRDS(list_tax, file.path(ROOT, "2016_stammler_microbiome_micehuman/results_09012026.RData"))
list_tax <- readRDS(file.path(ROOT, "2016_stammler_microbiome_micehuman/results_09012026.RData"))

## TRUTH. The two spiked OTUs, which are the only taxa known to differ between
## the pools -- the same reference the all-pairs block below and the CKF row
## above already use.
##
## This used to read `tax_gold = list_tax$ALDEx2_original`, i.e. ALDEx2's own
## CLR calls, so the table scored every method on agreement with ALDEx2 rather
## than on recovery of the spike-ins, and was not comparable with the per-pair
## tables a few hundred lines down that do use the spike-ins.
res_tpfp <- summarize_tpfp(
  tax_gold = data.frame(category = stammler_truth, stringsAsFactors = FALSE),
  list_tax)
res_tpfp


## ---- CKF -------------------------------------------------------------------
## Fitted on Y_0 / X_0, the SAME matrix and contrast handed to sparse_auto() a
## few lines above, so CKF is scored on exactly the taxa the other methods were
## scored on. The standalone script this replaces re-parsed the study and
## rebuilt the filter itself, which is how its taxon set drifted from this one.
## TRUTH SET. stammler_truth is the two spiked OTUs -- the same reference
## summarize_tpfp() is given below, so every row of this table is scored
## against the same thing.
ckf <- run_ckf(Y_0, X_0, truth = stammler_truth,
               R = CKF_R, q = CKF_Q, offset = CKF_OFFSET,
               se_target = CKF_SE_TARGET, R_min = CKF_R_MIN,
               derand_thresh = CKF_DERAND, seed = 1234,
               dataset = "Stammler2016", pair = "Pool6")
res_tpfp <- add_ckf_row(res_tpfp, ckf)
record_ckf(ckf, dir = file.path(ROOT, "2016_stammler_microbiome_micehuman"), file = "stammler_mrckf_out.RData")
## What R would have been enough, from the fits just done.
print(ckf_stability(ckf, se_target = 0.1))

saveRDS(res_tpfp, file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_09012026.RData"))
res_tpfp <- readRDS(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_09012026.RData"))


















################################################################################


##------------------------------------------------------------------------------------------
## Fully automatic
##------------------------------------------------------------------------------------------

out <- parse_2016_stammler_microbiome_micehuman()
set.seed(1234)

pool_names <- paste0("Pool", 1:6)
pairs <- combn(pool_names, 2)

for (i in 1:ncol(pairs)) {
  ref_pool <- pairs[1, i]
  target_pool <- pairs[2, i]
  
  message(paste("Running comparison:", ref_pool, "vs", target_pool))
  
  scale <- out$scale$`16S rDNA copies per sample`
  
  inds1 <- which(out$metadata$Treatment == ref_pool)
  id1 <- out$metadata$SampleID[inds1]
  inds6 <- which(out$metadata$Treatment == target_pool)
  id6 <- out$metadata$SampleID[inds6]
  
  # only include sample of interest (pool 1 and pool 6)
  OTU <- t(out$counts[,c(inds1,inds6)])
  taxa <- out$tax
  metadata = out$metadata[c(inds1,inds6),]
  scale = as.numeric(out$scale$`16S rDNA copies per sample`[c(inds1,inds6)])
  
  ##Creating a phyloseq object
  otu1 <- otu_table(OTU, taxa_are_rows=FALSE)
  sam1 <- sample_data(metadata) 
  tax1 <- tax_table(as.matrix(taxa))
  
  phylo <- phyloseq(otu1, sam1, tax1)
  
  ##Amalgamating to genus level
  phylo_genus <- tax_glom(phylo, taxrank = "Genus")
  
  ## Same rule as the main analysis above. MIN_MEAN is the global set there and
  ## is deliberately NOT re-declared here, so the all-pairs loop cannot drift
  ## onto a different filter from the headline comparison.
  ## Recomputed per pair, because the mean depends on which samples are in the
  ## subset -- 12 samples for the Pool1-vs-Pool6 headline, 2 pools' worth here.
  Mp <- as(otu_table(phylo_genus), "matrix")
  if (!taxa_are_rows(phylo_genus)) Mp <- t(Mp)       # taxa x samples
  keep_ids  <- taxa_names(phylo_genus)[rowMeans(Mp) >= MIN_MEAN]

  red_phylo <- prune_taxa(keep_ids, phylo_genus)

  ##Estimate of theta^\perp from the data
  inds <- which(sample_data(red_phylo)$Treatment == target_pool)
  
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
  X <- sample_data(red_phylo)$Treatment
  
  D = nrow(Y)
  N = ncol(Y)
  propY <- matrix(NA, nrow = D, ncol = N)
  ## Loop over `n`, NOT `i`: `i` is the pair counter of the enclosing loop, and
  ## overwriting it leaves `i` == N (10-12) at the prompt afterwards, which reads
  ## like the pair loop died partway through when in fact it ran to completion.
  for(n in 1:N){
    propY[,n] <- ((Y[,n]+1)/colSums(Y)[n]) # colSums = sequencing depth
  }
  LFC <- rowMeans(log2(propY[,X==target_pool])) - rowMeans(log2(propY[,X==ref_pool]))
  
  ## ALDEx2 models----------------------------------------------------------------
  aldex_fit <- aldex(Y_raw,X,mc.samples = 1000)
  tax_clr = aldex_fit %>% 
    rownames_to_column("category") %>%
    dplyr::select(category, effect, we.ep, we.eBH) %>%
    mutate(pval = we.ep) %>%
    mutate(padj = we.eBH) %>%
    dplyr::filter(padj < 0.05)
  
  aldex_fit <- aldex(Y_raw,X,mc.samples = 1000, gamma = .5)
  tax_default = aldex_fit %>% 
    rownames_to_column("category") %>%
    dplyr::select(category, effect, we.ep, we.eBH) %>%
    mutate(pval = we.ep) %>%
    mutate(padj = we.eBH) %>%
    dplyr::filter(padj < 0.05)
  
  ## DeSeq2
  coldata <- matrix(sample_data(red_phylo)$Treatment, ncol = 1)
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
  X <- sample_data(red_phylo)$Treatment
  conds <- sample_data(red_phylo)$Treatment 
  coldata <- ifelse(conds == ref_pool, 0, 1)
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
  X <- sample_data(red_phylo)$Treatment
  conds <- sample_data(red_phylo)$Treatment 
  design <- model.matrix(~0+conds)
  
  y <- DGEList(counts=as.matrix(Y_raw@.Data),group=conds)
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y,design,plot=FALSE)
  fit <- lmFit(v,design)
  
  contrast_str <- paste0("conds", target_pool, " - conds", ref_pool)
  contr <- makeContrasts(contrasts = contrast_str, levels = colnames(design))
  contr.fit <- contrasts.fit(fit, contr)
  contr.fit <- eBayes(contr.fit)
  tax_limma <- topTable(contr.fit, p.value = 1, number = 100) %>%
    rownames_to_column("category") %>%
    mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>%
    dplyr::filter(sig == TRUE)
  
  ##baySeq
  conds <- sample_data(red_phylo)$Treatment 
  coldata <- ifelse(conds == ref_pool, 0, 1)
  groups <- list(NDE=rep(1, length(conds)), 
                 DE=as.numeric(coldata))
  dds <- new("countData", data = as.matrix(Y_raw@.Data), groups = groups, replicates = conds)
  libsizes(dds) <- getLibsizes(dds)
  dds <- getPriors.NB(dds, cl = NULL)
  dds <- getLikelihoods(dds,cl=NULL,bootStraps=3,verbose=FALSE)
  res <- topCounts(dds,group="DE", number = 1600)
  rownames(res) <- taxa_names(red_phylo)
  
  post <- as.matrix(dds@posteriors)
  colnames(post) = c("NDE", "DEEst")
  tax_baySeq <- cbind(res, post) %>%
    rownames_to_column("category") %>%
    mutate(sig = ifelse(FWER.DE < 0.05, TRUE, FALSE)) %>%
    dplyr::filter(sig == TRUE)
  
  ### MiMIX
  x <- as.numeric(X == target_pool)
  # MINIMAL CHANGE: Added pool name to file paths to prevent overwriting during loop
  mimix_path <- paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX_"), target_pool)
  if(!dir.exists(mimix_path)) dir.create(mimix_path)
  
  write.csv(t(Y_raw), paste0(mimix_path, "/Y.csv"), row.names = F)
  write.csv(x, paste0(mimix_path, "/X.csv"), row.names = F)
  write.csv(1:length(x), paste0(mimix_path, "/Z.csv"), row.names = F)
  
  src_path <- MIMIX_SRC
  
  julia_command <- paste(
    shQuote(JULIA), shQuote(file.path(src_path, "scripts/fit-mcmc.jl")),
    "--hyper", shQuote(file.path(src_path, "nutnet-analysis/configs/hyper.yml")),
    "--monitor", shQuote(file.path(src_path, "nutnet-analysis/configs/monitor-mimix-no-factors.yml")),
    "--inits", shQuote(file.path(src_path, "nutnet-analysis/configs/inits.yml")),
    "--factors", "0",
    "--iters", "5000",
    "--chains", "1",
    shQuote(mimix_path),
    shQuote(mimix_path)
  )
  
  ## Delete the previous fit first, then insist Julia actually produced a new one.
  ## Without this, any failure of the julia call (missing binary, bad MIMIX_SRC,
  ## Julia error) leaves the old beta.tsv in place and the lines below happily
  ## report the *previous* pair's MiMIX result as if it were this pair's.
  beta_file <- file.path(mimix_path, "beta.tsv")
  if (file.exists(beta_file)) file.remove(beta_file)

  mimix_status <- system(julia_command)
  if (mimix_status != 0 || !file.exists(beta_file)) {
    stop("MiMIX/Julia failed for ", ref_pool, " vs ", target_pool,
         " (exit status ", mimix_status, "); no beta.tsv written to ", mimix_path)
  }

  mimix_beta <- readr::read_tsv(file.path(mimix_path, "beta.tsv"))
  mimix_out <- as.data.frame(cbind(t(apply(mimix_beta, 2, FUN = function(x){quantile(x, c(0.025, 0.975))}))))
  Sig <- rep("FALSE", D)
  Sig[(mimix_out$`2.5%` * mimix_out$`97.5%`) > 0] <- "TRUE"
  
  mimix_out <- cbind(mimix_out, Sig = Sig)
  rownames(mimix_out) <- rownames(Y)
  tax_mimix <- mimix_out[mimix_out$Sig == "TRUE",]
  tax_mimix$category <- rownames(tax_mimix)
  
  # Add ANCOMBC
  ancbcfit <- ancombc(data = red_phylo, assay_name = "counts",
                      tax_level = NULL, formula = "Treatment",
                      p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                      group = "Treatment", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                      max_iter = 100, conserve = TRUE, alpha = 0.05,
                      n_cl = 1, verbose = TRUE)
  
  # MINIMAL CHANGE: Dynamic column indexing for W and sig
  ancbc_out <- cbind(ancbcfit$res$lfc[,c("taxon",paste0("Treatment",target_pool))], 
                     W = ancbcfit$res$W[[paste0("Treatment", target_pool)]], 
                     sig = ancbcfit$res$diff_abn[[paste0("Treatment", target_pool)]])
  
  rownames(ancbc_out) <- rownames(Y)
  tax_ancbc <- ancbc_out[ancbc_out$sig == "TRUE",]
  tax_ancbc$category <- rownames(tax_ancbc)
  
  # Add ANCOM
  ancfit <- ancom(data = red_phylo, assay_name = "counts",
                  tax_level = NULL, p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                  main_var = "Treatment", rand_formula = NULL, lme_control = NULL, 
                  struc_zero = FALSE, neg_lb = TRUE, alpha = 0.05, n_cl = 2)
  
  anc_out <- data.frame(W = ancfit$res$W, sig = ancfit$res$detected_0.7)
  rownames(anc_out) <- rownames(Y)
  tax_anc <- as.data.frame(anc_out[anc_out$sig == "TRUE",])
  tax_anc$category <- rownames(tax_anc)
  
  ## linDA
  meta <- as.matrix(red_phylo@sam_data$Treatment)
  colnames(meta) <- "conds"
  rownames(meta) <- rownames((red_phylo@otu_table))
  linda.obj <- linda(t(red_phylo@otu_table), meta, formula = '~conds', alpha = 0.05,
                     prev.cut = 0, lib.cut = 0, winsor.quan = 0.97)
  
  linda_out <- linda.obj$output$conds
  tax_linda <- linda_out[linda_out$reject == "TRUE",]
  tax_linda$category <- rownames(tax_linda)
  
  # Add ANCOMBC2
  ## (A commented-out EU454563 exclusion and a self-assignment into @otu_table
  ## used to sit here. Both removed 2026-08-27 -- the assignment was a no-op, and
  ## the live version of it in the main analysis is what broke ANCOM-BC2 by
  ## leaving @tax_table out of step with @otu_table.)
  ancbc2fit <- ancombc2(data = red_phylo, assay_name = "counts",
                        fix_formula = "Treatment", rand_formula = NULL,
                        p_adj_method = "BH", pseudo_sens = TRUE,
                        prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                        group = "Treatment", struc_zero = FALSE, neg_lb = FALSE,
                        alpha = 0.05, n_cl = 2, verbose = TRUE,
                        global = FALSE, pairwise = TRUE, 
                        dunnet = FALSE, trend = FALSE)
  
  # MINIMAL CHANGE: Dynamic indexing for diff_Treatment
  ancbc2_out <- as.data.frame(cbind(taxon = ancbc2fit$res[,c("taxon")], 
                                    sig = ancbc2fit$res[[paste0("diff_Treatment", target_pool)]]))
  
  rownames(ancbc2_out) <- ancbc2fit$res$taxon
  tax_ancbc2 <- ancbc2_out[ancbc2_out$sig == "TRUE",]
  tax_ancbc2$category <- rownames(tax_ancbc2)
  
  ##--------- sparse SSRV
  ## sparse SSRV takes the RAW counts (Y_raw), NOT the pseudocounted Y that the
  ## other methods use. Y is Y_raw + 1, so passing it here would make the Dirichlet
  ## concentration Y_raw + 1.5 instead of Y_raw + 0.5. Corrected 2026-08-19.
  Y_0 <- Y_raw
  X_0 <- as.numeric(X == target_pool)
  S <- 10000
  D <- nrow(Y_0)
  N <- ncol(Y_0)
  
  resam.array <- array(NA, dim = c(S,D,N))
  for(n in 1:N){
    dirichlet_samples <- rdirichlet(S, as.vector(Y_0[,n]+0.5))
    resam.array[, , n] <- (log2(dirichlet_samples))
  }
  
  sparse_out <- sparse_auto(dat.array = resam.array, Y_0, X_0, D, N, S, 
                            filtered.ind = rep(FALSE, D))
  
  record_scale_model(sparse_out, dataset = "Stammler2016",
                     dir = file.path(ROOT, "2016_stammler_microbiome_micehuman"), pair = paste0(ref_pool, "_vs_", target_pool),
                     S = S, N = N)

  ## Both scale-model arms for THIS pair, keyed the same way the record is, so
  ## all_pairs_results.R can put a bootstrap and a laplace panel next to the
  ## selected one. See the headline comparison above for the full note.
  arm_report <- report_arm_features(sparse_out, "Stammler2016",
                                    pair = paste0(ref_pool, "_vs_", target_pool))
  tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")
  tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")

  sparse_sens <- check_sparse_bw_sens(dat.array = resam.array, Y_0, X_0, D, N, S, file_path = NULL,
                                      filtered.ind = rep(FALSE, D))
  saveRDS(sparse_sens, paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/sparse_sens_"),ref_pool,"_",target_pool,"_09012026.RData"))
  
  tax_ssrv_qval <- sparse_out[sparse_out$sig.q == "TRUE",]
  tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
  colnames(tax_ssrv_qval)[4] <- "Sig"
  
  ## Results Summary ----------------------------------------------------------------
  tax_truth <- list(category = c("AB076660XXXX", "AB247615XXXX"))
  
  list_tax <- list(ALDEx2_original = tax_clr, ALDEx2_def = tax_default, 
                   sparse_SSRV = tax_ssrv_qval, 
                   DESeq2 = tax_deseq, edgeR = tax_edgeR, limma = tax_limma, baySeq = tax_baySeq, MiMIX = tax_mimix, 
                   `ANCOM-BC` = tax_ancbc, ANCOM = tax_anc, `ANCOM-BC2` = tax_ancbc2, LinDA = tax_linda)
  
  saveRDS(list_tax, paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/results_"),ref_pool,"_",target_pool,"_09012026.RData"))
  
  res_tpfp <- summarize_tpfp(tax_gold = tax_truth, list_tax)

  ## ---- CKF -------------------------------------------------------------------
  ## Fitted on Y_0 / X_0, the SAME matrix and contrast handed to sparse_auto() a
  ## few lines above, so CKF is scored on exactly the taxa the other methods were
  ## scored on. The standalone script this replaces re-parsed the study and
  ## rebuilt the filter itself, which is how its taxon set drifted from this one.
  ## ALL 15 POOL PAIRS. legacy_index = i writes this pair into slot i of the
  ## list-of-matrices shape stammler_mrckf_out_all_pairs.RDS has always had, so
  ## the existing colMeans(ckf_stammler[[5]]) line keeps working while the
  ## keyed record in ckf_record.rds carries every pair by name.
  ckf <- run_ckf(Y_0, X_0, truth = tax_truth$category,
                 R = CKF_R, q = CKF_Q, offset = CKF_OFFSET,
                 se_target = CKF_SE_TARGET, R_min = CKF_R_MIN,
                 derand_thresh = CKF_DERAND, seed = 1234,
                 dataset = "Stammler2016", pair = paste0(ref_pool, "_vs_", target_pool))
  res_tpfp <- add_ckf_row(res_tpfp, ckf)
  record_ckf(ckf, dir = file.path(ROOT, "2016_stammler_microbiome_micehuman"), file = "stammler_mrckf_out_all_pairs.RDS", legacy_index = i)

  saveRDS(res_tpfp, paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_"),ref_pool,"_",target_pool,"_09012026.RData"))
}


# {
#   pool <- "Pool2"
#   scale <- out$scale$`16S rDNA copies per sample`
#   
#   inds1 <- which(out$metadata$Treatment == "Pool1")
#   id1 <- out$metadata$SampleID[inds1]
#   inds6 <- which(out$metadata$Treatment == pool)
#   id6 <- out$metadata$SampleID[inds6]
#   
#   # only include sample of interest (pool 1 and pool 6)
#   OTU <- t(out$counts[,c(inds1,inds6)])
#   taxa <- out$tax
#   metadata = out$metadata[c(inds1,inds6),]
#   scale = as.numeric(out$scale$`16S rDNA copies per sample`[c(inds1,inds6)])
#   
#   
#   ##Creating a phyloseq object
#   otu1 <- otu_table(OTU, taxa_are_rows=FALSE)
#   sam1 <- sample_data(metadata) 
#   tax1 <- tax_table(as.matrix(taxa))
#   
#   phylo <- phyloseq(otu1, sam1, tax1)
#   
#   ##From methods: "taxa unclassified at the genus level or present in less than 20% of samples were excluded from statistical analysis
#   ##Amalgamating to genus level
#   phylo_genus <- tax_glom(phylo, taxrank = "Genus")
#   
#   # #Filtering OTU table in phyloseq
#   # filter.tax <- function(vec){
#   #  return(sum(vec > 1) > .1*length(vec))
#   # }
#   filter.tax <- function(vec){
#     return(mean(vec) > 5)
#   }
#   red_phylo <- filter_taxa(phylo_genus, filter.tax, prune = TRUE)
#   
#   ##Estimate of theta^\perp from the data
#   inds <- which(sample_data(red_phylo)$Treatment == pool)
#   
#   ## Count matrices for this comparison, defined ONCE.
#   ##   Y_raw : counts as they come out of phyloseq. THIS is what every model is
#   ##           fitted on -- see count_input_policy.R for method-by-method reasons.
#   ##   Y     : Y_raw + 1, kept only for derived quantities that want a zero-free
#   ##           matrix (log-ratio summaries, prevalence filters written as Y > 1).
#   ##           No model is fitted on it.
#   ## Until 2026-08-19 every count-based method was fed Y. On the Jin data that
#   ## alone was the difference between reproducing the source paper's DESeq2 counts
#   ## (8 proximal / 1 distal) and not (9 / 5).
#   Y_raw <- t(otu_table(red_phylo))
#   Y     <- Y_raw + 1
#   X <- sample_data(red_phylo)$Treatment
#   
#   D = nrow(Y)
#   N = ncol(Y)
#   propY <- matrix(NA, nrow = D, ncol = N)
#   ## Loop over `n`, NOT `i`: `i` is the pair counter of the enclosing loop, and
#   ## overwriting it leaves `i` == N (10-12) at the prompt afterwards, which reads
#   ## like the pair loop died partway through when in fact it ran to completion.
#   for(n in 1:N){
#     propY[,n] <- ((Y[,n]+1)/colSums(Y)[n]) # colSums = sequencing depth
#   }
#   LFC <- rowMeans(log2(propY[,X==pool])) - rowMeans(log2(propY[,X=="Pool1"]))
#   
#   #-------------------------------------------------------------------------------
#   
#   ## ALDEx2 models----------------------------------------------------------------
#   ## First aldex2
#   aldex_fit <- aldex(Y_raw,X,mc.samples = 1000)
#   
#   tax_clr = aldex_fit %>% 
#     rownames_to_column("category") %>%
#     dplyr::select(category, effect, we.ep, we.eBH) %>%
#     mutate(pval = we.ep) %>%
#     mutate(padj = we.eBH) %>%
#     dplyr::filter(padj < 0.05)
#   
#   ##Default scale model
#   aldex_fit <- aldex(Y_raw,X,mc.samples = 1000, gamma = .5)
#   
#   tax_default = aldex_fit %>% 
#     rownames_to_column("category") %>%
#     dplyr::select(category, effect, we.ep, we.eBH) %>%
#     mutate(pval = we.ep) %>%
#     mutate(padj = we.eBH) %>%
#     dplyr::filter(padj < 0.05)
#   
#   ## Other methods ---------------------------------------------------------------
#   
#   ## DeSeq2
#   coldata <- matrix(sample_data(red_phylo)$Treatment, ncol = 1)
#   colnames(coldata) <- "Condition"
#   ## Raw counts, with the size-factor fallback described in
#   ## count_input_policy.R -- the default median-of-ratios estimator cannot
#   ## cope when every feature has a zero somewhere, which is what the old
#   ## +1 was masking.
#   dds <- deseq_dataset_raw(Y_raw@.Data,
#                            colData = coldata,
#                            design  = ~Condition)
#   dds <- DESeq(dds)
#   res <- deseq_results(dds)
#   tax_deseq <- res %>% 
#     as.data.frame() %>% 
#     rownames_to_column("category") %>% 
#     dplyr::select(category, log2FoldChange, padj, lfcSE) %>% 
#     mutate(sig = ifelse(padj < 0.05, TRUE, FALSE)) %>%
#     dplyr::filter(sig == TRUE)
#   
#   
#   ##EdgeR
#   
#   X <- sample_data(red_phylo)$Treatment
#   conds <- sample_data(red_phylo)$Treatment 
#   coldata <- ifelse(conds == "Pool1", 0, 1)
#   coldata <- matrix(coldata, ncol = 1)
#   
#   y <- DGEList(counts=as.matrix(Y_raw@.Data),group=c(coldata[,1]))
#   y <- calcNormFactors(y)
#   y <- estimateDisp(y)
#   et <- exactTest(y)
#   tax_edgeR <- topTags(et, n=1600, p.value = 1) %>%
#     as.data.frame() %>%
#     rownames_to_column("category") %>%
#     mutate(sig = ifelse(FDR < 0.05, TRUE, FALSE)) %>%
#     dplyr::filter(sig == TRUE)
#   
#   
#   ##limma
#   
#   X <- sample_data(red_phylo)$Treatment
#   conds <- sample_data(red_phylo)$Treatment 
#   coldata <- ifelse(conds == "Pool1", 0, 1)
#   coldata <- matrix(coldata, ncol = 1)
#   design <- model.matrix(~0+conds)
#   
#   y <- DGEList(counts=as.matrix(Y_raw@.Data),group=conds)
#   y <- calcNormFactors(y, method = "TMM")
#   v <- voom(y,design,plot=FALSE)
#   fit <- lmFit(v,design)
#   
#   contr <- makeContrasts(condsPool1 - condsPool2, levels = colnames(coef(fit)))
#   contr.fit <- contrasts.fit(fit, contr)
#   contr.fit <- eBayes(contr.fit)
#   tax_limma <- topTable(contr.fit, p.value = 1,number=100) %>%
#     rownames_to_column("category") %>%
#     mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>%
#     dplyr::filter(sig == TRUE)
#   
#   
#   ##baySeq
#   conds <- sample_data(red_phylo)$Treatment 
#   coldata <- ifelse(conds == "Pool1", 0, 1)
#   groups <- list(NDE=rep(1, length(conds)), 
#                  DE=as.numeric(coldata))
#   dds <- new("countData", data = as.matrix(Y_raw@.Data), groups = groups, replicates = conds)
#   libsizes(dds) <- getLibsizes(dds)
#   dds <- getPriors.NB(dds, cl = NULL)
#   dds <- getLikelihoods(dds,cl=NULL,bootStraps=3,verbose=FALSE)
#   res <- topCounts(dds,group="DE", number = 1600)
#   rownames(res) <- taxa_names(red_phylo)
#   
#   ##Merge to the data to see which taxa because the names are uninformative
#   post <- as.matrix(dds@posteriors)
#   colnames(post) = c("NDE", "DEEst")
#   tax_baySeq <- cbind(res, post) %>%
#     rownames_to_column("category") %>%
#     mutate(sig = ifelse(FWER.DE < 0.05, TRUE, FALSE)) %>%
#     dplyr::filter(sig == TRUE)
#   
#   ### MiMIX
#   x <-  as.numeric(X == pool)
#   write.csv(t(Y_raw), paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX"), "/Y.csv"), row.names = F)
#   write.csv(x, paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX"), "/X.csv"), row.names = F)
#   write.csv(1:length(x), paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX/"), "/Z.csv"), row.names = F)
#   
#   src_path <- MIMIX_SRC
#   res_path <- file.path(ROOT, "2016_stammler_microbiome_micehuman")
#   #dir.create(file.path(src_path, "SparseDOSSA-analysis-results/08042024/", paste("results",safe_param_string,i, sep = "_")),recursive=TRUE,showWarnings=TRUE) # change date accordingly
#   
#   julia_command <- paste(
#     shQuote(JULIA), shQuote(file.path(src_path, "scripts/fit-mcmc.jl")),
#     "--hyper", shQuote(file.path(src_path, "nutnet-analysis/configs/hyper.yml")),
#     "--monitor", shQuote(file.path(src_path, "nutnet-analysis/configs/monitor-mimix-no-factors.yml")),
#     "--inits", shQuote(file.path(src_path, "nutnet-analysis/configs/inits.yml")),
#     "--factors", "0",
#     "--iters", "5000",
#     "--chains", "1",
#     shQuote(file.path(res_path, "MiMIX")),
#     shQuote(file.path(res_path, "MiMIX"))
#   )
#   
#   julia.time <- system.time(system(julia_command))
#   
#   mimix_beta <- readr::read_tsv(file.path(ROOT, "2016_stammler_microbiome_micehuman/MiMIX/beta.tsv"))
#   mimix_out <- as.data.frame(cbind(t(apply(mimix_beta, 2, FUN = function(x){quantile(x, c(0.025, 0.975))}))))
#   Sig <- rep("FALSE", D)
#   Sig[(mimix_out$`2.5%` * mimix_out$`97.5%`) > 0] <- "TRUE"
#   
#   mimix_out <- cbind(mimix_out, Sig = Sig)
#   rownames(mimix_out) <- rownames(Y)
#   tax_mimix <- mimix_out[mimix_out$Sig == "TRUE",]
#   tax_mimix$category <- rownames(tax_mimix)
#   tax_mimix
#   
#   
#   # Add ANCOMBC
#   ancombc.time <- system.time(ancbcfit <- ancombc(data = red_phylo, assay_name = "counts",
#                                                   tax_level = NULL, #phyloseq = NULL,
#                                                   formula = "Treatment",
#                                                   p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
#                                                   group = "Treatment", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
#                                                   max_iter = 100, conserve = TRUE, alpha = 0.05, #global = TRUE,
#                                                   n_cl = 1, verbose = TRUE))
#   
#   ancbc_out <- cbind(ancbcfit$res$lfc[,c("taxon",paste0("Treatment",pool))], W = ancbcfit$res$W$TreatmentPool2, sig = ancbcfit$res$diff_abn[,paste0("Treatment",pool)])
#   
#   rownames(ancbc_out) <- rownames(Y)
#   tax_ancbc <- ancbc_out[ancbc_out$sig == "TRUE",]
#   tax_ancbc$category <- rownames(tax_ancbc)
#   tax_ancbc
#   
#   # Add ANCOM
#   ancom.time <- system.time(ancfit <- ancom(data = red_phylo, assay_name = "counts",
#                                             tax_level = NULL, #phyloseq = NULL,
#                                             p_adj_method = "BH", prv_cut = 0, lib_cut = 0, #prv_cut : omits taxa with prevalence less than given value
#                                             main_var = "Treatment", #adj_formula = "age + region",
#                                             rand_formula = NULL, lme_control = NULL, struc_zero = FALSE,#TRUE,
#                                             neg_lb = TRUE, alpha = 0.05, n_cl = 2))
#   
#   anc_out <- data.frame(W = ancfit$res$W, sig = ancfit$res$detected_0.7)
#   
#   rownames(anc_out) <- rownames(Y)
#   tax_anc <- as.data.frame(anc_out[anc_out$sig == "TRUE",])
#   tax_anc$category <- rownames(tax_anc)
#   tax_anc
#   
#   
#   ## linDA
#   meta <- as.matrix(red_phylo@sam_data$Treatment)
#   colnames(meta) <- "conds"
#   rownames(meta) <- rownames((red_phylo@otu_table))
#   linda.time <- system.time(linda.obj <- linda(t(red_phylo@otu_table), meta, formula = '~conds', alpha = 0.05,
#                                                prev.cut = 0, lib.cut = 0, winsor.quan = 0.97))
#   
#   linda_out <- linda.obj$output$conds
#   tax_linda <- linda_out[linda_out$reject == "TRUE",]
#   tax_linda$category <- rownames(tax_linda)
#   tax_linda
#   
#   
#   # Add ANCOMBC2
#   ancbc2.inds <- !colnames(red_phylo@otu_table) %in% c("EU454563")
#   red_phylo@otu_table <- red_phylo@otu_table[,ancbc2.inds]
#   ancombc2.time <- system.time(ancbc2fit <- ancombc2(data = red_phylo, assay_name = "counts",
#                                                      fix_formula = "Treatment", rand_formula = NULL,
#                                                      p_adj_method = "BH", pseudo_sens = TRUE,
#                                                      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
#                                                      group = "Treatment", struc_zero = FALSE, neg_lb = FALSE,
#                                                      alpha = 0.05, n_cl = 2, verbose = TRUE,
#                                                      global = FALSE, pairwise = TRUE, 
#                                                      dunnet = FALSE, trend = FALSE,
#                                                      iter_control = list(tol = 1e-5, max_iter = 20, 
#                                                                          verbose = FALSE),
#                                                      em_control = list(tol = 1e-5, max_iter = 100),
#                                                      lme_control = NULL, 
#                                                      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
#                                                      trend_control = NULL))
#   
#   ancbc2_out <- as.data.frame(cbind(taxon = ancbc2fit$res[,c("taxon")], sig = ancbc2fit$res$diff_TreatmentPool2))
#   
#   rownames(ancbc2_out) <- ancbc2fit$res$taxon
#   tax_ancbc2 <- ancbc2_out[ancbc2_out$sig == "TRUE",]
#   tax_ancbc2$category <- rownames(tax_ancbc2)
#   tax_ancbc2
#   
#   
#   ##--------- sparse SSRV
#   
#   ## sparse SSRV takes the RAW counts (Y_raw), NOT the pseudocounted Y that the
#   
#   ## other methods use. Y is Y_raw + 1, so passing it here would make the Dirichlet
#   
#   ## concentration Y_raw + 1.5 instead of Y_raw + 0.5. Corrected 2026-08-19.
#   
#   Y_0 <- Y_raw
#   X_0 <- as.numeric(X == pool)
#   
#   S <- 50000
#   D <- nrow(Y_0)
#   N <- ncol(Y_0)
#   
#   resam.array <- array(NA, dim = c(S,D,N))
#   sparse_time1 <- system.time(
#     for(n in 1:N){
#       dirichlet_samples <- rdirichlet(S, as.vector(Y_0[,n]+0.5))
#       resam.array[, , n] <- (log2(dirichlet_samples))
#     }
#   )
#   
#   sparse_time <- system.time(
#     sparse_out <- sparse_auto(dat.array = resam.array, Y_0, X_0, D, N, S, 
#                               filtered.ind = rep(FALSE, D))
#   )
#   
#   record_scale_model(sparse_out, dataset = "Stammler2016",
#                      dir = file.path(ROOT, "2016_stammler_microbiome_micehuman"), pair = paste0("Pool1_vs_", pool),
#                      S = S, N = N)
# 
#   ## Both scale-model arms for this Pool1 comparison, keyed as the record is.
#   arm_report <- report_arm_features(sparse_out, "Stammler2016",
#                                     pair = paste0("Pool1_vs_", pool))
#   tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")
#   tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")
# 
#   sparse_sens <- check_sparse_bw_sens(dat.array = resam.array, Y_0, X_0, D, N, S, file_path = NULL,
#                                       filtered.ind = rep(FALSE, D))
#   saveRDS(sparse_sens, paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/sparse_sens_Pool1_"),pool,"_09012026.RData"))
#   
#   tax_ssrv_qval <- sparse_out[sparse_out$sig.q == "TRUE",]
#   tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
#   colnames(tax_ssrv_qval)[4] <- "Sig"
#   tax_ssrv_qval
#   
#   
#   ## Results Summary ----------------------------------------------------------------
#   tax_truth <- c()
#   tax_truth$category <- c("AB076660XXXX", "AB247615XXXX")
#   
#   list_tax <- list(ALDEx2_original = tax_clr, ALDEx2_def = tax_default, 
#                    sparse_SSRV = tax_ssrv_qval,
#                    DESeq2 = tax_deseq, edgeR = tax_edgeR, limma = tax_limma, baySeq = tax_baySeq, MiMIX = tax_mimix, 
#                    `ANCOM-BC` = tax_ancbc, ANCOM = tax_anc, `ANCOM-BC2` = tax_ancbc2, LinDA = tax_linda)
#   
#   saveRDS(list_tax, paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/results_Pool1_"),pool,"_09012026.RData"))
#   list_tax <- readRDS(paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/results_Pool1_"),pool,"_09012026.RData"))
#   
#   res_tpfp <- summarize_tpfp(tax_gold = tax_truth, list_tax)
#   res_tpfp
#   
# 
#   ## ---- CKF -------------------------------------------------------------------
#   ## Fitted on Y_0 / X_0, the SAME matrix and contrast handed to sparse_auto() a
#   ## few lines above, so CKF is scored on exactly the taxa the other methods were
#   ## scored on. The standalone script this replaces re-parsed the study and
#   ## rebuilt the filter itself, which is how its taxon set drifted from this one.
#   ckf <- run_ckf(Y_0, X_0, truth = tax_truth$category,
#                  R = CKF_R, q = CKF_Q, offset = CKF_OFFSET,
#                  se_target = CKF_SE_TARGET, R_min = CKF_R_MIN,
#                  derand_thresh = CKF_DERAND, seed = 1234,
#                  dataset = "Stammler2016", pair = paste0("Pool1_vs_", pool))
#   res_tpfp <- add_ckf_row(res_tpfp, ckf)
#   ## legacy_index follows stammler_mrckf_out_filter5.RDS's own convention:
#   ## the standalone script stored Pool2..Pool6 at slots 1..5, so "Pool2" is 1.
#   record_ckf(ckf, dir = file.path(ROOT, "2016_stammler_microbiome_micehuman"),
#              file = "stammler_mrckf_out_filter5.RDS",
#              legacy_index = as.integer(sub("Pool", "", pool)) - 1L)
# 
#   saveRDS(res_tpfp, paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_Pool1_"),pool,"_09012026.RData"))
#   res_tpfp <- readRDS(paste0(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_Pool1_"),pool,"_09012026.RData"))
# }
# 
# 


## These are the Pool1-vs-PoolX comparisons, read under the "<prefix>_<comp>"
## short naming. Those short-name files are LEFTOVERS from an older naming
## scheme -- they still exist on disk but were last written 12 Aug, when the
## MiMIX/Julia call was silently failing. The loop above writes the same
## comparisons as "<prefix>_Pool1_<comp>", so read those instead; otherwise
## these two plots are built from stale results while everything else is current.
res_tpfp2 <- readRDS(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_Pool1_Pool2_09012026.RData"))
res_tpfp3 <- readRDS(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_Pool1_Pool3_09012026.RData"))
res_tpfp4 <- readRDS(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_Pool1_Pool4_09012026.RData"))
res_tpfp5 <- readRDS(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_Pool1_Pool5_09012026.RData"))
res_tpfp6 <- readRDS(file.path(ROOT, "2016_stammler_microbiome_micehuman/res_tpfp_Pool1_Pool6_09012026.RData"))

ckf_tpfp2_6 <- readRDS(file.path(ROOT, "2016_stammler_microbiome_micehuman/stammler_mrckf_out_filter5.RDS"))
res_ckf_2_6 <- lapply(ckf_tpfp2_6, function(x) colMeans(x))

## Dropped: identical(ckf_tpfp2_6, ckf_tpfp2_61) -- ckf_tpfp2_61 is never
## defined anywhere, so this scratch comparison halted the script here, before
## either of the two plots below was drawn.

res_2_6 <- list(res_tpfp2$aldex_gold, res_tpfp3$aldex_gold, res_tpfp4$aldex_gold, res_tpfp5$aldex_gold, res_tpfp6$aldex_gold)
res_2_6 <- base::mapply(rbind, res_2_6, CKF = res_ckf_2_6, SIMPLIFY = FALSE)
#colnames(res_2_6) <- paste0(rep(c("TP","FP"), 5), 2:6)

res_2_6 <- lapply(res_2_6, function(x) {
  rownames(x)[3] <- "Sparse_SSRV" 
  return(x)
})

library(dplyr)
combined_data <- bind_rows(
  lapply(seq_along(res_2_6), function(i) {
    as.data.frame(res_2_6[[i]]) %>%
      mutate(Method = rownames(res_2_6[[1]]), Pool = paste0("Pool ", (i+1)))
  })
)

# Make sure TP is numeric for plotting
combined_data$TP <- as.numeric(combined_data$TP)
combined_data$TP <- as.numeric(combined_data$TP)

combined_data$TDR <- combined_data$TP/2
combined_data$FDR <- combined_data$FP / (combined_data$TP + combined_data$FP)
combined_data$FDR[is.na(combined_data$FDR)] <- 0
combined_data$type1 <- combined_data$FP / 40 # 40 ; number of nulls


# Plot the data using ggplot2
library(ggplot2)

grey_color = "grey80"

colors <- c(
  "Sparse_SSRV" = "red", 
  #"ALDEx2_clr" = "skyblue1",           
  "ALDEx2_def" = "skyblue3",           
  "ANCOM" = "plum1", #"#A96A3E",            
  "ANCOM-BC" =  "plum1", #"#A96A3E",#"#e8b086",         
  "ANCOM-BC2" =  "plum3", #"#A96A3E",#"#894B33FF", 
  "LinDA" = "olivedrab3", 
  "MiMIX" = "yellow2", #"#fd9696",              
  "DESeq2" = grey_color,  # Set to grey  
  "baySeq" = grey_color,  # Set to grey
  "edgeR" = grey_color,  # Set to grey
  "limma" = grey_color,   # Set to grey
  "CKF" = grey_color
)

# Define distinct line types for grey methods
line_types <- c(
  #"ALDEx2_clr" = "dashed",
  "ALDEx2_def" = "dotdash",
  "ANCOM" = "dotdash",
  "ANCOM-BC" = "dotted",
  "ANCOM-BC2" = "longdash",
  "DESeq2" = "dashed",  # Unique line type for grey method
  "LinDA" = "dotdash",
  "MiMIX" = "dashed",
  "CKF" = "solid",
  "baySeq" = "dotted",  # Unique line type for grey method
  "edgeR" = "dotdash",  # Unique line type for grey method
  "limma" = "longdash",   # Unique line type for grey method
  "Sparse_SSRV" = "solid"
)

# ggplot(combined_data, aes(x = Pool, y = TP, group = Method, color = Method)) +
#   geom_line() +
#   labs(
#     title = "True Positives (TP) Across Methods",
#     x = "Matrix",
#     y = "True Positives (TP)",
#     color = "Method"
#   ) +
#   scale_color_manual(values = colors) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Dropped: a legend_override built from `method_colors`. That object is never
## defined (the palette here is called `colors`), so this halted the script --
## and legend_override was never referenced afterwards anyway; both plots below
## get their legend from scale_color_manual(colors) / scale_linetype_manual().

neworder <- c("Sparse_SSRV", "ALDEx2_original", "ALDEx2_def", "ANCOM", "ANCOM-BC", "ANCOM-BC2", "LinDA", "MiMIX", "DESeq2", "baySeq", "edgeR", "limma", "CKF")

combined_data$Method <- factor(combined_data$Method, levels = neworder)

library(ggplot2)
ggplot(combined_data, aes(x = Pool, y = TDR, group = Method, color = Method, linetype = Method)) +
  geom_line(size = 1.1, position =position_dodge(0.2)) +
  labs(
    title = "TDR across pools",
    x = "Matrix",
    y = "TDR",
    color = "Method"
  ) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(ROOT, "2016_stammler_microbiome_micehuman/tdr_dodge_09012026.pdf"), width = 8, height = 6)


ggplot(combined_data, aes(x = Pool, y = FDR, group = Method, color = Method, linetype = Method)) +
  geom_line(size = 1.1, position =position_dodge(0.2)) +
  labs(
    title = "FDR across pools",
    x = "Matrix",
    y = "FDR",
    color = "Method"
  ) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  geom_hline(yintercept = 0.05, color = "grey", linetype = "dashed", size = 0.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(ROOT, "2016_stammler_microbiome_micehuman/fdr_dodge_09012026.pdf"), width = 8, height = 6)


# ggplot(combined_data, aes(x = Pool, y = type1, group = Method, color = Method, linetype = Method)) +
#   geom_line(size = 1.1, position =position_dodge(0.2)) +
#   ylim(0, 0.15) +
#   labs(
#     title = "Type 1 Error Rate across pools",
#     x = "Matrix",
#     y = "Type 1 Error Rate",
#     color = "Method"
#   ) +
#   scale_color_manual(values = colors) +
#   scale_linetype_manual(values = line_types) +
#   geom_hline(yintercept = 0.05, color = "grey", linetype = "dashed", size = 0.8) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(file.path(ROOT, "2016_stammler_microbiome_micehuman/type1error_dodge.pdf"), width = 8, height = 6)
