################################################################################
## FILTER RULE: UNCHANGED -- this study is one of the two exceptions.
##                                                     (this run: 09012026)
##
## The project reverted to rowMeans(raw counts) >= 5 on 2026-09-01 for Prasse,
## Stammler, Zemb and selex. Vandeputte keeps its 20% prevalence rule, sum(vec > 1) > .2*length(vec),
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
library(driver)#not working yet
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

#-------------------------------------------------------------------------------

## Reading in the Vandputte data------------------------------------------------

### Metadata has disease status as well as avg cell count per data
metadata = read.csv(file.path(ROOT, "Vandeputte", "data","41586_2017_BFnature24460_MOESM10_ESM.csv"))

###Basic manipulation of the metadata
metadata = metadata %>%
  dplyr::rename("CellCount" = "Average.cell.count..per.gram.of.frozen.feces.")

###Quick plot of the cell counts by status
ggplot(metadata, aes(x = Health.status, y=CellCount)) + 
  geom_boxplot() +
  xlab("Health status") +
  ylab("Faecal cell count per gram") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

###Reading in the OTU table
OTU = fread(file.path(ROOT, "Vandeputte", "data","OTU_nochim.csv"))[-c(1:40),]

taxa = fread(file.path(ROOT, "Vandeputte", "data","taxa_assignments.csv"), header = TRUE) %>%
  column_to_rownames("V1")
taxa_id <- paste0("sp", 1:nrow(taxa))

rownames(taxa) = taxa_id
taxa <- as.matrix(taxa)
OTU <- as.data.frame(OTU) %>%
  column_to_rownames("V1")
colnames(OTU) <- taxa_id
sam_id <- paste0("sa", 1:nrow(OTU))
rownames(OTU) <- sam_id


##Creating a phyloseq object
otu1 <- otu_table(OTU, taxa_are_rows=FALSE)
sam1 <- sample_data(metadata) 
tax1 <- tax_table(taxa)

phylo <- phyloseq(otu1, sam1, tax1)

##From methods: "taxa unclassified at the genus level or present in less than 20% of samples were excluded from statistical analysis
##Amalgamating to genus level
phylo_genus <- tax_glom(phylo, taxrank = "Genus")

##Filtering OTU table in phyloseq
filter.tax <- function(vec){
  return(sum(vec > 1) > .2*length(vec))
}
red_phylo <- filter_taxa(phylo_genus, filter.tax, prune = TRUE)

##Estimate of theta^\perp from the data
inds <- which(sample_data(red_phylo)$Health.status == "CD")
normalized_GC <- apply(t(otu_table(red_phylo))+.5, 2, FUN = function(x){x/sum(x)})
gm_GC <- apply(normalized_GC, 2, FUN = function(x){mean(log2(x))})
mean(-1*gm_GC[inds]) - mean(-1*gm_GC[-inds])

## Value implied by the dirichlet samples
## CLR estimate
clr.samps <- aldex.clr(t(otu_table(red_phylo)), sample_data(red_phylo)$Health.status, mc.samples = 1000, denom = "all")
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



## Estimate of theta^\perp due to the flow cytometry measurements
mean(log2(sample_data(red_phylo)$CellCount[inds])) - mean(log2(sample_data(red_phylo)$CellCount[-inds]))

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
Y_raw <- t(otu_table(red_phylo))
Y     <- Y_raw + 1
X <- sample_data(red_phylo)$Health.status
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
scale_mean <- log2(sample_data(red_phylo)$CellCount)
gamma <-  c(1e-3,.1,.25,.5,0.7,1,1.5,2, 3,4,5,6,7,8,9,10)
gamma <- c(1e-3,.7)
for(j in 1:length(gamma)){
  scale_var <- rep(gamma[j], 95)
  
  scale_samples <- matrix(NA, nrow = 95, ncol = 1000)
  for(i in 1:95){
    scale_samples[i,] <- 2^rnorm(1000, scale_mean[i], scale_var[i])
  }
  clr <- aldex.clr(Y_raw,X, gamma = scale_samples,mc.samples = 1000)
  res <- aldex.ttest(clr)
  res <- cbind(res, aldex.lfc(clr))
  sen_res[[j]] <- res
}

## Choosing scale based off of the estimated technical variation from the other methods. See "Data" for the calculation in the excel file "41586_2017_BFnature24460_MOESM11_ESM.xlsx"
tax_scale = sen_res[[2]] %>% 
  rownames_to_column("category") %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)


## Informed model
scale.cd <- 2^matrix(rnorm(1000*29, mean = log2(.7), sd = .125), nrow = 29)
scale.control <- 2^matrix(rnorm(1000*66, mean = log2(1), sd = .125), nrow = 66)

scale.inf <- rbind(scale.cd, scale.control)
aldex_informed <- aldex(Y_raw,X,mc.samples = 1000, gamma = scale.inf)

tax_informed = aldex_informed %>% 
  rownames_to_column("category") %>%
  dplyr::select(category, effect, we.ep, we.eBH) %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)


#-------------------------------------------------------------------------------

## Other methods ---------------------------------------------------------------

## DeSeq2
coldata <- matrix(sample_data(red_phylo)$Health.status, ncol = 1)
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

X <- sample_data(red_phylo)$Health.status
conds <- sample_data(red_phylo)$Health.status 
coldata <- ifelse(conds == "CD", 0, 1)
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
X <- sample_data(red_phylo)$Health.status
conds <- sample_data(red_phylo)$Health.status 
coldata <- ifelse(conds == "CD", 0, 1)
coldata <- matrix(coldata, ncol = 1)
design <- model.matrix(~0+conds)

y <- DGEList(counts=as.matrix(Y_raw@.Data),group=conds)
y <- calcNormFactors(y, method = "TMM")
v <- voom(y,design,plot=FALSE)
fit <- lmFit(v,design)

contr <- makeContrasts(condsCD - condsControl, levels = colnames(coef(fit)))
contr.fit <- contrasts.fit(fit, contr)
contr.fit <- eBayes(contr.fit)
tax_limma <- topTable(contr.fit, p.value = 1,number=100) %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)


##baySeq
conds <- sample_data(red_phylo)$Health.status 
coldata <- ifelse(conds == "CD", 0, 1)
groups <- list(NDE=rep(1, length(conds)), 
               DE=as.numeric(coldata))
dds <- new("countData", data = as.matrix(Y_raw@.Data), groups = groups, replicates = conds)
libsizes(dds) <- getLibsizes(dds)
dds <- getPriors.NB(dds, cl = NULL)
dds <- getLikelihoods(dds,cl=NULL,bootStraps=1000,verbose=FALSE)
res <- topCounts(dds,group="DE", number = 1600)
rownames(res) <- taxa_names(red_phylo)

##Merge to the data to see which taxa because the names are uninformative
post <- as.matrix(dds@posteriors)
colnames(post) = c("NDE", "DEEst")
tax_baySeq <- cbind(res, post) %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(FWER.DE < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)

#-------------------------------------------------------------------------------

## vandpuette's qmp method------------------------------------------------------

# function from https://github.com/raeslab/QMP/blob/master/QMP.R
# with cnv_corrected_abundance_table: a copy number variation corrected abundance table with sample-identifiers as rows, copy number corrected taxa-abundances as columns
# with cell_counts_table: a table with sample-identifiers as rows, cell counts as columns 
rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) 
{
  try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE) stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names, Please check!"))
  cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
  cell_counts_table = t(cell_counts_table[row.names(cnv_corrected_abundance_table),]) # make sure the order of the samples is the same in both files    
  sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
  sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
  minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
  rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
  cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
  rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq), ncol = ncol(cnv_corrected_abundance_table_phyloseq), dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq), colnames(cnv_corrected_abundance_table_phyloseq)))
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
  {
    x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
    rarefied_matrix[i,] = x
  }
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,]
  return(QMP)
}
cellCounts <- data.frame(sample_data(red_phylo)) %>%
  dplyr::select(CellCount)
QMP_df <- rarefy_even_sampling_depth(otu_table(red_phylo),cellCounts)

## running a wilcoxon test on the qmp data
p.qmp <- rep(NA, ncol(QMP_df))
lfc.qmp <- rep(NA, ncol(QMP_df))
conds <- sample_data(red_phylo)$Health.status
for(i in 1:length(p.qmp)){
  x <- log2(QMP_df[which(conds == "CD"),i]+0.5)
  y <- log2(QMP_df[which(conds == "Control"),i]+0.5)
  lfc.qmp[i] <- mean(x) - mean(y)
  tmp.test <- t.test(x,y)
  p.qmp[i] <- tmp.test$p.value
}
padj.qmp <- p.adjust(p.qmp, method = "BH")

tax_QMP <- data.frame("category" = colnames(otu_table(red_phylo)),
                      "pval" = p.qmp,
                      "padj" = padj.qmp) %>%
  mutate(sig = ifelse(padj.qmp <= 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)
#-------------------------------------------------------------------------------

## Sparse model ----------------------------------------------------------------

## Y is already the +1 matrix defined above; do not re-derive it raw here.
X <- as.numeric(sample_data(red_phylo)$Health.status == "CD")

S <- 50000
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

sparse_time <- system.time(
  sparse_auto_out <- sparse_auto(dat.array = resam.array, Y_raw, X, D, N, S, #bandwidth = "nrd0"
                                 filtered.ind = rep(FALSE, D))
)

record_scale_model(sparse_auto_out, dataset = "Vandeputte2017",
                   dir = file.path(ROOT, "Vandeputte"),
                   S = S, N = N)

## Estimated theta^perp, scored against the flow cytometry -------------------
## theta^perp is the shift that turns this method's relative log fold changes
## into absolute ones; under the model it IS the log2 fold change of the total
## load between CD and healthy. Vandeputte measured that total by flow cytometry
## (cells per gram of frozen faeces), and the same difference is already computed
## by hand further up under "Estimate of theta^perp due to the flow cytometry
## measurements" -- this recomputes it against X so the two vectors are
## guaranteed to be in the same order and the same direction, and then records it
## next to the estimate instead of printing it and dropping it.
##
## X is 1 for CD, so this is treatment minus control, matching theta^perp.
cc <- sample_data(red_phylo)$CellCount
names(cc) <- sample_names(red_phylo)
stopifnot(length(cc) == N)
van_perp_obs <- mean(log2(cc[X == 1])) - mean(log2(cc[X == 0]))

## ...and the same measurement PROJECTED ONTO THE FITTED FEATURES, which is the
## quantity theta^perp actually estimates. sparse_auto() sees proportions
## renormalised within the retained genera, so log2 W_dn = log2 p^S_dn + log2
## T_S,n and the shift the mode identifies is Delta log2 T_S. The cytometry
## counts the whole community. With r_n = T_S,n / T_n,
## theta^perp_S = Delta log2 T_n + Delta log2 r_n, exactly.
##
## Y_all is `phylo` -- BEFORE tax_glom(), not phylo_genus. tax_glom() discards
## every OTU with no genus assignment (2041 OTUs -> 195 genera here), and those
## are real cells the cytometer counted, so they belong in the denominator. The
## consequence to state wherever obs_S is quoted: the matched reference now
## depends on how well the annotation database covers this community.
##
## Y_kept, NOT keep = taxa_names(red_phylo). tax_glom() names each genus after
## its most abundant constituent OTU, so those ids all resolve against `phylo`
## and would select 29 OTUs where the model was fitted on 29 GENERA -- a wrong
## retained fraction that no error would report.
##
## Caveat specific to this study: r_n is a READ fraction and the assay counts
## CELLS, so 16S copy-number variation and extraction bias enter the projection.
## Vandeputte is the weakest of the four studies on that footing, and it is also
## the one where the correction makes the bias slightly worse rather than better.
van_perp_matched <- perp_matched_obs(
  total = cc, X = X,
  Y_all = as(otu_table(phylo), "matrix"),
  taxa_are_rows = taxa_are_rows(phylo),
  Y_kept = Y_raw, kept_taxa_are_rows = TRUE)

perp_report <- record_theta_perp(sparse_auto_out, dataset = "Vandeputte2017",
                                 dir = file.path(ROOT, "Vandeputte"),
                                 S = S, N = N, obs = van_perp_obs,
                                 matched = van_perp_matched,
                                 notes = "flow cytometry: cells per gram of frozen faeces")

## All THREE scale-model results from this one fit ----------------------------
## sparse_auto_out is the hybrid rule's answer. sparse_auto() now also returns the
## laplace-only and bootstrap-only fits, so the run log says WHICH taxa the
## unselected scale model would have called, not just how many. Printed here and
## written to sparse_arm_features.csv by record_scale_model() above.
arm_report <- report_arm_features(sparse_auto_out, "Vandeputte2017")
tax_ssrv_L <- arm_sig_tax(arm_report, "laplace")
tax_ssrv_B <- arm_sig_tax(arm_report, "bootstrap")

sparse_sens <- check_sparse_bw_sens(dat.array = resam.array, Y_raw, X, D, N, S, bandwidth = c("nrd","SJ-ste"),
                                    filtered.ind = rep(FALSE, D))
saveRDS(sparse_sens, file.path(ROOT, "Vandeputte/sparse_sens_09012026.RData"))
#readRDS(file.path(ROOT, "Vandeputte/sparse_sens_09012026.RData"))

tax_ssrv_qval <- sparse_auto_out[sparse_auto_out$sig.q == "TRUE",]
tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
colnames(tax_ssrv_qval)[4] <- "Sig"

tax_ssrv <- sparse_auto_out[sparse_auto_out$sig == "TRUE",]
tax_ssrv$category <- rownames(tax_ssrv)
colnames(tax_ssrv)[4] <- "Sig"


### MiMIX
src_path <- MIMIX_SRC
res_path <- file.path(ROOT, "Vandeputte")

## Write this run's design/counts rather than relying on whatever X/Y/Z.csv
## happened to be sitting in data/ from a previous session. Y is D x N here, and
## MiMIX wants samples x taxa, hence the transpose.
mimix_in <- file.path(res_path, "data")
dir.create(mimix_in, showWarnings = FALSE, recursive = TRUE)
write.csv(t(Y_raw), file.path(mimix_in, "Y.csv"), row.names = FALSE)
write.csv(X, file.path(mimix_in, "X.csv"), row.names = FALSE)
write.csv(1:length(X), file.path(mimix_in, "Z.csv"), row.names = FALSE)
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
Sig <- rep("FALSE", D)
Sig[(mimix_out$`2.5%` * mimix_out$`97.5%`) > 0] <- "TRUE"

mimix_out <- cbind(mimix_out, Sig = Sig)
rownames(mimix_out) <- rownames(Y)
mimix_out <- cbind(mimix_out, taxa = taxa[,"Genus"][match(rownames(mimix_out), rownames(taxa))])
tax_mimix <- mimix_out[mimix_out$Sig == "TRUE",]
tax_mimix$category <- rownames(tax_mimix)
tax_mimix


library(LinDA)
meta <- as.matrix(red_phylo@sam_data$Health.status)
colnames(meta) <- "Health.status"
linda.time <- system.time(linda.obj <- linda(t(red_phylo@otu_table), meta, formula = '~Health.status', alpha = 0.05,
                                             prev.cut = 0, lib.cut = 0, winsor.quan = 0.97))

linda_out <- linda.obj$output$Health.statusControl

rownames(linda_out) <- rownames(Y)
linda_out <- cbind(linda_out, taxa = taxa[,"Genus"][match(rownames(linda_out), rownames(taxa))])
tax_linda <- linda_out[linda_out$reject == "TRUE",]
tax_linda$category <- rownames(tax_linda)
tax_linda

# Add ANCOMBC
ancombc.time <- system.time(ancbcfit <- ancombc(data = red_phylo, assay_name = "counts",
                                                tax_level = NULL, #phyloseq = NULL,
                                                formula = "Health.status",
                                                p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                                                group = "Health.status", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                max_iter = 100, conserve = TRUE, alpha = 0.05, #global = TRUE,
                                                n_cl = 1, verbose = TRUE))

ancbc_out <- cbind(ancbcfit$res$lfc[,c("taxon","Health.statusControl")], W = ancbcfit$res$W$Health.statusControl, adj.p = ancbcfit$res$q_val$Health.statusControl, sig = ancbcfit$res$diff_abn[,"Health.statusControl"])

rownames(ancbc_out) <- rownames(Y)
ancbc_out <- cbind(ancbc_out, taxa = taxa[,"Genus"][match(rownames(ancbc_out), rownames(taxa))])
tax_ancbc <- ancbc_out[ancbc_out$sig == "TRUE",]
tax_ancbc$category <- rownames(tax_ancbc)
tax_ancbc

# Add ANCOM
ancom.time <- system.time(ancfit <- ancom(data = red_phylo, assay_name = "counts",
                                          tax_level = NULL, #phyloseq = NULL,
                                          p_adj_method = "BH", prv_cut = 0, lib_cut = 0, #prv_cut : omits taxa with prevalence less than given value
                                          main_var = "Health.status", #adj_formula = "age + region",
                                          rand_formula = NULL, lme_control = NULL, struc_zero = FALSE,#TRUE,
                                          neg_lb = TRUE, alpha = 0.05, n_cl = 2))

#write.csv(cbind(ancfit$res, effect_size = effect_sizes, newW), paste0(output_dir,"/ANCOM.out N",n,"D",d,".csv"))
anc_out <- data.frame(sig = ancfit$res$detected_0.7)

rownames(anc_out) <- rownames(Y)
anc_out <- cbind(anc_out, taxa = taxa[,"Genus"][match(rownames(anc_out), rownames(taxa))])
tax_anc <- anc_out[anc_out$sig == "TRUE",]
tax_anc$category <- rownames(tax_anc)
tax_anc

# Add ANCOMBC2
ancombc2.time <- system.time(ancbc2fit <- ancombc2(data = red_phylo, assay_name = "counts",
                                                   fix_formula = "Health.status", rand_formula = NULL,
                                                   p_adj_method = "BH", pseudo_sens = TRUE,
                                                   prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                                   group = "Health.status", struc_zero = FALSE, neg_lb = FALSE,
                                                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                   global = FALSE, pairwise = TRUE, 
                                                   dunnet = FALSE, trend = FALSE,
                                                   iter_control = list(tol = 1e-5, max_iter = 20, 
                                                                       verbose = FALSE),
                                                   em_control = list(tol = 1e-5, max_iter = 100),
                                                   lme_control = NULL, 
                                                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                                                   trend_control = NULL))

ancbc2_out <- as.data.frame(cbind(taxon = ancbc2fit$res[,c("taxon")], sig = ancbc2fit$res$diff_Health.statusControl))

rownames(ancbc2_out) <- ancbc2fit$res$taxon
tax_ancbc2 <- ancbc2_out[ancbc2_out$sig == "TRUE",]
tax_ancbc2$category <- rownames(tax_ancbc2)
tax_ancbc2

# exe.time <- rbind(sparse_time1, sparse_time2, sparse_time3,
#                   aldex.time, deseq2.time, ancom.time, ancombc.time, linda.time)
# write.csv(exe.time, paste0(output_dir,"/execution time.csv"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


list_tax <- list(ALDEx2_gold = tax_scale, ALDEx2_informed = tax_informed, ALDEx2_original = tax_clr, ALDEx2_def = tax_default, 
                 sparse_SSRV = tax_ssrv_qval,
                 DESeq2 = tax_deseq, 
                 edgeR = tax_edgeR, limma = tax_limma, baySeq = tax_baySeq, MiMIX = tax_mimix, 
                 `ANCOM-BC` = tax_ancbc, ANCOM = tax_anc, `ANCOM-BC2` = tax_ancbc2, LinDA = tax_linda)

saveRDS(list_tax, file.path(ROOT, "Vandeputte/results_09012026.RData"))
list_tax <- readRDS(file.path(ROOT, "Vandeputte/results_09012026.RData"))

res_tpfp <- summarize_tpfp(tax_gold = tax_scale, list_tax)
res_tpfp


## ---- CKF -------------------------------------------------------------------
## Fitted on Y_raw / X, the SAME matrix and contrast handed to sparse_auto() a
## few lines above, so CKF is scored on exactly the taxa the other methods were
## scored on. The standalone script this replaces re-parsed the study and
## rebuilt the filter itself, which is how its taxon set drifted from this one.
## TRUTH SET. tax_scale is the flow-cytometry-calibrated ALDEx2 gold this script
## just built, and is the same object summarize_tpfp() was handed above -- so
## here the CKF row and every other row of the table are on the same reference.
## The standalone script read a results_02252025.RData instead, and wrote its
## output only into ~/Research/datarepo, which is why this study's CKF result
## looked missing from the project tree for months.
ckf <- run_ckf(Y_raw, X, truth = tax_scale$category,
               R = CKF_R, q = CKF_Q, offset = CKF_OFFSET,
               se_target = CKF_SE_TARGET, R_min = CKF_R_MIN,
               derand_thresh = CKF_DERAND, seed = 1234,
               dataset = "Vandeputte2017", pair = NA_character_)
res_tpfp <- add_ckf_row(res_tpfp, ckf)
record_ckf(ckf, dir = file.path(ROOT, "Vandeputte"), file = "vandeputte_mrckf_out.RData")
## What R would have been enough, from the fits just done.
print(ckf_stability(ckf, se_target = 0.1))

saveRDS(res_tpfp, file.path(ROOT, "Vandeputte/res_tpfp_09012026.RData"))
#res_tpfp <- readRDS(file.path(ROOT, "Vandeputte/res_tpfp_09012026.RData"))


