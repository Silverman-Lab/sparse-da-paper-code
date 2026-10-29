################################################################################
## count_input_policy.R -- which matrix each differential-abundance method is
## fed, and why. Sourced by all six real-data scripts so the decision lives in
## one place instead of being re-made in each of them.
##
## THE POLICY
##
##   Y_raw is the default. A method is handed the +1-pseudocounted matrix Y only
##   if it cannot run on raw counts.
##
## As of 2026-08-19 that second category is EMPTY: every method used in this
## project either takes raw counts natively or applies its own zero handling
## internally. The scripts still build Y, because a few derived quantities want
## it (log-ratio summaries, and Zemb's prevalence filter, which is written as
## `Y > 1` and therefore means `raw > 0`), but no model is fitted on it.
##
## WHY THIS CHANGED
##
## Every count-based method used to receive Y. On the Jin cOTU data that single
## choice was the whole difference between reproducing the source paper and not:
## DESeq2 called 8 proximal / 1 distal on raw counts, exactly as published, and
## 9 / 5 on Y. The distal number moving 1 -> 5 matters more than it looks,
## because the paper's claim IS the proximal/distal asymmetry; the pseudocount
## flattened it into the pattern the paper attributes to gene-abundance data.
## edgeR and limma moved too, less far and not always in the same direction.
##
## WHAT EACH METHOD ACTUALLY WANTS
##
##   DESeq2   raw. Negative-binomial GLM on counts, with median-of-ratios size
##            factors estimated from them. Adding 1 to every cell perturbs the
##            size factors and erases the structural zeros the dispersion fit
##            depends on. See deseq_size_factors() below for the one wrinkle.
##   edgeR    raw. TMM handles zeros.
##   limma    raw. voom adds its own 0.5 offset internally.
##   baySeq   raw. NB model with its own library-size handling.
##   ALDEx2   raw. aldex.clr() draws Dirichlet Monte Carlo instances with a 0.5
##            prior, which IS the zero handling; feeding it counts+1 silently
##            changes that prior to 1.5. This applies to the scale-model arms
##            too -- aldex.clr(gamma = scale_samples) is the gold standard in
##            four of the six studies, so it moved when this changed.
##   LinDA    raw. Has its own pseudo.cnt argument, default 0.5. Five of the six
##            scripts already passed it raw counts; Prasse passed Y, which was
##            simply an inconsistency between scripts.
##   ANCOM-BC / ANCOM-BC2 / ANCOM
##            raw. They are handed the phyloseq object, so they were never
##            affected either way.
##   MiMIX    raw. Already written t(Y_raw) to disk.
##   sparse_SSRV
##            raw. Already Y_0 <- Y_raw, with its own +0.5 Dirichlet prior.
################################################################################

## DESeq2's default size-factor estimator takes the median of the ratios to the
## per-feature geometric mean, and a feature with a zero in ANY sample has a
## geometric mean of zero. If that is true of every feature -- routine in sparse
## microbiome tables -- estimateSizeFactors() stops with "every gene contains at
## least one zero". Handing it counts+1 makes the error go away, which is very
## likely why the pseudocount was there in the first place, but it fixes the
## error by changing the data for every method at once.
##
## The supported fix is the "poscounts" estimator, which computes each feature's
## geometric mean over its nonzero samples only. This applies it ONLY when the
## default genuinely fails, so datasets where the default works are bit-for-bit
## unchanged, and it says so out loud when it kicks in rather than switching
## estimators silently.
deseq_size_factors <- function(dds, label = NULL) {
  ok <- tryCatch({
    dds <- DESeq2::estimateSizeFactors(dds)
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    message("  [DESeq2] default size factors failed on raw counts",
            if (is.null(label)) "" else paste0(" (", label, ")"),
            " -- every feature has at least one zero. Using type=\"poscounts\".")
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  }
  dds
}

## Build the DESeqDataSet from RAW counts and resolve its size factors up front,
## so DESeq() further down finds them already set and does not re-estimate with
## the default. Returns the object; the caller still runs DESeq() with whatever
## arguments that study uses (Jin passes betaPrior = TRUE to match its source
## paper; the others do not).
deseq_dataset_raw <- function(counts, colData, design, label = NULL) {
  counts <- as.matrix(counts)
  ## DESeq2 requires integers and refuses non-integer input. phyloseq otu_tables
  ## come back as doubles even when every value is a whole number, so round
  ## rather than let DESeqDataSetFromMatrix() do it with a warning -- and assert
  ## the values really were whole, so a genuinely non-integer table (a
  ## normalised matrix reaching this function by mistake) fails loudly instead
  ## of being quietly truncated.
  stopifnot(all(abs(counts - round(counts)) < 1e-8, na.rm = TRUE))
  storage.mode(counts) <- "integer"

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData   = colData,
                                        design    = design)
  deseq_size_factors(dds, label = label)
}

## DESeq2's results() applies INDEPENDENT FILTERING by default: it picks a
## baseMean quantile that maximises the number of rejections at alpha, and sets
## padj to NA for every feature below it. Those features are not dropped from the
## table -- they come back with a p-value and no adjusted p-value -- so a script
## that counts `padj < 0.05` silently scores DESeq2 on a SUBSET of the feature
## set every other method was given, with no error and no message.
##
## It was left on until 2026-09-01 on the grounds that it was measured to filter
## 0 features (the blocks pass counts with a +1 offset in the simulation, and the
## real-data tables clear the auto threshold), and that DESeq2's recommended
## usage keeps it on. Both are still true. It is off anyway, because "0 features
## on the runs we happened to measure" is a weaker guarantee than "the twelve
## methods are scored on the same features by construction", and this project's
## whole comparison rests on the second one.
##
## THE COST is real and worth stating: independent filtering exists to raise
## power, and turning it off costs DESeq2 some rejections it would otherwise
## make. That is the intended trade -- a like-for-like comparison, not the
## strongest possible DESeq2.
##
## `...` passes through whatever the study needs (Jin's betaPrior fit and the
## Prasse subestuary contrast both call results() with extra arguments).
deseq_results <- function(dds, ...) {
  DESeq2::results(dds, ..., independentFiltering = FALSE)
}

################################################################################
## SECOND POLICY: no method applies its own filter on top of the script's.
##
## Four of the twelve methods take filtering arguments of their own:
##
##   ANCOM-BC, ANCOM, ANCOM-BC2   prv_cut  (features, by prevalence)
##                                lib_cut  (samples, by library size)
##   LinDA                        prev.cut, lib.cut  (same two)
##
## As of 2026-09-01 ALL FOUR arguments are 0 in ALL SIX study scripts. The
## feature set and the sample set are decided ONCE, by the script, and every
## method sees the same data. There is no longer a study where a non-zero
## argument is carried on the grounds that it happens to be a no-op.
##
## WHY THAT CHANGED, and it is not tidying. Four of the six studies reverted on
## 2026-09-01 to rowMeans(raw counts) >= 5 (Prasse, Stammler, Zemb, selex;
## Vandeputte and Jin keep their papers' own rules). That rule is an ABUNDANCE
## rule and implies NO PREVALENCE FLOOR AT ALL -- a genus holding every one of
## its reads in a single sample of 12 sits at prevalence 0.083 and still clears a
## mean of 5. Under the previous prevalence and relative-abundance rules the
## 0.10 / 0.1 arguments were provable no-ops; under this rule they are not, and
## they would have quietly handed ANCOM-BC2 and LinDA a smaller feature set than
## the other ten methods saw.
##
## prv_cut / prev.cut
##   Both packages compare with >= (LinDA: `rowSums(Y > 0)/n >= prev.cut`;
##   ANCOM-BC: `which(prevalence >= prv_cut)`). The minimum prevalence each
##   script's own filter GUARANTEES, which is what a non-zero argument would have
##   to clear to be a no-op:
##
##     Prasse       rowMeans(raw) >= 5                          -> NONE (>= 1/57)
##     Stammler     rowMeans(raw) >= 5                          -> NONE (>= 1/12)
##     selex        rowMeans(raw) >= 5                          -> NONE (>= 1/14)
##     Zemb         rowMeans(raw) >= 5                          -> NONE (>= 1/8)
##     Vandeputte   20% prevalence (its paper's rule)           -> >= 0.20
##     Jin          >= MIN_COUNT reads in >= MIN_SAMPLES of 8   -> >= 0.375
##
##   Only the last two rows still guarantee anything, and both are 0 anyway so
##   the guarantee is not load-bearing. Each mean-filtered script PRINTS the
##   prevalence its survivors actually reached, so the realised floor is measured
##   per run rather than argued from the rule.
##
##   If any script ever returns to a prevalence-shaped filter, these arguments
##   may go back to 0.10 -- but there is no reason to: 0 is correct under every
##   rule, and a non-zero value is only ever correct by coincidence.
##
## lib_cut / lib.cut
##   These were 1000 in the Prasse scripts until 2026-09-01 and are the more
##   dangerous of the two, because they drop SAMPLES and say nothing when they
##   do. Measured minimum library size under the August filters:
##
##     Prasse (raw)      1127        Prasse (QSeq)   2.3e6
##     Stammler         14959        selex          127372
##     Zemb             10748        Vandeputte       2226
##
##   So 1000 was a no-op in every study -- but Prasse cleared it by only 13%, and
##   that margin is a function of the taxon filter: at prev >= 0.20 the minimum
##   was 1027, and at prev >= 0.50 five samples fell under 1000. Those five would
##   have vanished from ANCOM-BC, ANCOM, ANCOM-BC2 and LinDA while the other
##   eight methods kept them -- a 9% change to the design, affecting a third of
##   the table, with no error and no message.
##
##   The mean filter changes which taxa are kept and therefore these totals, so
##   the numbers above are NOT current. They are left as the record of how thin
##   the margin was, which is the argument for 0 rather than for a smaller
##   cutoff. 0 is the no-filter value for both APIs (`lib_size >= lib_cut`,
##   `colSums(otu.tab) >= lib.cut`) and is already the package default for the
##   three ANCOM-BC functions; LinDA's default is 1.
##
## If a genuinely degenerate sample ever needs excluding, exclude it in the
## script, where it is visible to all twelve methods and to the reader.
################################################################################
