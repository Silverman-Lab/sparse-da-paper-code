# run_all_methods() -- the 12-method differential-abundance block, lifted from
# src/05052025/sparse_analysis_trycatch.R:146-402 so it can be run more than once per
# replicate (the oracle study runs it on three data arms) without copy-paste divergence.
#
# Deliberate differences from the original, all documented in the plan:
#
#  * Truth is keyed by taxon NAME, not by position. The original built
#    `effect_sizes <- rep(0, D-1)` and indexed positionally while sim_sparse_auto()
#    built a length-D version that included the aggregated "other" row -- so Sparse
#    SSRV's type-1 error denominator counted "other" as a true null and no other
#    method's did. Name-keying makes that class of bug impossible.
#  * ALDEx2's `diff.btw` (the actual log2 CLR fold-change) is retained. summary_aldex2()
#    keeps only the standardized `effect`, which is not an LFC and cannot be compared
#    against the spiked effect size.
#  * `offsets` allows the true scale to be injected, for the scale-informed arm.
#  * library(stringr) is loaded -- the ANCOMBC2 retry handler calls str_match() and the
#    original never loaded it, so a zero-variance error silently lost ANCOMBC2.
#
# The naive path (offsets = NULL) is otherwise faithful to the original.

library(stringr)

################################################################################
## COUNT INPUT POLICY -- no method applies its own filter on top of the script's.
##
## The mirror of the SECOND POLICY block in count_input_policy.R, which states
## this for the real-data analysis. The feature set and the sample set are
## decided ONCE, by apply_filter() in simulate_one.R -- FILTER_TYPE = "mean",
## rowMeans(counts) >= 5 as of 2026-09-01 -- and every method sees the same Y.
##
## Four of the thirteen methods take filtering arguments of their own. All four
## are set to 0, the no-filter value for both APIs (`lib_size >= lib_cut`,
## `colSums(otu.tab) >= lib.cut`) and already the package default for the three
## ANCOM-BC functions:
##
##   ANCOM-BC, ANCOM, ANCOM-BC2   prv_cut, lib_cut
##   LinDA                        prev.cut, lib.cut
##
## WHY, measured 2026-08-29 rather than assumed.
##
## lib_cut / lib.cut -- these were 1000 and they were NOT no-ops. They drop
##   SAMPLES, silently, with no error and no message. Minimum library size after
##   the taxon filter, 5 seeds per cell:
##
##     lambda    min depth   samples dropped by lib_cut = 1000
##     5e3          595      3.8 / 100  (max 6, group-imbalanced by up to 3)
##     1e4         1190      0.4 / 100  (max 1)
##     5e4         5949      0
##
##   So at the leftmost point of the overLambda sweep, four methods fit on ~96
##   samples while the other nine fit on 100. That is a 4% change to the design
##   affecting a third of the table, and it is exactly the failure mode
##   count_input_policy.R:142-158 flags on the real data.
##
## prv_cut / prev.cut -- LinDA was 0.1 and ANCOM-BC2 was 0.10. Both are 0, and
##   they have now been unsafe under BOTH filter rules, for opposite reasons:
##
##     under rowMeans >= 5 (the rule again, since 2026-09-01) a mean threshold
##       implies no prevalence floor at all, so a retained taxon can sit below
##       0.10 and would be cut by those arguments but not by the script. The
##       0.20 worst-cell minimum measured before 2026-08-29 was a property of
##       that grid, not a guarantee of the rule.
##     under prevalence >= 0.10 the script cut at exactly the same 0.10 and taxa
##       landed on the boundary -- measured minimum prevalence among retained
##       taxa was 0.1000 at lambda = 5e3 and at D = 1000, so floating-point
##       equality would have decided them.
##
##   0 removes the question under either rule, which is why it is not worth
##   re-deriving a "safe" non-zero value the next time the filter moves.
##
## The other nine methods were audited and filter nothing:
##
##   DESeq2     run_deseq2() now calls results(dds, independentFiltering = FALSE)
##              -- turned OFF on 2026-09-01, at the same time as the real-data
##              scripts (deseq_results() in count_input_policy.R).
##
##              Independent filtering does not drop rows; it sets padj to NA
##              below a fitted baseMean quantile, so a block counting padj <
##              alpha scores DESeq2 on fewer features than the other twelve saw,
##              with no error and no message. It was measured at 0 features
##              filtered at lambda 5e5 and 5e3 (auto threshold baseMean 2.4-2.9,
##              identical discovery counts with it off) under FILTER_TYPE =
##              "prevalence", and that holds a fortiori under "mean", which is
##              harsher on low-abundance taxa. So this is expected to change
##              nothing -- it is off so the guarantee is structural rather than
##              measured on the cells someone happened to check.
##
##              It costs DESeq2 whatever power the filtering was buying. That is
##              the intended trade: a like-for-like comparison, not the strongest
##              possible DESeq2.
##   edgeR      no filterByExpr call.
##   ALDEx2, limma-voom, baySeq, MRCKF, MIMIX, Sparse SSRV
##              no filter arguments at all.
##
## If a genuinely degenerate sample ever needs excluding, exclude it in
## simulate_one(), where it is visible to all thirteen methods and to the reader.
################################################################################

# mrckf_fit() lives in src/common/mrckf_fit.R -- consumers must source that first.

# Y                      taxa x samples count matrix, rownames = taxa, last row "other"
# x                      n x 1 matrix of 0/1 group labels
# output_dir             per-arm directory; filenames match the original layout
# filtered_metadata_df   feature_spiked / effect_size for taxa surviving the filter
# n, d                   sample size and feature count, used only in output filenames
# offsets                NULL for the naive run, or a list for the scale-informed run:
#                          size_factor     geometric-mean-1 scale factor (DESeq2)
#                          lib_size        same, scaled to count magnitude (edgeR etc)
#                          theta_perp_true scalar, true log2 load difference (treat-cont)
# methods                which methods to run; the informed arm only has offset APIs for
#                        a subset, so it passes the shorter vector
run_all_methods <- function(Y, x, output_dir, filtered_metadata_df, n, d,
                            S = 2000, pval = 0.05, offsets = NULL,
                            methods = c("sparse", "aldex", "deseq2", "ancombc", "ancom",
                                        "linda", "edger", "limma", "bayseq", "ancombc2"),
                            mrckf_reseed = FALSE, keep_raw = FALSE,
                            # Supplementary variant grid: the NON-default KDE
                            # bandwidths to refit Sparse SSRV at. c() disables it.
                            ssrv_bandwidths = c("nrd", "SJ-ste"),
                            # passed straight through to sim_sparse_auto(); Study B
                            # forces the arm on the absolute-abundance data, Study A
                            # always uses "hybrid" and reads the three tables instead
                            scale_model = c("hybrid", "laplace", "bootstrap")) {

  scale_model <- match.arg(scale_model)

  informed <- !is.null(offsets)
  D <- nrow(Y)
  N <- ncol(Y)
  cont.ind <- x == 0
  treat.ind <- x == 1

  # --- name-keyed truth -------------------------------------------------------
  truth <- setNames(rep(0, D), rownames(Y))
  if (!is.null(filtered_metadata_df) && nrow(filtered_metadata_df) > 0) {
    hit <- as.character(filtered_metadata_df$feature_spiked)
    hit <- hit[hit %in% names(truth)]
    truth[hit] <- filtered_metadata_df$effect_size[
      match(hit, as.character(filtered_metadata_df$feature_spiked))]
  }
  # "other" is an aggregate of filtered-out taxa, not a real feature; every method's
  # metrics exclude it.
  keep <- rownames(Y) != "other"
  effect_sizes <- unname(truth[keep])

  timings <- list()
  results <- list()
  failures <- list()
  raw <- list()

  # keep_raw stashes each method's NATIVE result table -- every column, before any
  # extraction -- so calibrate_lfc_scales.R can show what each method actually returns
  # rather than only what we picked out of it. Off by default: the production path must
  # be byte-identical to before.
  stash_raw <- function(name, obj) if (isTRUE(keep_raw)) raw[[name]] <<- obj

  # store one method's per-taxon result in a common shape
  stash <- function(name, taxon, est, sig, pvalue = NA, qvalue = NA) {
    results[[name]] <<- data.frame(
      taxon    = as.character(taxon),
      est      = as.numeric(est),
      sig      = as.logical(sig),
      pval     = as.numeric(pvalue),
      qval     = as.numeric(qvalue),
      true_lfc_loge = unname(truth[as.character(taxon)]),
      stringsAsFactors = FALSE
    )
  }

  ##### Sparse SSRV ############################################################
  if ("sparse" %in% methods) {
    resam.array <- array(NA, dim = c(S, D, N))
    timings$sparse_resample <- system.time(
      for (j in 1:N) {
        dirichlet_samples <- rdirichlet(S, Y[, j] + 0.5)
        resam.array[, , j] <- (log2(dirichlet_samples))
      }
    )

    if (informed) {
      timings$sparse <- system.time(
        sparse_out <- sim_sparse_oracle(dat.array = resam.array, Y, x, D, N, S, output_dir,
                                        load = offsets$load,
                                        theta_perp_true = offsets$theta_perp_true,
                                        scale_uncertainty = offsets$scale_uncertainty,
                                        metadata_df = filtered_metadata_df,
                                        filtered.ind = rep(FALSE, D))
      )
      sparse_file <- file.path(output_dir, "sparse_oracle_filterFALSE.csv")
    } else {
      timings$sparse <- system.time(
        invisible(sim_sparse_auto(dat.array = resam.array, Y, x, D, N, S, output_dir,
                                  metadata_df = filtered_metadata_df,
                                  filtered.ind = rep(FALSE, D),
                                  scale_model = scale_model))
      )
      sparse_file <- file.path(output_dir, "sparse_auto_filterFALSE.csv")
    }

    if (file.exists(sparse_file)) {
      so <- read.csv(sparse_file, row.names = 1)
      so <- so[rownames(so) != "other", , drop = FALSE]
      stash_raw(if (informed) "Sparse_SSRV_oracle" else "Sparse_SSRV", so)
      stash(if (informed) "Sparse_SSRV_oracle" else "Sparse_SSRV", rownames(so), so$mean, so$sig.q)
    }

    ## simulation_rerun EDIT (c): stash the Laplace-only and bootstrap-only arms
    ## as two additional METHODS.
    ##
    ## sim_sparse_auto() computes both arms on every call and now writes both
    ## tables (edit (a)), so this costs NO extra model fits -- the three curves in
    ## goal 2 are three views of ONE fit, which is exactly what makes them
    ## comparable. Everything downstream is automatic: tidy_replicate() and
    ## summarise_replicate() treat them as ordinary methods and produce
    ## FDR / TDR / F0.5 for each. shared.R adds their METHOD_SCALE entries.
    ##
    ## Only on the naive path. The informed/oracle variant has a single arm by
    ## construction, so there is nothing to decompose.
    if (!informed) {
      arm_files <- c(Sparse_SSRV_L = "sparse_L_filterFALSE.csv",
                     Sparse_SSRV_B = "sparse_B_filterFALSE.csv")
      for (nm in names(arm_files)) {
        f <- file.path(output_dir, arm_files[[nm]])
        if (file.exists(f)) {
          sa <- read.csv(f, row.names = 1)
          sa <- sa[rownames(sa) != "other", , drop = FALSE]
          stash_raw(nm, sa)
          stash(nm, rownames(sa), sa$mean, sa$sig.q)
        } else {
          # sparse_L_filter*.csv is absent when the Laplace branch collapsed
          # (fewer than 2 usable draws). Record it as a FAILURE rather than
          # letting the method silently vanish from that replicate -- a missing
          # method is indistinguishable from one that discovered nothing once
          # the results are aggregated.
          # `<-`, not `<<-`: this runs directly in run_all_methods()'s own frame,
          # so `<<-` would skip the local `failures` and write to the global env.
          # (stash()/stash_raw() use `<<-` correctly because they are nested.)
          failures[[nm]] <- "arm table not written (branch collapsed)"
        }
      }
    }

    ## ------------------------------------------------------------------------
    ## SUPPLEMENTARY: the nine-variant grid, 3 arms x 3 bandwidths.
    ## ------------------------------------------------------------------------
    ## The same grid as "Task 3 sparse SSRV variant tables.R" reports on the real
    ## data, and for the same reason: the Laplace variance is O(h^2) in the
    ## bandwidth, so reading a bandwidth sweep at the hybrid rule alone hides the
    ## rule SWITCHING ARMS as h changes, and reading the arms at one bandwidth
    ## alone hides how much of an arm's behaviour is just the bandwidth.
    ##
    ## Cheap because `resam.array` is REUSED -- the Dirichlet resample is the
    ## expensive half (7.2 s of the 13.3 s pair at the baseline) and it does not
    ## depend on the bandwidth. Same reuse pattern as check_sparse_bw_sens() in
    ## ../../sparse_SSRV.R. Measured marginal cost at the baseline: +12.1 s on a
    ## ~433 s replicate.
    ##
    ## Each call writes its own scale_model / sparse_auto / sparse_L / sparse_B /
    ## theta_perp tables under a "_bw_<tag>" suffix. Nothing is stashed as a
    ## METHOD -- these must not enter the main FDR/TPR curves. aggregate_main.R
    ## reads the tagged files into a separate supplementary frame.
    ##
    ## Naive path only: the oracle arm has one arm by construction.
    if (!informed && length(ssrv_bandwidths)) {
      for (bw in ssrv_bandwidths) {
        sfx <- paste0("_bw_", bw)
        timings[[paste0("sparse", sfx)]] <- system.time(tryCatch({
          invisible(sim_sparse_auto(dat.array = resam.array, Y, x, D, N, S, output_dir,
                                    metadata_df = filtered_metadata_df,
                                    filtered.ind = rep(FALSE, D),
                                    bandwidth = bw, file_suffix = sfx,
                                    scale_model = scale_model))
        }, error = function(e) {
          # A bandwidth selector can fail on its own (SJ needs enough distinct
          # points). That must not take down the replicate's MAIN result, which
          # is already written by the call above.
          failures[[paste0("Sparse_SSRV", sfx)]] <<- conditionMessage(e)
        }))
      }
    }
  }

  ##### ALDEx2 #################################################################
  # Two configurations: the original CLR ("clr") and the default scale model ("def").
  # Both are scale-invariant and have no usable offset API here -- the matrix-gamma
  # branch of aldex.clr.function() has a live name mismatch (scale_samples vs
  # scale.samples) -- so they run naive only.
  if ("aldex" %in% methods && !informed) {
    timings$aldex_clr <- system.time(
      afit_ori_raw <- aldex(Y, as.vector(x), denom = "all", mc.samples = 1000, gamma = NULL)
    )
    timings$aldex_def <- system.time(
      afit_def_raw <- aldex(Y, as.vector(x), mc.samples = 1000, gamma = .5)
    )

    for (nm in c("ALDEx2_clr", "ALDEx2_def")) {
      # NB: named aldex_raw, not raw -- `raw` is the keep_raw accumulator in this scope
      # and a local of that name would clobber it through stash_raw's <<-.
      aldex_raw <- if (nm == "ALDEx2_clr") afit_ori_raw else afit_def_raw
      s <- append_sig(summary_aldex2(aldex_raw), function(z) sig_aldex2(z, pval = pval))
      s <- s[!s$category %in% "other", , drop = FALSE]
      # diff.btw is the log2 CLR fold-change; `effect` is standardized and not an LFC.
      dbtw <- aldex_raw[match(s$category, rownames(aldex_raw)), "diff.btw"]
      out <- cbind(s, effect_sizes, diff.btw = dbtw)
      write.csv(out, file.path(output_dir,
        sprintf("%s.out N%sD%s.csv", if (nm == "ALDEx2_clr") "ALDEx2_ori" else "ALDEx2_def", n, d)))
      stash_raw(nm, out)
      stash(nm, s$category, dbtw, s$sig)
    }
  }

  ##### DESeq2 #################################################################
  if ("deseq2" %in% methods) {
    dat <- cbind(Condition = as.numeric(x), (t(Y) + 1))
    timings$deseq2 <- system.time(
      dfit_raw <- if (informed) run_deseq2_offset(dat, size_factors = offsets$size_factor) else run_deseq2(dat)
    )
    # summary_deseq2() already carries log2FoldChange (and sets mean = log2FoldChange).
    dfit <- append_sig(summary_deseq2(dfit_raw), function(z) sig_deseq2(z, pval = pval))
    dfit <- dfit[!dfit$category %in% "other", , drop = FALSE]
    out <- cbind(dfit, effect_sizes)
    write.csv(out, file.path(output_dir, sprintf("deseq2.out N%sD%s.csv", n, d)))
    stash_raw(if (informed) "DESeq2_oracle" else "DESeq2", as.data.frame(dfit_raw))
    stash(if (informed) "DESeq2_oracle" else "DESeq2",
          dfit$category, dfit$log2FoldChange, dfit$sig, qvalue = dfit$padj)
  }

  ##### phyloseq object shared by ANCOM / ANCOM-BC / ANCOM-BC2 #################
  if (any(c("ancom", "ancombc", "ancombc2") %in% methods) && !informed) {
    samples_df <- as.data.frame(cbind(sample = paste0("sample", 1:length(x)), x))
    samples_df <- samples_df %>% tibble::column_to_rownames("sample")
    otu_mat <- as.matrix(Y)
    colnames(otu_mat) <- rownames(samples_df)
    OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
    ANCOM.dat <- phyloseq(OTU, sam_data = sample_data(samples_df))
  }

  if ("ancombc" %in% methods && !informed) {
    tryCatch({
      timings$ancombc <- system.time(ancbcfit <- ancombc(
        data = ANCOM.dat, assay_name = "counts",
        tax_level = NULL, formula = "disease",
        p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
        group = "disease", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
        max_iter = 100, conserve = TRUE, alpha = 0.05, n_cl = 1, verbose = TRUE))

      # as.vector() is load-bearing. With struc_zero = TRUE, ANCOM-BC merges structural
      # zero information into the results, which leaves res$se, res$p_val and res$q_val
      # holding their `disease1` entry as a 1-COLUMN MATRIX whose dimname is
      # "structural_zero (disease = 1)" -- res$lfc, res$W and res$diff_abn stay plain
      # vectors. cbind() honours a matrix's own dimname over the tag you supply, so
      # `adj.p = <that matrix>` silently produced a column named
      # "structural_zero (disease = 1)". Nothing errored there; ancbc.out$adj.p was just
      # NULL, and the failure surfaced later in stash() as
      #     "arguments imply differing number of rows: <n>, 1, 0"
      # (n taxa, the length-1 pval default, the length-0 NULL qval).
      #
      # This is inherited from sparse_analysis_trycatch.R:232-237, which is why the
      # checked-in src/05052025/run_failed already lists ANCOMBC. Stripping dim/dimnames
      # is applied to all three defensively -- se and p_val have the same shape and
      # would break identically if anyone reached for them.
      ancbc.out <- cbind(
        ancbcfit$res$lfc[, c("taxon", "disease1")],
        W = as.vector(ancbcfit$res$W$disease1),
        adj.p = as.vector(ancbcfit$res$q_val$disease1),
        sig = as.vector(ancbcfit$res$diff_abn[, "disease1"])
      )
      stopifnot(all(c("taxon", "disease1", "W", "adj.p", "sig") %in% names(ancbc.out)))
      ancbc.out <- ancbc.out[!ancbc.out$taxon %in% "other", , drop = FALSE]
      write.csv(cbind(ancbc.out, effect_size = effect_sizes),
                file.path(output_dir, sprintf("ANCOMBC.out N%sD%s.csv", n, d)))
      stash_raw("ANCOM-BC", ancbc.out)
      stash("ANCOM-BC", ancbc.out$taxon, ancbc.out$disease1, ancbc.out$sig,
            qvalue = ancbc.out$adj.p)
    }, error = function(e) failures[["ANCOMBC"]] <<- conditionMessage(e))
  }

  if ("ancom" %in% methods && !informed) {
    tryCatch({
      timings$ancom <- system.time(ancfit <- ancom(
        data = ANCOM.dat, assay_name = "counts",
        tax_level = NULL, p_adj_method = "BH",
        prv_cut = 0, lib_cut = 0,
        main_var = "disease", rand_formula = NULL,
        lme_control = NULL, struc_zero = FALSE,
        neg_lb = TRUE, alpha = 0.05, n_cl = 1))

      ares <- ancfit$res[!ancfit$res$taxon %in% "other", , drop = FALSE]
      write.csv(cbind(ares, effect_size = effect_sizes),
                file.path(output_dir, sprintf("ANCOM.out N%sD%s.csv", n, d)))
      # ANCOM reports only the W statistic -- not an LFC, so est is NA by design.
      stash_raw("ANCOM", ares)
      stash("ANCOM", ares$taxon, rep(NA_real_, nrow(ares)), ares$detected_0.7)
    }, error = function(e) failures[["ANCOM"]] <<- conditionMessage(e))
  }

  ##### LinDA ##################################################################
  if ("linda" %in% methods && !informed) {
    tryCatch({
      meta <- as.matrix(x); colnames(meta) <- "disease"
      otu_mat <- as.matrix(Y)
      timings$linda <- system.time(linda.obj <- linda(
        otu_table(otu_mat, taxa_are_rows = TRUE), meta, formula = '~disease', alpha = 0.05,
        prev.cut = 0, lib.cut = 0, winsor.quan = 0.97))

      lo <- linda.obj$output$disease
      lo <- lo[!rownames(lo) %in% "other", , drop = FALSE]
      write.csv(cbind(lo, effect_size = effect_sizes),
                file.path(output_dir, sprintf("LinDA.out N%sD%s.csv", n, d)))
      stash_raw("LinDA", lo)
      stash("LinDA", rownames(lo), lo$log2FoldChange, lo$reject, qvalue = lo$padj)
    }, error = function(e) failures[["LinDA"]] <<- conditionMessage(e))
  }

  ##### edgeR ##################################################################
  if ("edger" %in% methods) {
    tryCatch({
      timings$edgeR <- system.time({
        y <- DGEList(counts = as.matrix(Y + 1), group = c(x))
        if (informed) {
          # Hand edgeR the true scale instead of letting TMM re-estimate one.
          y$samples$lib.size <- offsets$lib_size
          y$samples$norm.factors <- rep(1, ncol(y))
        } else {
          y <- calcNormFactors(y)
        }
        y <- estimateDisp(y)
        et <- exactTest(y)
      })

      # FDR, not the raw p (2026-08-26). Every real-data script calls
      #     topTags(et, n = nrow(Y), p.value = 1) %>% mutate(sig = FDR < 0.05)
      # and topTags() sets FDR = p.adjust(PValue, "BH") over the whole tested
      # table. This block used to threshold PValue directly, so edgeR ran with no
      # multiplicity control in the simulation and with BH on the real data -- the
      # simulation's edgeR FDR was measuring the missing correction, not the method.
      #
      # BH is applied BEFORE "other" is dropped, so the adjustment family is the
      # full table exactly as topTags() would form it. That also matches DESeq2 a
      # few blocks up, whose padj comes from results() over every row.
      et_tab <- et$table
      et_tab$FDR <- p.adjust(et_tab$PValue, method = "BH")
      tb <- et_tab[!rownames(et_tab) %in% "other", , drop = FALSE]
      write.csv(cbind(tb, sig = tb$FDR < 0.05, effect_size = effect_sizes),
                file.path(output_dir, sprintf("edgeR N%sD%s.csv", n, d)))
      stash_raw(if (informed) "edgeR_oracle" else "edgeR", tb)
      stash(if (informed) "edgeR_oracle" else "edgeR", rownames(tb), tb$logFC, tb$FDR < 0.05,
            pvalue = tb$PValue, qvalue = tb$FDR)
    }, error = function(e) failures[["edgeR"]] <<- conditionMessage(e))
  }

  ##### limma-voom #############################################################
  if ("limma" %in% methods) {
    tryCatch({
      X <- x; X[x == 1] <- "T"; X[x == 0] <- "C"
      timings$limma <- system.time({
        design <- model.matrix(~0 + X)
        y <- DGEList(counts = as.matrix(Y + 1), group = X)
        if (informed) {
          y$samples$lib.size <- offsets$lib_size
          y$samples$norm.factors <- rep(1, ncol(y))
        } else {
          y <- calcNormFactors(y, method = "TMM")
        }
        v <- voom(y, design, plot = FALSE)
        fit <- lmFit(v, design)
        contr <- makeContrasts(XT - XC, levels = colnames(coef(fit)))
        contr.fit <- eBayes(contrasts.fit(fit, contr))
      })

      # adj.P.Val, not the raw p (2026-08-26). Every real-data script calls
      #     topTable(contr.fit, p.value = 1, number = nrow(Y)) %>% mutate(sig = adj.P.Val < 0.05)
      # and topTable()'s adj.P.Val is BH over every tested feature. This block
      # thresholded contr.fit$p.value directly -- no multiplicity control at all.
      # Adjusted before "other" is dropped, for the same reason as edgeR above.
      adj.p <- p.adjust(as.vector(contr.fit$p.value), method = "BH")
      limma.out <- as.data.frame(cbind(mean = contr.fit$coefficients,
                                       p.val = contr.fit$p.value,
                                       adj.P.Val = adj.p,
                                       sig = as.logical(adj.p < 0.05)))
      limma.out <- limma.out[!rownames(limma.out) %in% "other", , drop = FALSE]
      colnames(limma.out) <- c("estimate", "p.val", "adj.P.Val", "sig")
      limma.out$sig <- as.logical(limma.out$sig)
      write.csv(cbind(limma.out, effect_sizes = effect_sizes),
                file.path(output_dir, sprintf("limma N%sD%s.csv", n, d)))
      stash_raw(if (informed) "limma_oracle" else "limma", limma.out)
      stash(if (informed) "limma_oracle" else "limma", rownames(limma.out), limma.out$estimate, limma.out$sig,
            pvalue = limma.out$p.val, qvalue = limma.out$adj.P.Val)
    }, error = function(e) failures[["limma"]] <<- conditionMessage(e))
  }

  ##### baySeq #################################################################
  if ("bayseq" %in% methods) {
    tryCatch({
      coldata <- x
      groups <- list(NDE = rep(1, length(coldata)), DE = as.numeric(coldata))
      timings$baySeq <- system.time({
        dds <- new("countData", data = as.matrix(Y + 1), groups = groups, replicates = coldata)
        # Carry taxon names through as annotation. topCounts() returns rows sorted by
        # posterior likelihood, so without this there is no way to map results back to
        # taxa -- the original code did `rownames(res) <- rownames(Y)` on the sorted
        # table, silently mislabelling every taxon.
        dds@annotation <- data.frame(taxon = rownames(Y), stringsAsFactors = FALSE)
        libsizes(dds) <- if (informed) offsets$lib_size else getLibsizes(dds)
        dds <- getPriors.NB(dds, cl = NULL)
        dds <- getLikelihoods(dds, cl = NULL, bootStraps = 3, verbose = FALSE)
        res <- topCounts(dds, group = "DE", number = nrow(Y))
      })

      # Keep `if (...) a else b` on ONE line. Split across lines it parses inside a
      # brace block but is a syntax error at top level, so the moment this block is
      # lifted out to run standalone -- which is what debugging one method looks like --
      # it fails with "unexpected 'else'".
      ord <- if ("taxon" %in% colnames(res)) match(rownames(Y), as.character(res$taxon)) else match(rownames(Y), rownames(res))
      stopifnot(!anyNA(ord))
      res <- res[ord, , drop = FALSE]
      bayseq.out <- data.frame(p.val = res$FWER.DE, sig = res$FWER.DE < 0.05,
                               row.names = rownames(Y))
      bayseq.out <- bayseq.out[rownames(bayseq.out) != "other", , drop = FALSE]
      write.csv(cbind(bayseq.out, effect_size = effect_sizes),
                file.path(output_dir, sprintf("baySeq N%sD%s.csv", n, d)))
      # baySeq reports a posterior probability, not an LFC.
      stash_raw(if (informed) "baySeq_oracle" else "baySeq", bayseq.out)
      stash(if (informed) "baySeq_oracle" else "baySeq", rownames(bayseq.out), rep(NA_real_, nrow(bayseq.out)),
            bayseq.out$sig, pvalue = bayseq.out$p.val)
    }, error = function(e) failures[["baySeq"]] <<- conditionMessage(e))
  }

  ##### ANCOM-BC2 (retry loop for the zero-variance error) #####################
  if ("ancombc2" %in% methods && !informed) {
    tryCatch({
      ancbc2fit <- NULL
      attempts <- 0
      ancombc2_pruned <- character(0)
      while (is.null(ancbc2fit) && attempts < 10) {
        attempts <- attempts + 1
        timings$ancombc2 <- system.time({
          ancbc2fit <- tryCatch({
            ancombc2(data = ANCOM.dat, assay_name = "counts",
                     fix_formula = "disease", rand_formula = NULL,
                     p_adj_method = "BH", pseudo_sens = TRUE,
                     prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                     group = "disease", struc_zero = FALSE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 1, verbose = TRUE,
                     global = FALSE, pairwise = TRUE,
                     dunnet = FALSE, trend = FALSE,
                     iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     lme_control = NULL,
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                     trend_control = NULL)
          }, error = function(e) {
            message <- conditionMessage(e)
            match <- str_match(message,
              "Zero variances have been detected for the following taxa:\\s*(.*)\\n")
            if (!is.na(match[2])) {
              taxa_to_remove <- unlist(strsplit(match[2], ",\\s*"))
              ## LOG IT. This is the one place a method ends up testing a smaller
              ## feature set than the other twelve, which the count input policy
              ## above otherwise rules out -- and it is silent by default. It is
              ## an error path, not a filter, but a run where it fires has an
              ## ANCOM-BC2 column that is not comparable to the rest of the row,
              ## and `cbind(a2, effect_size = effect_sizes)` below will recycle
              ## against the shortened table rather than error.
              message(sprintf(
                "ANCOMBC2: pruning %d zero-variance taxa and retrying (attempt %d): %s",
                length(taxa_to_remove), attempts, paste(taxa_to_remove, collapse = ", ")))
              ancombc2_pruned <<- c(ancombc2_pruned, taxa_to_remove)
              ANCOM.dat <<- prune_taxa(!(taxa_names(ANCOM.dat) %in% taxa_to_remove), ANCOM.dat)
              return(NULL)
            } else {
              stop("Unexpected error: ", message)
            }
          })
        })
      }

      a2 <- ancbc2fit$res[!ancbc2fit$res$taxon %in% "other", , drop = FALSE]
      ## effect_sizes is keyed to the FULL feature set. If the retry loop pruned
      ## anything, a2 is shorter and cbind() would recycle it silently -- so key
      ## by taxon name instead, which is the same fix the header describes for
      ## the truth vector generally.
      write.csv(cbind(a2, effect_size = unname(truth[as.character(a2$taxon)])),
                file.path(output_dir, sprintf("ANCOMBC2.out N%sD%s.csv", n, d)))
      if (length(ancombc2_pruned))
        failures[["ANCOMBC2_pruned"]] <- paste(unique(ancombc2_pruned), collapse = ", ")
      stash_raw("ANCOM-BC2", a2)
      stash("ANCOM-BC2", a2$taxon, a2$lfc_disease1, a2$diff_disease1,
            qvalue = a2$q_disease1)
    }, error = function(e) failures[["ANCOMBC2"]] <<- conditionMessage(e))
  }

  ##### MRCKF ##################################################################
  # Structurally different from every other method here: taxa are the PREDICTORS and
  # disease is the RESPONSE, so it answers "which taxa jointly predict disease" rather
  # than "how much does taxon j change between groups". Included for reference, on the
  # detection tier and the per-taxon-sensitivity tier only -- never on the bias tier.
  # Naive arms only; there is no offset API.
  if ("mrckf" %in% methods && !informed) {
    tryCatch({
      if (!requireNamespace("knockoff", quietly = TRUE) ||
          !requireNamespace("compositions", quietly = TRUE)) {
        stop("knockoff and compositions are required for MRCKF")
      }
      timings$mrckf <- system.time(mr <- mrckf_fit(Y, x, q = pval))
      sel <- setdiff(mr$selected, "other")
      write.csv(data.frame(taxon = names(mr$W), W = unname(mr$W),
                           selected = names(mr$W) %in% sel,
                           effect_size = unname(truth[names(mr$W)])),
                file.path(output_dir, sprintf("MRCKF N%sD%s.csv", n, d)), row.names = FALSE)
      keepn <- names(mr$W)[names(mr$W) != "other"]
      stash_raw("MRCKF", data.frame(taxon = names(mr$W), W = unname(mr$W), selected = names(mr$W) %in% sel))
      stash("MRCKF", keepn, unname(mr$W[keepn]), keepn %in% sel)

      # Second draw on the SAME data. knockoff.filter defaults to create.second_order,
      # which draws knockoff copies at random, so W differs run-to-run on identical
      # input. This gives the noise floor that MRCKF's between-arm movement must be read
      # against -- no other method needs one.
      if (isTRUE(mrckf_reseed)) {
        mr2 <- mrckf_fit(Y, x, q = pval)
        sel2 <- setdiff(mr2$selected, "other")
        stash("MRCKF_reseed", keepn, unname(mr2$W[keepn]), keepn %in% sel2)
      }
    }, error = function(e) failures[["MRCKF"]] <<- conditionMessage(e))
  }

  ##### bookkeeping ############################################################
  # Failures must be LOUD. Every method sits inside a tryCatch, so without this a
  # method that throws simply vanishes from the results -- indistinguishable from one
  # that was never requested, and the aggregator would silently report a smaller method
  # set rather than an error. (The existing pipeline already hits this: the checked-in
  # src/05052025/run_failed lists ANCOMBC2, ANCOM and ANCOMBC.)
  failed_file_path <- file.path(output_dir, "run_failed")
  if (length(failures) > 0) writeLines(names(failures), failed_file_path) else file.create(failed_file_path)

  # machine-readable, with the message -- this is what tells you WHY on Roar
  write.csv(
    if (length(failures) > 0)
      data.frame(method = names(failures), message = unlist(failures), stringsAsFactors = FALSE)
    else data.frame(method = character(0), message = character(0)),
    file.path(output_dir, "method_failures.csv"), row.names = FALSE)

  # surface it in the SLURM .err log too, so a systematic failure is visible without
  # having to go digging through per-replicate directories
  for (nm in names(failures)) warning("method FAILED: ", nm, " -- ", failures[[nm]], call. = FALSE)

  missing <- setdiff(
    c(if ("sparse" %in% methods) if (informed) "Sparse_SSRV_oracle" else "Sparse_SSRV",
      # the two extra arms are expected on every naive sparse run (edit (c)); if
      # one is absent without an entry in `failures`, something wrote no table
      # and reported nothing, which must not pass silently
      if ("sparse" %in% methods && !informed) c("Sparse_SSRV_L", "Sparse_SSRV_B"),
      if ("deseq2" %in% methods) if (informed) "DESeq2_oracle" else "DESeq2"),
    names(results))
  if (length(missing) > 0 && length(failures) == 0)
    warning("method produced no result but reported no error: ",
            paste(missing, collapse = ", "), call. = FALSE)

  if (length(timings) > 0) {
    exe.time <- do.call(rbind, lapply(timings, function(z) as.numeric(z)[1:3]))
    colnames(exe.time) <- c("user.self", "sys.self", "elapsed")
    write.csv(exe.time, file.path(output_dir, "execution time.csv"))
  }

  list(results = results, truth = truth, failures = failures, timings = timings, raw = raw)
}
