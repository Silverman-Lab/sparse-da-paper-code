# tidy_metrics.R -- turn run_all_methods() output into the two per-replicate tables the
# aggregator consumes, and compute the bias decomposition.
#
# Scale handling. `effect_size` from SparseDOSSA2 is a shift in the NATURAL-LOG mean of
# absolute abundance (spike_oneA_metadata computes mu = param["mu"] + metadata %*% eff),
# so the true log2 fold-change is effect_size / log(2). Methods report on mixed scales;
# METHOD_SCALE below records the multiplier that maps each onto log2. Verify these
# empirically with calibrate_lfc_scales.R before trusting them -- the ANCOM-BC family in
# particular is version-dependent.

# multiplier applied to a method's native estimate to obtain log2 units.
# NA means the method reports no fold-change at all and is excluded from bias metrics.
METHOD_SCALE <- c(
  "Sparse_SSRV"        = 1,           # theta is already log2
  "Sparse_SSRV_oracle" = 1,
  "DESeq2"             = 1,           # log2FoldChange
  "DESeq2_oracle"      = 1,
  "edgeR"              = 1,           # logFC
  "edgeR_oracle"       = 1,
  "limma"              = 1,           # contrast coefficient, log2 via voom
  "limma_oracle"       = 1,
  "LinDA"              = 1,           # log2FoldChange
  "ALDEx2_clr"         = 1,           # diff.btw, log2 CLR difference
  "ALDEx2_def"         = 1,
  "ANCOM-BC"           = 1/log(2),    # lfc reported in natural log
  "ANCOM-BC2"          = 1/log(2),
  "MiMIX"              = 1/log(2),    # posterior mean of beta, natural log CLR
  "ANCOM"              = NA,          # W statistic, not an LFC
  "baySeq"             = NA,          # posterior probability
  "baySeq_oracle"      = NA,
  # MRCKF's knockoff W is |beta_j| - |beta_knockoff_j|: an importance statistic, not a
  # coefficient and not an effect size, with no directional sign. It has no coordinate
  # on the log2 LFC axis, so tier 3 skips it and classify() waives the sign check. It
  # DOES participate in tier 2, where only relative movement in its own units matters.
  "MRCKF"              = NA,
  "MRCKF_reseed"       = NA
)

# Methods excluded from tier 2 (per-taxon sensitivity) despite having a continuous
# output. Theirs are bounded or ordinal -- a rank-based W and a probability -- so a
# change standardized by an across-taxa SD is not comparable to the unbounded
# statistics. They appear in tier 1 only.
TIER2_EXCLUDE <- c("ANCOM", "baySeq", "baySeq_oracle")

# TP/FP/TN/FN under ONE rule for every method. The existing aggregators use two
# different conventions (`effect_size == 0` for ALDEx2/DESeq2 vs `lfc*effect <= 0` for
# ANCOM-BC), which puts type-1 error on different denominators per method.
#   TP = significant, truly non-null, and the sign matches
#   FP = significant but not a TP
#   TN = not significant and truly null
#   FN = not significant and truly non-null
classify <- function(sig, est, true_lfc) {
  sig <- ifelse(is.na(sig), FALSE, sig)
  nonnull <- true_lfc != 0
  sign_ok <- if (all(is.na(est))) rep(TRUE, length(sig)) else (est * true_lfc) > 0
  sign_ok <- ifelse(is.na(sign_ok), FALSE, sign_ok)
  tp <- sig & nonnull & sign_ok
  list(tp = tp, fp = sig & !tp, tn = !sig & !nonnull, fn = !sig & nonnull)
}

f.score <- function(f, precision, recall) {
  (1 + f^2) * (precision * recall) / (f^2 * precision + recall)
}

# One row per (method, taxon) for a single replicate/arm.
tidy_replicate <- function(res, rep_id, arm, params, theta_perp_true) {
  out <- list()
  for (m in names(res$results)) {
    r <- res$results[[m]]
    sc <- if (m %in% names(METHOD_SCALE)) METHOD_SCALE[[m]] else NA_real_
    out[[m]] <- data.frame(
      rep = rep_id, arm = arm,
      N = params$n, D_kept = params$d_kept, lambda = params$lambda, pi = params$pi,
      method = m, taxon = r$taxon,
      est = r$est,
      est_log2_abs = r$est * sc,
      pval = r$pval, qval = r$qval, sig = r$sig,
      true_lfc_loge = r$true_lfc_loge,
      true_lfc_log2 = r$true_lfc_loge / log(2),
      is_spiked = r$true_lfc_loge != 0,
      theta_perp_true = theta_perp_true,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# Collapse a tidy replicate table to one row per (method, arm):
# detection metrics plus the bias decomposition.
#
#   raw_bias = shift_hat + centered_bias
#              ^scale err  ^shape err
#
# shift_hat is the mean estimate over TRUE NULL taxa. A null has absolute LFC exactly 0,
# so any method that recovers the absolute scale has shift_hat = 0, while a purely
# compositional method has shift_hat = -theta_perp_true (every taxon is shifted by the
# same constant, log2 of the between-group total-load ratio). Subtracting it puts all
# methods on a common footing using no oracle information.
summarise_replicate <- function(tidy) {
  do.call(rbind, lapply(split(tidy, tidy$method), function(z) {
    cl <- classify(z$sig, z$est_log2_abs, z$true_lfc_log2)
    ntp <- sum(cl$tp); nfp <- sum(cl$fp); ntn <- sum(cl$tn); nfn <- sum(cl$fn)

    ppv <- ntp / (ntp + nfp)
    tdr <- ntp / (ntp + nfn)

    spiked <- z$is_spiked
    nullish <- !spiked
    has_lfc <- !all(is.na(z$est_log2_abs))

    shift_hat <- if (has_lfc && any(nullish)) mean(z$est_log2_abs[nullish], na.rm = TRUE) else NA_real_
    resid_c <- if (has_lfc) z$est_log2_abs[spiked] - shift_hat - z$true_lfc_log2[spiked] else NA_real_
    resid_r <- if (has_lfc) z$est_log2_abs[spiked] - z$true_lfc_log2[spiked] else NA_real_

    slope <- NA_real_
    if (has_lfc && sum(spiked) > 2) {
      fit <- try(stats::lm(z$est_log2_abs[spiked] ~ z$true_lfc_log2[spiked]), silent = TRUE)
      if (!inherits(fit, "try-error")) slope <- unname(coef(fit)[2])
    }

    tpt <- z$theta_perp_true[1]
    data.frame(
      rep = z$rep[1], arm = z$arm[1], method = z$method[1],
      N = z$N[1], D_kept = z$D_kept[1], lambda = z$lambda[1], pi = z$pi[1],
      # detection
      fdr = nfp / (ntp + nfp),
      tdr = tdr,
      typ1 = nfp / (ntn + nfp),
      ppv = ppv,
      f0.5 = f.score(0.5, ppv, tdr),
      f1 = f.score(1, ppv, tdr),
      # n_sig matters: a method making zero discoveries yields FDR = 0/0 = NaN, and the
      # existing aggregators coerce that to 0, making "discovered nothing" look perfect.
      n_sig = ntp + nfp,
      n_spiked = sum(spiked),
      # bias decomposition, all in log2 units
      shift_hat = shift_hat,
      centered_bias = mean(resid_c, na.rm = TRUE),
      centered_rmse = sqrt(mean(resid_c^2, na.rm = TRUE)),
      raw_bias = mean(resid_r, na.rm = TRUE),
      slope = slope,
      theta_perp_true = tpt,
      # 0 = fully compositional (no scale recovered), 1 = perfect absolute recovery
      scale_recovered = if (is.na(shift_hat) || is.na(tpt) || tpt == 0) NA_real_
                        else 1 - abs(shift_hat) / abs(tpt),
      stringsAsFactors = FALSE
    )
  }))
}


########################################################################################
## TIER 2 -- per-taxon sensitivity
########################################################################################
#
# How much does an INDIVIDUAL TAXON's estimate move between two arms, in the method's
# own units? The primitive is the per-taxon paired difference; everything reported is a
# summary of the distribution of those differences across taxa. This is NOT a comparison
# of two distributions -- that would answer a different (and weaker) question.
#
#   s        = sd(e1) across taxa      the method's own spread on the reference arm
#   delta_j  = (e2_j - e1_j) / s       one value PER TAXON
#
# Standardizing by the across-taxa SD rather than by |e1_j| is deliberate: null taxa sit
# near zero, so a per-taxon relative change would blow up. It also makes the whole tier
# unit-free, which is the point -- no METHOD_SCALE, so a method with an exotic estimand
# (MRCKF's knockoff W) can be compared with DESeq2 without pretending they share a scale.
#
# mean_delta vs sd_delta is the informative split:
#   mean large, sd ~ 0   every taxon moved by the same constant -- a re-normalization
#   mean ~ 0, sd large   individual estimates are unstable -- genuine per-taxon sensitivity
# Reporting only rms_delta would conflate the two.

# Summaries of a delta vector. Separated out so the spiked/null split reuses it exactly.
.delta_summary <- function(delta, prefix = "") {
  delta <- delta[is.finite(delta)]
  nm <- function(x) if (nzchar(prefix)) paste0(prefix, "_", x) else x
  if (length(delta) == 0) {
    out <- as.list(rep(NA_real_, 6)); names(out) <- nm(c("n","mean_delta","sd_delta",
      "rms_delta","q90_abs_delta","max_abs_delta"))
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  out <- list(
    length(delta),
    mean(delta),
    if (length(delta) > 1) sd(delta) else NA_real_,
    sqrt(mean(delta^2)),
    unname(quantile(abs(delta), 0.90)),
    max(abs(delta))
  )
  names(out) <- nm(c("n","mean_delta","sd_delta","rms_delta","q90_abs_delta","max_abs_delta"))
  as.data.frame(out, stringsAsFactors = FALSE)
}

# tidy_a1 / tidy_a2: two tidy_replicate() tables for the SAME replicate, different arms.
# Methods are matched after stripping an "_oracle" suffix, so the informed arm's
# DESeq2_oracle pairs with the naive arm's DESeq2.
change_metrics <- function(tidy_a1, tidy_a2) {
  base_name <- function(m) sub("_oracle$", "", m)
  tidy_a1$.m <- base_name(tidy_a1$method)
  tidy_a2$.m <- base_name(tidy_a2$method)
  common <- intersect(unique(tidy_a1$.m), unique(tidy_a2$.m))
  common <- setdiff(common, base_name(TIER2_EXCLUDE))
  if (length(common) == 0) return(NULL)

  do.call(rbind, lapply(common, function(m) {
    a <- tidy_a1[tidy_a1$.m == m, ]
    b <- tidy_a2[tidy_a2$.m == m, ]

    # match by taxon NAME, never by position -- the arms can order taxa differently and
    # a method may drop taxa (prevalence filters inside ANCOM-BC etc.)
    idx <- match(a$taxon, b$taxon)
    ok <- !is.na(idx)
    e1 <- a$est[ok]
    e2 <- b$est[idx[ok]]
    spiked <- a$is_spiked[ok]

    s <- if (sum(is.finite(e1)) > 1) sd(e1, na.rm = TRUE) else NA_real_
    no_signal <- all(is.na(e1)) || all(is.na(e2)) || is.na(s) || s == 0

    if (no_signal) {
      # no continuous output, or a degenerate spread -- return NA rather than a
      # spurious zero, mirroring how classify() waives the sign check
      base <- .delta_summary(numeric(0))
      sp   <- .delta_summary(numeric(0), "spiked")
      nl   <- .delta_summary(numeric(0), "null")
      extra <- data.frame(sd_e1 = NA_real_, cor_pearson = NA_real_,
                          cor_spearman = NA_real_, mean_raw_delta = NA_real_)
    } else {
      delta <- (e2 - e1) / s
      base <- .delta_summary(delta)
      sp   <- .delta_summary(delta[spiked], "spiked")
      nl   <- .delta_summary(delta[!spiked], "null")
      extra <- data.frame(
        sd_e1 = s,
        # secondary and explicitly DISTRIBUTION-level, not per-taxon: does the ordering
        # of taxa survive the arm change at all
        cor_pearson  = suppressWarnings(cor(e1, e2, use = "complete.obs")),
        cor_spearman = suppressWarnings(cor(e1, e2, use = "complete.obs", method = "spearman")),
        # unstandardized, so raw movement in the method's own units is recoverable
        mean_raw_delta = mean(e2 - e1, na.rm = TRUE)
      )
    }

    cbind(data.frame(rep = a$rep[1], arm_from = a$arm[1], arm_to = b$arm[1],
                     method = m, lambda = a$lambda[1], N = a$N[1], D_kept = a$D_kept[1],
                     pi = a$pi[1], stringsAsFactors = FALSE),
          base, sp, nl, extra)
  }))
}


########################################################################################
## TIER 1 -- selection-set concordance
########################################################################################
#
# The only cross-arm comparison a pure variable-selection method (MRCKF) can join, and
# meaningful for everything else too: how much does the SET of discoveries change?
concordance <- function(tidy_a1, tidy_a2) {
  base_name <- function(m) sub("_oracle$", "", m)
  tidy_a1$.m <- base_name(tidy_a1$method)
  tidy_a2$.m <- base_name(tidy_a2$method)
  common <- intersect(unique(tidy_a1$.m), unique(tidy_a2$.m))
  if (length(common) == 0) return(NULL)

  do.call(rbind, lapply(common, function(m) {
    a <- tidy_a1[tidy_a1$.m == m, ]
    b <- tidy_a2[tidy_a2$.m == m, ]
    s1 <- a$taxon[which(a$sig)]
    s2 <- b$taxon[which(b$sig)]
    uni <- length(union(s1, s2))
    data.frame(rep = a$rep[1], arm_from = a$arm[1], arm_to = b$arm[1],
               method = m, lambda = a$lambda[1],
               n_sig1 = length(s1), n_sig2 = length(s2),
               n_shared = length(intersect(s1, s2)),
               # both arms discovering nothing is perfect agreement, not 0/0
               jaccard = if (uni == 0) 1 else length(intersect(s1, s2)) / uni,
               stringsAsFactors = FALSE)
  }))
}
