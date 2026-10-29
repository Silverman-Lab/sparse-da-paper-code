# mrckf_fit() -- shared by the main worker (sparse_analysis_trycatch.R) and by
# run_all_methods() in the oracle study, so there is ONE definition of how MRCKF is
# invoked. Requires the MRCKF, knockoff and compositions packages.

# MRCKF single-response path, reproduced so the knockoff W statistic survives.
#
# MRCKF() itself returns only selected taxon NAMES: runSingleKF() calls
# knockoff.filter(...) and keeps $selected, discarding $statistic. But knockoff.filter
# returns structure(list(..., statistic = W, threshold = t, selected = selected)), so
# calling it once here yields both the selection MRCKF would have made -- same knockoff
# draw, no divergence -- and a continuous per-taxon quantity for tier 2.
#
# The preprocessing below is MRCKF()'s, reproduced VERBATIM including a defect: the
# package binds `ref` to the RAW last column while the correctly-closed reference is
# computed into `res` and never used, so
#     X[i,t] = alr_true[i,t] - log(rowSum_i)
# i.e. every design column of row i carries the same -log(n_i) offset, which the later
# column-wise scale() cannot remove. MRCKF is therefore NOT scale-invariant as installed.
# It is reproduced rather than fixed on purpose: the oracle study must use the same
# estimator as sparse_MRCKF.R, and a silent divergence would be worse than the defect.
#
# Note MRCKF drops the last column of W as the log-ratio reference. Our Y has "other"
# last, so "other" is the reference and can never be selected -- harmless here since
# "other" is excluded from every metric anyway, but in ckf_simulation.R the same code
# silently makes a real taxon undetectable.
#
# W_j = |beta_j| - |beta_knockoff_j|: an importance statistic, NOT a coefficient and not
# an effect size. It has no directional sign, so METHOD_SCALE is NA and classify()
# waives the sign check.
mrckf_fit <- function(Y, x, q = 0.05) {
  W <- t(Y + 1)                                  # matches sparse_MRCKF.R:111
  pstar <- ncol(W)
  p <- pstar - 1
  ref <- W[, pstar]                              # raw, per the package
  Z <- as.matrix(compositions::acomp(W))
  X <- matrix(0, nrow = nrow(W), ncol = p)
  for (t in seq_len(p)) X[, t] <- log(Z[, t] / ref)
  predictorNames <- colnames(Z)[-pstar]

  # runSingleKF() calls knockoff.filter TWICE -- offset = 1 first, then offset = 0 --
  # and the single-response branch of MRCKF() returns the offset = 0 result. The first
  # call must be reproduced even though its result is discarded, because knockoff
  # construction is random: dropping it would advance the RNG stream differently and
  # yield a different selection from the same seed. Verified by an equivalence test
  # against MRCKF() itself.
  invisible(knockoff::knockoff.filter(X, as.numeric(x),
                                      statistic = knockoff::stat.glmnet_coefdiff,
                                      fdr = q, offset = 1))
  kf <- knockoff::knockoff.filter(X, as.numeric(x),
                                  statistic = knockoff::stat.glmnet_coefdiff,
                                  fdr = q, offset = 0)
  # back onto Y's taxon order; the reference taxon has no statistic
  Wstat <- setNames(rep(NA_real_, nrow(Y)), rownames(Y))
  Wstat[predictorNames] <- kf$statistic
  list(selected = predictorNames[kf$selected], W = Wstat)
}
