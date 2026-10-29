################################################################################
## Measurement uncertainty on Jin's observed theta^perp.
##
## Sourced by both manuscript_tables.R and manuscript_figures.R so the interval
## in the scale-recovery table and the band in the LFC figure are the same
## numbers, computed once.
##
## Jin's "Absolute total abundance" sheet reports, per cell-sample, the mean
## total abundance in copies/mg, the standard deviation over its N = 4 technical
## replicates, and N. The observed
##
##     theta^perp = mean_{VD} log2(s_i) - mean_{VS} log2(s_i)
##
## is therefore a function of eight quantities each known only up to replicate
## noise, and this propagates that noise to an interval on theta^perp.
##
## WHAT THE INTERVAL COVERS. The standard error of each reported mean,
## sd_i / sqrt(N_i) -- how well the BarBIQ assay pinned down that sample's total
## abundance. It does NOT cover between-mouse variation: theta^perp is an average
## over four mice per group, and how much those four mice vary is a separate
## question this sheet does not answer. So the band says "the measurement is at
## least this uncertain", not "the true load shift lies in here".
##
## The lognormal parameterisation matches the one the gold reference model uses
## in sparse_DA_<datestamp>.R, so the same assay noise is described the same way
## in both places rather than two ways.
################################################################################

## One draw of theta^perp per Monte Carlo replicate. `mu` and `se` are on the
## natural (copies/mg) scale; the returned draws are log2 fold changes.
perp_measurement_draws <- function(mu, se, grp, B = 20000) {
  stopifnot(length(mu) == length(se), length(mu) == length(grp))
  ok <- is.finite(mu) & mu > 0
  if (!all(ok)) stop("non-positive or missing total abundance in ", sum(!ok), " sample(s)")
  se[!is.finite(se)] <- 0
  ## Lognormal with the reported mean and standard error. A sample whose se is 0
  ## contributes its point value on every draw.
  sig_log <- sqrt(log(1 + (se^2 / mu^2)))
  mu_log  <- log(mu) - 0.5 * sig_log^2
  d <- vapply(seq_along(mu), function(i)
    if (sig_log[i] > 0) rnorm(B, mu_log[i], sig_log[i]) else rep(log(mu[i]), B),
    numeric(B))                                    # B x n, natural log
  d <- d / log(2)                                  # to log2
  rowMeans(d[, grp == 1, drop = FALSE]) - rowMeans(d[, grp == 0, drop = FALSE])
}

## The observed theta^perp and its measurement interval, one row per cecal
## location. `dir_jin` is Jin's study directory.
jin_scale_uncertainty <- function(dir_jin, level = 0.95, B = 20000, seed = 20260824) {
  f <- file.path(dir_jin, "Supplementary Data 8.xlsx")
  if (!file.exists(f)) return(NULL)
  s <- as.data.frame(readxl::read_xlsx(f, sheet = "Absolute total abundance", skip = 1))
  names(s)[1:4] <- c("sample", "abund", "sd", "n")

  diet <- substr(s$sample, 1, 2)
  locn <- substr(s$sample, 4, 7)
  ## cell-samples only: the "-ecDNA" rows are a different measurement.
  cell <- !grepl("-", s$sample, fixed = TRUE)

  set.seed(seed)
  a <- (1 - level) / 2
  out <- lapply(c("prox", "dist"), function(loc) {
    k <- cell & locn == loc & diet %in% c("VS", "VD")
    if (!any(k)) return(NULL)
    mu <- as.numeric(s$abund[k]); sd_ <- as.numeric(s$sd[k]); n <- as.numeric(s$n[k])
    grp <- as.numeric(diet[k] == "VD")
    se <- sd_ / sqrt(ifelse(is.finite(n) & n > 0, n, 1))
    obs <- mean(log2(mu[grp == 1])) - mean(log2(mu[grp == 0]))
    dr  <- perp_measurement_draws(mu, se, grp, B = B)
    data.frame(dataset = "Jin2022", pair = loc, obs = obs,
               lower = unname(quantile(dr, a)), upper = unname(quantile(dr, 1 - a)),
               n_samples = sum(k), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  rownames(out) <- NULL
  out
}
