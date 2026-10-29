################################################################################
## Closed-form Gaussian-KDE derivatives.
##
## DUPLICATED from sparse_SSRV.R on purpose -- this file is sourced on Roar with no
## source() of its own (see simulation/README.md), so it must stay self-contained.
## Keep the two copies in step; the full rationale lives in sparse_SSRV.R.
##
## Short version: the Laplace scale model reads the flatness of the LFC density at
## its mode as the uncertainty in theta^perp, sigma^2 = -1/(log f)''(m). f is the
## Gaussian KDE, which is C-infinity, so that curvature is exact in closed form.
## The previous code took numDeriv::hessian() of approx() -- a LINEAR interpolant,
## whose second derivative is 0 inside a segment and a kink at each knot, so it had
## no curvature to measure. What came back scaled as 1/step, and since numDeriv's
## step is 0.1*|mode| the answer depended on how far the mode sat from zero:
## measured against the closed form, sigma^2 was understated 4.6x on Stammler2016
## and 28x on Jin2022. If every LFC is identical, this returns sigma^2 = h^2 exactly.
kde_derivs <- function(x, X, h){
  u  <- (x - X) / h
  w  <- dnorm(u)
  n  <- length(X)
  c(f  =  sum(w)            / (n * h),
    f1 = -sum(u * w)        / (n * h^2),
    f2 =  sum((u^2 - 1) * w)/ (n * h^3))
}

kde_mode_var <- function(v, dens = NULL){
  v <- v[!is.na(v)]
  if(is.null(dens)) dens <- density(v, na.rm = TRUE)
  h   <- dens$bw
  gap <- dens$x[2] - dens$x[1]
  m0  <- dens$x[which.max(dens$y)]
  m   <- optimize(function(z) kde_derivs(z, v, h)[["f"]],
                  interval = c(m0 - 2 * gap, m0 + 2 * gap), maximum = TRUE)$maximum
  d     <- kde_derivs(m, v, h)
  logf2 <- d[["f2"]] / d[["f"]] - (d[["f1"]] / d[["f"]])^2
  list(mode = m, variance = -1 / logf2, logf2 = logf2)
}

## This file used to also define its own sparse_auto(), a third variant distinct from
## both sparse_auto() and sparse_auto_old() in sparse_SSRV.R: mean-sd selection, SJ
## bandwidth, KDE-grid mode. It was removed as a stale duplicate that shadowed the
## real-data function of the same name.
##
## The only callers were 2015_prasse.../sparse_DA_02252025.R:486 and
## 2016_stammler.../sparse_DA_02252025.R:499,910, both superseded by their
## sparse_DA_03212026.R versions, which source sparse_SSRV.R instead. Those two
## scripts will now error on sparse_auto(); point them at sparse_SSRV.R if they ever
## need to run again, accepting that they will no longer reproduce their old numbers.
##
## sim_sparse_auto() below is the function the simulation actually uses.


################################################################################
## theta^perp SUMMARY -- COPIED VERBATIM from ../../sparse_SSRV.R (2026-08-29).
##
## perp_mix_quantile(), ssrv_perp_table() and ssrv_perp_summary(), unchanged. The
## simulation has to report theta^perp coverage as the SAME quantity the
## real-data analysis reports (theta_perp_check.R scores `covered` off exactly
## these intervals), and re-deriving the summary here is how the two would come
## to disagree while both looking correct.
##
## COPIED, not sourced: this file is self-contained by contract (see shared.R),
## and sourcing sparse_SSRV.R would redefine sparse_auto(), sparse_auto_sens()
## and the rest on top of the simulation's own versions. IF THE ORIGINALS CHANGE,
## CHANGE THESE.
##
## The laplace arm is summarised from the MIXTURE of per-draw (center, var)
## pairs rather than from the jittered draws, so its interval carries no Monte
## Carlo error beyond the S modes -- which is also why it is unaffected by
## scale_draw, and why every theta^perp number recorded before the rnorm(1) fix
## already described the shared-scalar model.
################################################################################

perp_mix_quantile <- function(p, center, var, draws = NULL){
  emp <- function(pp) if (is.null(draws)) NA_real_ else unname(quantile(draws, pp, na.rm = TRUE))
  ok  <- is.finite(center) & is.finite(var) & var > 0
  if(!any(ok)) return(vapply(p, emp, numeric(1)))
  m <- center[ok]; s <- sqrt(var[ok])
  Fmix <- function(q) mean(pnorm(q, m, s))
  lo <- min(m - 10 * s); hi <- max(m + 10 * s)
  vapply(p, function(pp)
    tryCatch(uniroot(function(q) Fmix(q) - pp, c(lo, hi), tol = 1e-9)$root,
             error = function(e) emp(pp)),
    numeric(1))
}


ssrv_perp_table <- function(draws, arm = NA_character_, center = NULL, var = NULL,
                            level = 0.95){
  a  <- (1 - level) / 2
  pr <- c(a, 0.5, 1 - a)
  d  <- draws[is.finite(draws)]
  na_row <- data.frame(arm = arm, n_draws = length(d), mean = NA_real_, sd = NA_real_,
                       median = NA_real_, lower = NA_real_, upper = NA_real_,
                       width = NA_real_, fold_change = NA_real_, fc_lower = NA_real_,
                       fc_upper = NA_real_, mc_se = NA_real_, stringsAsFactors = FALSE)

  ok  <- if (is.null(center) || is.null(var)) logical(0) else
           is.finite(center) & is.finite(var) & var > 0
  mix <- any(ok)

  if(mix){
    cen <- center[ok]
    m   <- mean(cen)
    ## law of total variance over the mixture: E[Var] + Var[E]
    v   <- mean(var[ok]) + mean((cen - m)^2)
    q   <- perp_mix_quantile(pr, center, var, d)
    n   <- sum(ok)
    mc  <- stats::sd(cen) / sqrt(n)
  } else {
    if(length(d) < 2) return(na_row)
    m  <- mean(d); v <- stats::var(d); q <- unname(quantile(d, pr)); n <- length(d)
    mc <- sqrt(v) / sqrt(n)
  }

  data.frame(arm = arm, n_draws = n, mean = m, sd = sqrt(v), median = q[2],
             lower = q[1], upper = q[3], width = q[3] - q[1],
             ## the same interval read as a fold change on the total load, which
             ## is the unit the external measurements are reported in
             fold_change = 2^m, fc_lower = 2^q[1], fc_upper = 2^q[3],
             mc_se = mc, stringsAsFactors = FALSE)
}


ssrv_perp_summary <- function(perp, level = 0.95){
  one <- function(a) ssrv_perp_table(perp[[a]]$draws, a, perp[[a]]$center,
                                     perp[[a]]$var, level = level)
  L <- one("laplace"); B <- one("bootstrap")
  H <- if (identical(perp$scale_model, "bootstrap")) B else L
  H$arm <- "hybrid"
  out <- rbind(H, L, B)
  out$selected <- perp$scale_model
  rownames(out) <- NULL
  out
}

sparse_auto2 <- function(dat.array, Y, x, D, N, S, output_dir, filtered.ind = NULL){
  
  theta_perpB <- matrix(NA, nrow = D, ncol = S)
  theta_perpL <- matrix(NA, nrow = D, ncol = S)
  theta_pll <- matrix(NA, nrow = D, ncol = S)
  thetaB <- matrix(NA, nrow = D, ncol = S)
  thetaL <- matrix(NA, nrow = D, ncol = S)
  theta <- matrix(NA, nrow = D, ncol = S)
  
  log_wn <- matrix(NA, nrow = D, ncol = N)
  theta_pllB <- matrix(NA, nrow = D, ncol = S)
  
  na.ind <- c()
  
  for(s in 1:S){
    
    # remove samples for filtered taxa
    log_wn <- dat.array[s, !filtered.ind, ]
    
    ##### bootstrap method
    ### bootstrap over N
    n.index <- sample(1:N, replace = T)
    
    # keep sampling until both groups have at least one sample
    while(sum(unique(x) %in% x[n.index]) < length(unique(x))){
      n.index <- sample(1:N, replace = T)
    }
    log_wn_bootstrap <- log_wn[,n.index]
    treat.ind <- x[n.index] == 1
    cont.ind <- x[n.index] == 0 
    
    # rowMeans if a group has more than one sample. No rowMeans needed when there's only one sample
    treat_bootstrap <- if (is.matrix(log_wn_bootstrap[,treat.ind])) rowMeans(log_wn_bootstrap[,treat.ind]) else log_wn_bootstrap[,treat.ind]
    cont_bootstrap <- if (is.matrix(log_wn_bootstrap[,cont.ind])) rowMeans(log_wn_bootstrap[,cont.ind]) else log_wn_bootstrap[,cont.ind]
    theta_pll[!filtered.ind,s] = treat_bootstrap - cont_bootstrap
    theta_pll[filtered.ind,s] = NA
    
    ### bootstrap over D
    d.index <- sample(1:D, replace = T)
    theta_pllB[,s] <- theta_pll[d.index,s]
    
    dens.LFC <- density(theta_pllB[,s], na.rm = T)
    
    # find mode
    theta_perpB[,s] <- -dens.LFC$x[which.max(dens.LFC$y)]
    
    # get theta using bootstrap method
    thetaB[,s] <- theta_pll[,s] + theta_perpB[,s]
    
    
    ##### laplace approximation
    theta_pll[!filtered.ind,s] = rowMeans(log_wn[,x==1]) - rowMeans(log_wn[,x==0])
    theta_pll[filtered.ind,s] = NA
    
    # get density
    dens.LFC <- density(theta_pll[,s], na.rm = T)
    log_density <- function(x) {
      approx(dens.LFC$x, log(dens.LFC$y), xout = x, rule = 2)$y
    }
    
    # ## find mode
    # find mode solely using KDE
    # can oversmooth - dirichlet noise might be erased due to this
    # mode[s] <- optimize(f = log_density, interval = range(theta_pll[,s], na.rm = T), maximum = TRUE)$maximum
    # #mode <- dens.LFC$x[which.max(dens.LFC$y)]
    # theta_perpL[,s] <- -mode[s]
    # # variance at the mode
    # hessian <- hessian(func = log_density, x = mode[s])
    # variance[s] <- -1 / hessian
    
    # find mode 2
    mode <- dens.LFC$x[which.max(dens.LFC$y)]
    # more sensitive to dirichlet noise
    # but can stuck at local mode
    # hence find mode using KDE first, and then find mode more precisely within certain range
    mode <- optimize(f = log_density, interval = range(c(mode-1,mode+1), na.rm = T), maximum = TRUE)$maximum
    theta_perpL[,s] <- -mode
    hessian <- hessian(func = log_density, x = mode)
    variance <- -1 / hessian
    
    
    # in case hessian have numerical issue
    if((!is.na(variance)) & variance > 0){
      theta_perpL[,s] <- theta_perpL[,s] + rnorm(D, 0, sd = sqrt(variance))
      thetaL[,s] <- theta_pll[,s] + theta_perpL[,s]
    }else{
      thetaL[,s] <- NA
      na.ind <- c(na.ind, s)
    }
  }  
  
  if(length(na.ind)>0){
    theta <- theta[,-na.ind]
  }
  
  ### compare uncertainty added from Bootstrap method and Laplace approximation
  # get average standard deviation of each rows (sd of theta for each subject)
  uncertaintyB <- mean(apply(thetaB, MARGIN = 1, function(x) sd(x, na.rm = T)), na.rm = T)
  uncertaintyL <- mean(apply(thetaL, MARGIN = 1, function(x) sd(x, na.rm = T)), na.rm = T)
  
  ##### choose the method with bigger uncertainty as final output
  if(uncertaintyB > uncertaintyL){
    print("boot")
    theta <- thetaB
    theta_perp <- theta_perpB
  }else{
    theta <- thetaL
    theta_perp <- theta_perpL
  }
  
  # calculate posterior predictive p-value
  p.val <- apply(theta, 1, FUN = function(x){min(sum(x < 0,na.rm = T), sum(x > 0,na.rm = T)) * 2 / ncol(theta) })
  q.val <- p.adjust(p.val, method = "BH")
  
  # get 95% credible interval and mean
  out <- apply(theta, 1, FUN = function(x){quantile(x, c(0.025, 0.975), na.rm = T)})
  mean <- apply(theta, 1, FUN = function(x){mean(x)})
  out <- cbind(as.data.frame(t(out)), mean)
  
  # decide significance based on multiple testing adjusted p-value
  sig.q <- (q.val < 0.05)
  out <- cbind(out, p.val, q.val, sig.q)
  
  rownames(out) <- rownames(Y)
  
  # in case a list of filtered taxa is given
  if(sum(filtered.ind) == 0){
    filter = "FALSE"
  }else{
    filter = "TRUE"
  }
  
  return(out)
}

################################################################################
## sim_sparse_auto: simulation counterpart of sparse_auto() in sparse_SSRV.R.
##
## SELECTION RULE UPDATED to match sparse_auto(). It used to keep whichever branch
## had the larger mean per-taxon sd of theta; it now keeps whichever scale model
## CONTRIBUTES more variance, measured per taxon d across the S draws as
##
##     v(d) = Var(theta_d) - Var(theta^||_d)
##          = Var(theta^perp) + 2Cov(theta^||_d, theta^perp)
##
## an exact identity, not an approximation. Subtracting Var(theta^||) is the point:
## the bootstrap branch resamples the N subjects while the laplace branch uses the
## plug-in LFC, so the old comparison of theta alone charged the bootstrap arm for
## resampling variance that has nothing to do with its scale model. scale_rule sets
## how v(d) is summarised across taxa -- see sparse_SSRV.R for the full note.
##
## The KDE BANDWIDTH also changed, from bw = "SJ" to R's default nrd0, to match
## sparse_auto(). That divergence predates the rule change and was never deliberate.
## It is not part of the selection rule, but it does feed it: the Laplace variance is
## -1/hessian of the log-KDE at the mode, so bandwidth sets the curvature that
## variance is read off.
##
## Direction, measured rather than assumed. On these LFC densities -- unimodal, so
## nrd0's 0.9 shrinkage factor bites and SJ has no multimodality to resolve -- nrd0
## is the NARROWER of the two, by roughly 20-25% across D and effect spread. The
## sharper peak shrinks the Laplace variance by 1.2-2.6x, so this change tends to
## favour the BOOTSTRAP arm. (SJ is narrower than nrd0 on multimodal or heavy-tailed
## densities, which is the usual rule of thumb; it does not apply here.)
##
## The KDE-grid mode is left alone -- sparse_auto() refines it with optimize(), this
## does not. That is the one remaining difference between the two functions.
##
## Defaults scale_rule = "majority" and var_max = 100 are the ones the real-data
## analysis and simulation/scale_choice/run_scale_choice.R already use.
##
## Two fixes ride along, both required for the new rule to be computable at all:
##   * failed laplace draws are now dropped from thetaL (the old `theta[,-na.ind]`
##     acted on the still-empty `theta`, so the NA columns survived into thetaL);
##   * the sign-based p-value now passes na.rm, so one failed draw no longer turns
##     every p-value into NA.
################################################################################
## simulation_rerun EDIT (b): scale_model forces the arm instead of selecting it.
##
##   "hybrid"    the selection rule decides -- identical to the original behaviour
##   "laplace"   always use thetaL
##   "bootstrap" always use thetaB
##
## Both thetaB and thetaL are computed unconditionally either way, so forcing costs
## no extra draws and consumes no extra RNG: a forced run is bit-identical to the
## hybrid run of the same seed, differing only in which matrix is carried forward.
## That is what makes "laplace" reproduce the hybrid run's sparse_L_filter*.csv.
##
## scale_stat and the hybrid decision are still computed and recorded under a forced
## value, so scale_model.csv always says what the rule WOULD have chosen.
sim_sparse_auto <- function(dat.array, Y, x, D, N, S, output_dir, metadata_df = NULL, filtered.ind = NULL,
                            scale_rule = c("width", "majority", "median", "q10", "mean"),
                            scale_draw = c("shared", "independent"),
                            var_max = 100,
                            bandwidth = NULL, file_suffix = "",
                            scale_model = c("hybrid", "laplace", "bootstrap")){

  scale_rule  <- match.arg(scale_rule)
  ## "shared" is the model: ONE theta^perp per draw, matching sparse_auto() in
  ## ../../sparse_SSRV.R. "independent" reproduces the pre-2026-08-28 runs, kept
  ## only so those numbers can be regenerated, not because it is a choice.
  scale_draw  <- match.arg(scale_draw)
  scale_model <- match.arg(scale_model)

  ## bandwidth = NULL is R's density() default, nrd0 -- BYTE-IDENTICAL to the
  ## previous behaviour, since density(v) and density(v, bw = NULL) take the same
  ## branch. Pass "nrd", "SJ-ste" etc. for the supplementary variant grid.
  ##
  ## It matters most in the laplace branch: that arm's variance is -1/(log f)''
  ## at the mode, i.e. read straight off the curvature the bandwidth sets, and it
  ## is O(h^2). The bootstrap branch only uses the KDE's argmax, which is far less
  ## sensitive. That asymmetry is why the two knobs are crossed rather than swept
  ## separately -- the rule can switch arms AS h changes.
  ##
  ## file_suffix tags the output CSVs so several bandwidths can share one
  ## output_dir. "" reproduces the existing filenames exactly, so the main path
  ## and everything reading it are untouched.
  dens_bw <- function(v) if (is.null(bandwidth)) density(v, na.rm = TRUE)
                         else                    density(v, bw = bandwidth, na.rm = TRUE)

  theta_pll <- matrix(NA, nrow = D, ncol = S)
  theta_perpB <- matrix(NA, nrow = D, ncol = S)
  theta_perpL <- matrix(NA, nrow = D, ncol = S)
  thetaB <- matrix(NA, nrow = D, ncol = S)
  thetaL <- matrix(NA, nrow = D, ncol = S)
  theta <- matrix(NA, nrow = D, ncol = S)

  # theta_pll is written twice per draw (bootstrapped for thetaB, then plug-in for
  # thetaL). Keep both so Var(theta^||) can be subtracted from the matching branch.
  theta_pllB_store <- matrix(NA, nrow = D, ncol = S)
  theta_pllL_store <- matrix(NA, nrow = D, ncol = S)

  ## Per-draw laplace mode and variance, for ssrv_perp_summary(). Left NA on any
  ## draw the var_max guard rejects, so the summarised mixture is exactly the set
  ## of draws the fit used. n_na_draws / n_var_reject in scale_model.csv say how
  ## many were dropped.
  perp_centerL <- rep(NA_real_, S)
  perp_varL    <- rep(NA_real_, S)

  log_wn <- matrix(NA, nrow = D, ncol = N)
  theta_pllB <- matrix(NA, nrow = D, ncol = S)

  na.ind <- c()
  var.rej <- c()   # subset of na.ind rejected for an implausible laplace variance
  n.resample.rej <- 0   # single-arm subject resamples redrawn; see the guard below

  for(s in 1:S){

    log_wn <- dat.array[s, !filtered.ind, ]

    ## bootstrap over N
    ## simulation_rerun EDIT (d): reject a resample that lost an entire arm.
    ##
    ## sample(1:N, replace = TRUE) can return indices from ONE group only. Then
    ## log_wn_bootstrap[, treat.ind2] is a D x 0 matrix, rowMeans() of it returns
    ## NaN rather than erroring, theta_pll[, s] becomes all-NaN, and the failure
    ## surfaces two lines later as
    ##     density.default: need at least 2 points to select a bandwidth
    ## which takes down the whole replicate -- sim_sparse_auto() is called outside
    ## the per-method tryCatch, so every other method in that replicate dies too.
    ##
    ## This is the guard sparse_oracle_scale.R:92-99 already documents as missing
    ## here; the code below is that same guard, so the two paths now agree.
    ##
    ## Why it matters now and did not before. The chance a replicate hits this is
    ## 1 - (1 - 2*(1/2)^N)^S, which at S = 2000 is
    ##     N = 10  ~98%      N = 16  ~6%      N = 20  ~0.4%      N >= 40  negligible
    ## The HPC sweep ran N >= 40 and never saw it. This re-run deliberately moves
    ## to the real-data regime (N = 10-16 in four of the six studies), where
    ## without this guard the method simply cannot be run.
    ##
    ## RNG note: the redraw consumes extra sample() calls, so the stream shifts
    ## from that point on -- but ONLY in the draws that would otherwise have
    ## crashed. At N >= 40 no redraw ever fires and results are unchanged.
    ## n_resample_reject is recorded per replicate so the frequency is auditable
    ## rather than assumed.
    n.index <- sample(1:N, replace = T)
    while(sum(unique(x) %in% x[n.index]) < length(unique(x))){
      n.resample.rej <- n.resample.rej + 1
      n.index <- sample(1:N, replace = T)
    }
    log_wn_bootstrap <- log_wn[,n.index]
    treat.ind2 <- x[n.index] == 1
    cont.ind2 <- x[n.index] == 0 
    
    treat_bootstrap <- if (is.matrix(log_wn_bootstrap[,treat.ind2])) rowMeans(log_wn_bootstrap[,treat.ind2]) else log_wn_bootstrap[,treat.ind2]
    cont_bootstrap <- if (is.matrix(log_wn_bootstrap[,cont.ind2])) rowMeans(log_wn_bootstrap[,cont.ind2]) else log_wn_bootstrap[,cont.ind2]
    theta_pll[!filtered.ind,s] = treat_bootstrap - cont_bootstrap
    theta_pll[filtered.ind,s] = NA
    theta_pllB_store[,s] <- theta_pll[,s]      # theta^|| going into thetaB

    ## bootstrap over D
    d.index <- sample(1:D, replace = T)
    theta_pllB[,s] <- theta_pll[d.index,s]

    # Bandwidth is R's default nrd0, matching sparse_auto() in sparse_SSRV.R. This
    # used to pass bw = "SJ", which had diverged from the real-data analysis.
    # dens_bw() keeps that default and lets the variant grid override it.
    dens.LFC <- dens_bw(theta_pllB[,s])
    theta_perpB[,s] <- -dens.LFC$x[which.max(dens.LFC$y)]

    thetaB[,s] <- theta_pll[,s] + theta_perpB[,s]

    ## laplace approximation

    theta_pll[!filtered.ind,s] = rowMeans(log_wn[,x==1]) - rowMeans(log_wn[,x==0])
    theta_pll[filtered.ind,s] = NA
    theta_pllL_store[,s] <- theta_pll[,s]      # theta^|| going into thetaL

    # nrd0 again -- see the note in the bootstrap branch above. This one matters more:
    # the Laplace variance is -1/hessian of the log-KDE at the mode, so the bandwidth
    # sets the curvature that variance is read off. nrd0 is the narrower choice here,
    # so it shrinks that variance relative to the old SJ -- see the header note.
    dens.LFC <- dens_bw(theta_pll[,s])

    # mode and Laplace variance from the KDE's exact derivatives; replaces the
    # grid max + numDeriv::hessian(approx(...)) pair -- see kde_mode_var() above
    mv <- kde_mode_var(theta_pll[,s], dens.LFC)
    mode <- mv$mode
    theta_perpL[,s] <- -mode
    variance <- mv$variance

    ## Keep the (center, var) pair. ssrv_perp_table() summarises the laplace arm
    ## from the MIXTURE of these, not from the jittered draws, so the interval it
    ## reports is exact given the S modes. They were computed and thrown away
    ## before, which is the only reason the simulation had no theta^perp estimate
    ## to score coverage against.
    ##
    ## Recorded only for draws the var_max guard ACCEPTS -- see the guard below.
    ## A rejected draw's variance is finite and positive (that is the whole point
    ## of var_max: it catches a locally flat log-density, not a negative one), so
    ## ssrv_perp_table()'s own `is.finite(var) & var > 0` filter would NOT drop it
    ## and one numerically degenerate draw would widen the reported interval by
    ## orders of magnitude. The summary must describe the draws the fit kept.

    # # different way of finding mode
    # # first find mode using KDE
    # mode <- dens.LFC$x[which.max(dens.LFC$y)]
    # # but KDE can oversmooth density - can flush away dirichlet noise
    # # use optimization - which is more sensitive, but can stuck at local mode, not global mode
    # # hence, use optimization but using interval that include the mode found via KDE
    # mode <- optimize(f = log_density, interval = range(c(mode-0.5,mode+0.5), na.rm = T), maximum = TRUE)$maximum
    # theta_perpL[,s] <- -mode
    # hessian <- hessian(func = log_density, x = mode)
    # variance <- -1 / hessian
    
    
    # in case hessian have numerical issue.
    # The var_max arm matters as much as the other two: when the numeric hessian at
    # the mode is ~0 the log-density is locally flat and -1/hessian explodes. That is
    # a numerical failure, not uncertainty -- and unlike the old mean-sd rule, a
    # variance-based comparison is not robust to it. One such draw can multiply
    # Var(theta^perp) by orders of magnitude while shifting the sign-based p-value by
    # at most 1/S, handing the scale model to the laplace branch for no real reason.
    # var_max = 100 is sd = 10 on the log2 scale, a 2^10 fold change. Pass Inf to
    # recover the old accept-any-positive-variance behaviour.
    if((!is.na(variance)) & variance > 0 & variance < var_max){
      ## rnorm(1), NOT rnorm(D). theta^perp is ONE number -- the log2 fold change
      ## of the total load -- shared by every taxon, so the draw is a scalar that
      ## recycles across the D rows. Drawing D independent values leaves each
      ## taxon's MARGINAL variance unchanged, which is why it went unnoticed, but
      ## it zeroes the cross-taxon correlation the scale model exists to induce
      ## (cor 0.996 shared vs 0.001 independent at D = 200) and shrinks the
      ## community-wide shift the posterior represents by sqrt(D). The bootstrap
      ## branch has always used one scalar per draw.
      ##
      ## Ported from sparse_auto() in ../../sparse_SSRV.R (2026-08-28) so the
      ## simulation fits the model the real-data analysis fits. On the real data
      ## this changed nothing on five of six studies; here it is expected to
      ## matter more, because the "width" rule compares the two arms by their
      ## DELIVERED interval width and the laplace arm's width is exactly what a
      ## shared draw changes.
      theta_perpL[,s] <- theta_perpL[,s] +
        rnorm(if (scale_draw == "shared") 1 else D, 0, sd = sqrt(variance))
      thetaL[,s] <- theta_pll[,s] + theta_perpL[,s]
      perp_centerL[s] <- -mode          # accepted draws only; see the note above
      perp_varL[s]    <- variance
    }else{
      thetaL[,s] <- NA
      na.ind <- c(na.ind, s)
      if(!is.na(variance) && variance >= var_max) var.rej <- c(var.rej, s)
    }
  }

  # Old mean-sd statistic, kept only as a diagnostic so a re-run can be compared
  # against the previous results without re-fitting. It selects nothing now.
  uncertaintyB <- mean(apply(thetaB, MARGIN = 1, function(x) sd(x, na.rm = T)), na.rm = T)
  uncertaintyL <- mean(apply(thetaL, MARGIN = 1, function(x) sd(x, na.rm = T)), na.rm = T)
  chosen_old   <- if(!is.na(uncertaintyB) && !is.na(uncertaintyL) && uncertaintyB > uncertaintyL) "Bootstrap" else "Laplace"

  # Drop failed draws from thetaL AND from its theta^|| store, so the two stay
  # column-aligned and the variances below are taken over the same draws.
  if(length(na.ind)>0){
    thetaL <- thetaL[,-na.ind, drop = FALSE]
    theta_pllL_store <- theta_pllL_store[,-na.ind, drop = FALSE]
    theta_perpB_keep <- theta_perpB[,-na.ind, drop = FALSE]
  }else{
    theta_perpB_keep <- theta_perpB
  }
  laplace_usable <- ncol(thetaL) >= 2

  ### variance theta^perp contributes to theta, taxon by taxon
  ##
  ## COMMON theta^||. Both branches are evaluated against
  ## the SAME theta^|| -- the plug-in one -- so the statistic reflects only the
  ## scale model:
  ##
  ##     v_arm(d) = Var_s( theta^||_plug[d,] + theta^perp_arm[d,] )
  ##                - Var_s( theta^||_plug[d,] )
  ##
  ## Why this changed. The branches are built on DIFFERENT theta^||: bootstrap
  ## resamples subjects, laplace uses the plug-in. Subtracting each branch's OWN
  ## Var(theta^||) removed exactly the term in which they differ, so the
  ## comparison was not like for like -- the bootstrap branch's conservativeness
  ## lives in its subject-resampling theta^||, and that was discarded, while the
  ## laplace branch's independent rnorm survived (being uncorrelated with
  ## theta^||). Measured over 1300 replicates under the old form: v_L > v_B in
  ## 100%, including cells where the bootstrap posterior was 20% WIDER and made
  ## fewer discoveries.
  ##
  ## Only the SELECTION changes. thetaB is still fitted on the bootstrapped
  ## theta^|| and thetaL on the plug-in, so a given choice returns the same
  ## posterior as before.
  ##
  ## Note the columns: theta^perp_B is subset by na.ind to match, so both arms
  ## are compared over the SAME draws.
  ##
  ## The fallback when the laplace branch has collapsed is the OLD form. It is
  ## unreachable in practice -- laplace_usable = FALSE forces the bootstrap arm a
  ## few lines below without consulting the statistic -- but vB has to be defined
  ## for the record, and Var(theta^||_plug) is not usable then.
  varPll_plug <- if(laplace_usable) apply(theta_pllL_store, 1, var, na.rm = TRUE) else rep(NA_real_, D)
  vB <- if(laplace_usable){
    apply(theta_pllL_store + theta_perpB_keep, 1, var, na.rm = TRUE) - varPll_plug
  }else{
    apply(thetaB, 1, var, na.rm = TRUE) - apply(theta_pllB_store, 1, var, na.rm = TRUE)
  }
  vL <- if(laplace_usable){
    apply(thetaL, 1, var, na.rm = TRUE) - varPll_plug
  }else{
    rep(NA_real_, D)
  }

  ## "width" (the DEFAULT): compare the arms by the MEDIAN 95% credible-interval
  ## width of the DELIVERED posterior and take the wider -- i.e. the more
  ## conservative -- one.
  ##
  ## Why this replaced the v-based default. The v rules compare theta^perp
  ## variance against a common plug-in theta^||, and there the two arms are not
  ## on the same footing: theta_perpL is an INDEPENDENT draw, so vL is a clean
  ## var(theta_perpL), while theta_perpB is correlated with theta^|| and vB is
  ## deflated by a covariance term that vL has no counterpart for. Measured over
  ## 1150 replicates: v_L / v_B ~ 2.79, while the delivered interval widths
  ## differ by only ~7% and in the OPPOSITE direction (bootstrap wider in 68%).
  ## The v rules therefore chose laplace in 98.8% of replicates and left 11 of 23
  ## settings above 0.05 FDR (worst 0.206 at pi = 0.02); this rule leaves none
  ## (worst 0.021). All four v rules are kept and still selectable -- they agree
  ## with each other on 100% of replicates, so they are one alternative, not four.
  ci_width <- function(th) {
    if (is.null(dim(th)) || ncol(th) < 2) return(rep(NA_real_, D))
    apply(th, 1, function(x) {
      q <- quantile(x, c(0.025, 0.975), na.rm = TRUE); unname(q[2] - q[1])
    })
  }
  wB <- ci_width(thetaB); wL <- ci_width(thetaL)

  scale_stat <- switch(scale_rule,
    # the wider delivered posterior wins; see the note above
    width    = c(bootstrap = median(wB, na.rm = TRUE),
                 laplace   = median(wL, na.rm = TRUE)),
    # paired: for what fraction of taxa does each branch contribute more?
    majority = c(bootstrap = mean(vB >  vL, na.rm = TRUE),
                 laplace   = mean(vL >= vB, na.rm = TRUE)),
    median   = c(bootstrap = median(vB, na.rm = TRUE),
                 laplace   = median(vL, na.rm = TRUE)),
    q10      = c(bootstrap = unname(quantile(vB, 0.1, na.rm = TRUE)),
                 laplace   = unname(quantile(vL, 0.1, na.rm = TRUE))),
    mean     = c(bootstrap = mean(vB, na.rm = TRUE),
                 laplace   = mean(vL, na.rm = TRUE)))

  ## ALL FOUR scale_rule variants, scored from this one fit. The rule only picks
  ## between two arms that are both already computed and both already written to
  ## disk, so an alternative rule needs no refitting -- but vB and vL are
  ## per-taxon VECTORS that nothing else saves, so the summaries have to be taken
  ## here or they are gone. Purely additive: nothing below reads these, and they
  ## consume no RNG, so the selected arm and every fitted value are unchanged.
  ## width90mean: the rule the AOAS manuscript's simulation actually used --
  ## MEAN over taxa of the 5%-95% (90%) interval width, wider arm wins. It differs
  ## from the live `width` rule on two axes at once (90% vs 95% interval, mean vs
  ## median across taxa), so neither `width` nor the v-rules stand in for it.
  ## Recorded only; nothing below reads it and it consumes no RNG.
  w90 <- function(th) {
    if (is.null(dim(th)) || ncol(th) < 2) return(rep(NA_real_, D))
    apply(th, 1, function(x) {
      q <- quantile(x, c(0.05, 0.95), na.rm = TRUE); unname(q[2] - q[1])
    })
  }
  w90B <- w90(thetaB); w90L <- w90(thetaL)

  all_rule_stats <- list(
    width90mean = c(mean(w90B, na.rm = TRUE), mean(w90L, na.rm = TRUE)),
    majority = c(mean(vB >  vL, na.rm = TRUE), mean(vL >= vB, na.rm = TRUE)),
    median   = c(median(vB, na.rm = TRUE), median(vL, na.rm = TRUE)),
    q10      = c(unname(quantile(vB, 0.1, na.rm = TRUE)),
                 unname(quantile(vL, 0.1, na.rm = TRUE))),
    mean     = c(mean(vB, na.rm = TRUE), mean(vL, na.rm = TRUE)))
  ## Same fallback as chosen_hybrid below: an unusable laplace branch is Bootstrap
  ## under every rule, so the recorded alternatives stay comparable to the live one.
  chosen_by_rule <- vapply(all_rule_stats, function(st)
    if(!laplace_usable || is.na(st[2]) || st[1] > st[2]) "Bootstrap" else "Laplace",
    character(1))

  ##### choose the scale model contributing the larger variance
  # Unlike sparse_auto(), a laplace branch with fewer than 2 usable draws falls back
  # to the bootstrap rather than stopping: this runs unattended over hundreds of
  # replicates, and sim_sparse_auto() is called outside the tryCatch blocks, so an
  # error here would take every other method in the replicate down with it.
  chosen_hybrid <- if(!laplace_usable || is.na(scale_stat[2]) || scale_stat[1] > scale_stat[2]) "Bootstrap" else "Laplace"

  # EDIT (b). A forced value overrides the rule, EXCEPT that a collapsed laplace
  # branch cannot be forced -- thetaL has fewer than 2 usable columns, so the
  # quantiles and the p-value denominator below are not defined. Falling back
  # silently would mislabel a bootstrap fit as a laplace one, so the fallback is
  # recorded in scale_model.csv (forced_ok = FALSE) rather than hidden.
  forced_ok <- TRUE
  if(scale_model == "hybrid"){
    chosen <- chosen_hybrid
  }else if(scale_model == "bootstrap"){
    chosen <- "Bootstrap"
  }else{
    if(laplace_usable){
      chosen <- "Laplace"
    }else{
      chosen <- "Bootstrap"
      forced_ok <- FALSE
      warning("scale_model = 'laplace' requested but the laplace branch collapsed; used bootstrap")
    }
  }
  theta <- if(chosen == "Bootstrap") thetaB else thetaL
  print(chosen)
  print(c(rule = scale_rule, round(scale_stat, 4)))
  if(!laplace_usable){
    message(sprintf("note: laplace branch retained %d usable draws of %d; falling back to bootstrap",
                    ncol(thetaL), S))
  }else if(length(var.rej) > 0){
    message(sprintf("note: %d of %d draws rejected for variance >= var_max (%g); %d rejected in total",
                    length(var.rej), S, var_max, length(na.ind)))
  }

  # Record which scale model won. One row per replicate: the comparison is a single
  # summary over taxa, so there is one decision per call. margin says whether the
  # call was decisive or a coin-flip -- under scale_rule = "majority" it is a
  # difference of proportions, so it lives on (-1, 1) and is 0 at a dead tie.
  # chosen_old / uncertainty* carry the superseded mean-sd rule for comparison.
  write.csv(
    data.frame(
      chosen         = chosen,
      # what the RULE would have picked, regardless of scale_model. Under
      # scale_model = "hybrid" these two agree by construction; under a forced
      # value chosen_hybrid is the column the goal-1 proportions are computed from.
      chosen_hybrid  = chosen_hybrid,
      scale_model    = scale_model,
      forced_ok      = forced_ok,
      rule           = scale_rule,
      stat_bootstrap = unname(scale_stat[1]),
      stat_laplace   = unname(scale_stat[2]),
      margin         = unname(scale_stat[1] - scale_stat[2]),
      med_vB         = median(vB, na.rm = TRUE),
      med_vL         = median(vL, na.rm = TRUE),
      # the four scale_rule variants scored on this same fit (see above)
      chosen_width90mean = chosen_by_rule[["width90mean"]],
      mean_w90B          = mean(w90B, na.rm = TRUE),
      mean_w90L          = mean(w90L, na.rm = TRUE),
      chosen_majority = chosen_by_rule[["majority"]],
      chosen_median   = chosen_by_rule[["median"]],
      chosen_q10      = chosen_by_rule[["q10"]],
      chosen_mean     = chosen_by_rule[["mean"]],
      mean_vB        = mean(vB, na.rm = TRUE),
      mean_vL        = mean(vL, na.rm = TRUE),
      q10_vB         = unname(quantile(vB, 0.1, na.rm = TRUE)),
      q10_vL         = unname(quantile(vL, 0.1, na.rm = TRUE)),
      chosen_old     = chosen_old,
      uncertaintyB   = uncertaintyB,
      uncertaintyL   = uncertaintyL,
      n_na_draws     = length(na.ind),
      n_var_reject   = length(var.rej),
      # how often the subject bootstrap lost an arm and had to be redrawn --
      # 0 at N >= 40, rises steeply as N falls. See EDIT (d).
      n_resample_reject = n.resample.rej,
      var_max        = var_max,
      S              = S
    ),
    paste0(output_dir, "/scale_model", file_suffix, ".csv"),
    row.names = FALSE
  )

  ## ---------------------------------------------------------------------------
  ## theta^perp, one row per arm -- the estimate the coverage study scores.
  ## ---------------------------------------------------------------------------
  ## Built with ssrv_perp_summary(), the SAME function the real-data analysis
  ## uses, so simulated coverage and measured coverage are the same quantity.
  ##
  ## theta_perpB[1, ] and not colMeans: the bootstrap arm's theta^perp is ONE
  ## scalar per draw recycled across the D rows, so row 1 IS the draw vector.
  ## Under scale_draw = "shared" the laplace arm is now the same shape, but it is
  ## summarised from (center, var) regardless -- see the header on the ported
  ## block above.
  ##
  ## The hybrid row is a COPY of whichever arm the rule chose, never a
  ## re-summary, so the two can never disagree.
  perp_summary <- tryCatch(
    ssrv_perp_summary(list(
      laplace   = list(draws  = as.vector(theta_perpL[1, ]),
                       center = perp_centerL, var = perp_varL),
      bootstrap = list(draws  = as.vector(theta_perpB[1, ]),
                       center = NULL,        var = NULL),
      scale_model = if (chosen == "Bootstrap") "bootstrap" else "laplace")),
    error = function(e) {
      warning("theta^perp summary failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  if (!is.null(perp_summary)) {
    perp_summary$scale_model_arg <- scale_model
    perp_summary$bandwidth       <- if (is.null(bandwidth)) "nrd0" else bandwidth
    write.csv(perp_summary, paste0(output_dir, "/theta_perp", file_suffix, ".csv"),
              row.names = FALSE)
  }

  # calculate posterior predictive p-value
  p.val <- apply(theta, 1, FUN = function(x){min(sum(x < 0, na.rm = T), sum(x > 0, na.rm = T)) * 2 / ncol(theta) })
  q.val <- p.adjust(p.val, method = "BH")

  out <- apply(theta, 1, FUN = function(x){quantile(x, c(0.025, 0.975), na.rm = T)})
  mean <- apply(theta, 1, FUN = function(x){mean(x)})
  out <- cbind(as.data.frame(t(out)), mean)
  sig <- (out$`2.5%` * out$`97.5%` > 0)
  sig.q <- (q.val < 0.05)
  out <- cbind(out, sig, sig.q)
  rownames(out) <- rownames(Y)
  
  if(!is.null(metadata_df)){
    #spiked.ind <- as.numeric(str_remove_all(metadata_df$feature_spiked,"Feature"))
    matched_indices <- match(as.character(metadata_df$feature_spiked), rownames(Y), nomatch=0)
    #spiked.ind <- rownames(Y) %in% metadata_df$feature_spiked
    effect_sizes <- rep(0, D)
    effect_sizes[matched_indices] <- metadata_df$effect_size
    
    out <- cbind(out, effect_sizes)
  }
  
  if(sum(filtered.ind) == 0){
    filter = "FALSE"
  }else{
    filter = "TRUE"
  }
  
  #return(out)
  write.csv(out, paste0(output_dir, "/sparse_auto_filter", filter, file_suffix, ".csv"))
  
  
  ################
  # Laplace-only table. thetaL has had its failed draws removed by now, so ncol() is
  # the number of USABLE draws and is the right p-value denominator. Skipped entirely
  # when the branch collapsed -- writing a table of NAs would look like a result.
  if(laplace_usable){
    # calculate posterior predictive p-value
    p.val <- apply(thetaL, 1, FUN = function(x){min(sum(x < 0, na.rm = T), sum(x > 0, na.rm = T)) * 2 / ncol(thetaL) })
    q.val <- p.adjust(p.val, method = "BH")

    out <- apply(thetaL, 1, FUN = function(x){quantile(x, c(0.025, 0.975), na.rm = T)})
    mean <- apply(thetaL, 1, FUN = function(x){mean(x)})
    out <- cbind(as.data.frame(t(out)), mean)
    sig <- (out$`2.5%` * out$`97.5%` > 0)
    sig.q <- (q.val < 0.05)
    out <- cbind(out, sig, sig.q)
    rownames(out) <- rownames(Y)

    if(!is.null(metadata_df)){
      #spiked.ind <- as.numeric(str_remove_all(metadata_df$feature_spiked,"Feature"))
      matched_indices <- match(as.character(metadata_df$feature_spiked), rownames(Y), nomatch=0)
      #spiked.ind <- rownames(Y) %in% metadata_df$feature_spiked
      effect_sizes <- rep(0, D)
      effect_sizes[matched_indices] <- metadata_df$effect_size

      out <- cbind(out, effect_sizes)
    }
    write.csv(out, paste0(output_dir, "/sparse_L_filter", filter, file_suffix, ".csv"))
  }

  
  ################
  # calculate posterior predictive p-value
  p.val <- apply(thetaB, 1, FUN = function(x){min(sum(x < 0, na.rm = T), sum(x > 0, na.rm = T)) * 2 / ncol(thetaB) })
  q.val <- p.adjust(p.val, method = "BH")

  out <- apply(thetaB, 1, FUN = function(x){quantile(x, c(0.025, 0.975), na.rm = T)})
  mean <- apply(thetaB, 1, FUN = function(x){mean(x)})
  out <- cbind(as.data.frame(t(out)), mean)
  sig <- (out$`2.5%` * out$`97.5%` > 0)
  sig.q <- (q.val < 0.05)
  out <- cbind(out, sig, sig.q)
  rownames(out) <- rownames(Y)
  
  if(!is.null(metadata_df)){
    #spiked.ind <- as.numeric(str_remove_all(metadata_df$feature_spiked,"Feature"))
    matched_indices <- match(as.character(metadata_df$feature_spiked), rownames(Y), nomatch=0)
    #spiked.ind <- rownames(Y) %in% metadata_df$feature_spiked
    effect_sizes <- rep(0, D)
    effect_sizes[matched_indices] <- metadata_df$effect_size
    
    out <- cbind(out, effect_sizes)
  }

  # simulation_rerun EDIT (a): the bootstrap-only table is now written, not just
  # returned. It was already fully built above; the original's write.csv sat AFTER
  # return(out) and was commented out, so the file never appeared. Writing it makes
  # the three arms symmetric on disk -- sparse_auto_filter (hybrid),
  # sparse_L_filter (laplace-only), sparse_B_filter (bootstrap-only) -- which is
  # what lets goal 2's three curves come from ONE fit.
  #
  # Purely additive: no draw is taken and no RNG is consumed here, so a replicate
  # run with this file is bit-identical to one run with the original.
  write.csv(out, paste0(output_dir, "/sparse_B_filter", filter, file_suffix, ".csv"))

  # NOTE: what is RETURNED is the BOOTSTRAP-arm table, not the selected one. That
  # is pre-existing behaviour (calibrate_lfc_scales.R depends on it); the selected
  # table is the sparse_auto_filter*.csv written above. The decision rides along as
  # attributes, mirroring sparse_auto(), so callers need not re-read scale_model.csv.
  attr(out, "scale_model")    <- chosen
  attr(out, "scale_model_arg")<- scale_model
  attr(out, "chosen_hybrid")  <- chosen_hybrid
  attr(out, "forced_ok")      <- forced_ok
  attr(out, "scale_rule")     <- scale_rule
  attr(out, "scale_stat")     <- scale_stat
  attr(out, "v_contribution") <- list(bootstrap = vB, laplace = vL)
  attr(out, "n_laplace_fail") <- length(na.ind)
  attr(out, "n_var_reject")   <- length(var.rej)
  attr(out, "n_resample_reject") <- n.resample.rej
  attr(out, "var_max")        <- var_max
  return(out)
}


#different way of choosing bootstrap vs laplace
sim_sparse_auto2 <- function(dat.array, Y, x, D, N, S, output_dir, metadata_df = NULL, filtered.ind = NULL){
  
  theta_pll <- matrix(NA, nrow = D, ncol = S)
  theta_perpB <- matrix(NA, nrow = D, ncol = S)
  theta_perpL <- matrix(NA, nrow = D, ncol = S)
  thetaB <- matrix(NA, nrow = D, ncol = S)
  thetaL <- matrix(NA, nrow = D, ncol = S)
  theta <- matrix(NA, nrow = D, ncol = S)
  
  log_wn <- matrix(NA, nrow = D, ncol = N)
  theta_pllB <- matrix(NA, nrow = D, ncol = S)
  
  na.ind <- c()
  
  for(s in 1:S){
    
    log_wn <- dat.array[s, !filtered.ind, ]
    
    ## bootstrap over N
    n.index <- sample(1:N, replace = T)
    log_wn_bootstrap <- log_wn[,n.index]
    treat.ind2 <- x[n.index] == 1
    cont.ind2 <- x[n.index] == 0 
    
    treat_bootstrap <- if (is.matrix(log_wn_bootstrap[,treat.ind2])) rowMeans(log_wn_bootstrap[,treat.ind2]) else log_wn_bootstrap[,treat.ind2]
    cont_bootstrap <- if (is.matrix(log_wn_bootstrap[,cont.ind2])) rowMeans(log_wn_bootstrap[,cont.ind2]) else log_wn_bootstrap[,cont.ind2]
    theta_pll[!filtered.ind,s] = treat_bootstrap - cont_bootstrap
    theta_pll[filtered.ind,s] = NA
    
    ## bootstrap over D
    d.index <- sample(1:D, replace = T)
    theta_pllB[,s] <- theta_pll[d.index,s]
    
    dens.LFC <- density(theta_pllB[,s], na.rm = T)
    theta_perpB[,s] <- -dens.LFC$x[which.max(dens.LFC$y)]
    
    thetaB[,s] <- theta_pll[,s] + theta_perpB[,s]
    
    ## laplace approximation
    
    theta_pll[!filtered.ind,s] = rowMeans(log_wn[,x==1]) - rowMeans(log_wn[,x==0])
    theta_pll[filtered.ind,s] = NA
    
    dens.LFC <- density(theta_pll[,s], na.rm = T)
    
    log_density <- function(x) {
      approx(dens.LFC$x, log(dens.LFC$y), xout = x, rule = 2)$y
    }
    
    mode <- optimize(f = log_density, interval = range(theta_pll[,s], na.rm = T), maximum = TRUE)$maximum
    theta_perpL[,s] <- -mode
    hessian <- hessian(func = log_density, x = mode)
    variance <- -1 / hessian 
    
    
    # in case hessian have numerical issue
    if((!is.na(variance)) & variance > 0){
      theta_perpL[,s] <- theta_perpL[,s] + rnorm(D, 0, sd = sqrt(variance))
      thetaL[,s] <- theta_pll[,s] + theta_perpL[,s]
    }else{
      thetaL[,s] <- NA
      na.ind <- c(na.ind, s)
    }
  }  
  
  if(length(na.ind)>0){
    theta <- theta[,-na.ind]
  }
  
  # Bootstrap
  # calculate posterior predictive p-value
  p.val <- apply(thetaB, 1, FUN = function(x){min(sum(x < 0), sum(x > 0)) * 2 / ncol(thetaB) })
  q.val <- p.adjust(p.val, method = "BH")
  
  out <- apply(thetaB, 1, FUN = function(x){quantile(x, c(0.025, 0.975), na.rm = T)})
  mean <- apply(thetaB, 1, FUN = function(x){mean(x)})
  out <- cbind(as.data.frame(t(out)), mean)
  sig <- (out$`2.5%` * out$`97.5%` > 0)
  sig.q <- (q.val < 0.05)
  outB <- cbind(out, sig, sig.q)
  rownames(outB) <- rownames(Y)
  
  if(!is.null(metadata_df)){
    matched_indices <- match(as.character(metadata_df$feature_spiked), rownames(Y), nomatch=0)
    effect_sizes <- rep(0, D)
    effect_sizes[matched_indices] <- metadata_df$effect_size
    outB <- cbind(outB, effect_sizes)
  }
  
  # Laplace
  # calculate posterior predictive p-value
  p.val <- apply(thetaL, 1, FUN = function(x){min(sum(x < 0), sum(x > 0)) * 2 / ncol(thetaL) })
  q.val <- p.adjust(p.val, method = "BH")
  
  out <- apply(thetaL, 1, FUN = function(x){quantile(x, c(0.025, 0.975), na.rm = T)})
  mean <- apply(thetaL, 1, FUN = function(x){mean(x)})
  out <- cbind(as.data.frame(t(out)), mean)
  sig <- (out$`2.5%` * out$`97.5%` > 0)
  sig.q <- (q.val < 0.05)
  outL <- cbind(out, sig, sig.q)
  rownames(outL) <- rownames(Y)
  
  if(!is.null(metadata_df)){
    matched_indices <- match(as.character(metadata_df$feature_spiked), rownames(Y), nomatch=0)
    effect_sizes <- rep(0, D)
    effect_sizes[matched_indices] <- metadata_df$effect_size
    outL <- cbind(outL, effect_sizes)
  }
  
  if(sum(outL$sig.q) > sum(outB$sig.q)){
    out <- outB
    print("Bootstrap is more conservative")
  }else{
    out <- outL
    print("Laplace is more conservative")
  }
  
  if(sum(filtered.ind) == 0){
    filter = "FALSE"
  }else{
    filter = "TRUE"
  }
  
  return(list(out, outL, outB))
  #write.csv(out, paste0(output_dir, "/sparse_auto2_filter",filter,".csv"))
  #write.csv(outL, paste0(output_dir, "/sparse_L_filter",filter,".csv"))
  #write.csv(outB, paste0(output_dir, "/sparse_B_filter",filter,".csv"))
}

