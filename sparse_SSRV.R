## Previous version, kept verbatim for reference and reproducibility.
## Selection rule: larger mean 90% interval width of theta.
## Superseded by sparse_auto() below -- see the note there for why.
sparse_auto_old <- function(dat.array, Y, x, D, N, S, output_dir, filtered.ind = NULL){

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
    
    # find mode
    mode <- optimize(f = log_density, interval = range(theta_pll[,s], na.rm = T), maximum = TRUE)$maximum
    theta_perpL[,s] <- -mode
    
    # variance at the mode
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
    thetaL <- thetaL[,-na.ind]
  }
  
  ### compare uncertainty added from Bootstrap method and Laplace approximation
  # get average standard deviation of each rows (sd of theta across samples)
  uncertaintyB_interval <- t(apply(thetaB, MARGIN = 1, function(x) quantile(x,c(0.05, 0.95), na.rm = TRUE)))
  uncertaintyB <- mean(apply(uncertaintyB_interval, MARGIN = 1, function(x) {x[2]-x[1]}), na.rm = TRUE)
  uncertaintyL_interval <- t(apply(thetaL, MARGIN = 1, function(x) quantile(x,c(0.05, 0.95), na.rm = TRUE)))
  uncertaintyL <- mean(apply(uncertaintyL_interval, MARGIN = 1, function(x) {x[2]-x[1]}), na.rm = TRUE)
  
  ##### choose the method with bigger uncertainty as final output
  if(uncertaintyB > uncertaintyL){
    theta <- thetaB
    theta_perp <- theta_perpB
    print("bootstrap")
  }else{
    theta <- thetaL
    theta_perp <- theta_perpL
    print("laplace")
    print(c(uncertaintyB, uncertaintyL))
  }
  
  # calculate posterior predictive p-value
  p.val <- apply(theta, 1, FUN = function(x){min(sum(x < 0,na.rm = T), sum(x > 0,na.rm = T)) * 2 / (ncol(theta)) })
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
## Closed-form Gaussian-KDE derivatives.
##
## The Laplace scale model reads the FLATNESS of the LFC density at its mode as the
## uncertainty in theta^perp: log f(x) ~ log f(m) - (x-m)^2/(2 sigma^2), so
## sigma^2 = -1/(log f)''(m). The object whose flatness that refers to is the KDE
## itself, and R's density(kernel = "gaussian", bw = h) is
##
##     f(x) = (1/(n h)) sum phi((x - X_i)/h)
##
## which is C-infinity. So the curvature exists in closed form and needs no numerical
## differentiation at all:
##
##     f (x) =  (1/(n h))   sum phi(u_i)          u_i = (x - X_i)/h
##     f'(x) = -(1/(n h^2)) sum u_i phi(u_i)
##     f"(x) =  (1/(n h^3)) sum (u_i^2 - 1) phi(u_i)
##     (log f)" = f"/f - (f'/f)^2
##
## WHY THIS REPLACED numDeriv::hessian() ON approx().
## The old code built log_density() with approx(), i.e. LINEAR interpolation of the
## 512-point log-KDE grid. A piecewise-linear function has second derivative 0 inside
## each segment and a kink at each knot, so it has no curvature to measure. What
## numDeriv returned was an artefact of its step size meeting a kink: the second
## difference at the peak knot diverges like 1/h_step (measured: -13 at h_step = 2
## grid cells, -1301 at 0.01 cells) and is exactly 0 strictly inside a segment.
## Since numDeriv's step is 0.1*|mode|, the reported variance depended on how far the
## mode happened to sit from zero. Measured against the closed form at the same mode,
## the old path understated sigma^2 by 4.6x on Stammler2016 and 28x on Jin2022, and
## was accurate on McMurrough2014 only because its mode near 8 gave numDeriv a step
## of ~19 grid cells, wide enough to average over the real peak.
##
## Sanity check that pins the scale: if every LFC is identical, f is a single Gaussian
## of width h and sigma^2 = h^2 exactly. The Laplace uncertainty is therefore O(h^2) --
## the bandwidth is not a nuisance detail here, it sets the scale uncertainty.
kde_derivs <- function(x, X, h){
  u  <- (x - X) / h
  w  <- dnorm(u)
  n  <- length(X)
  c(f  =  sum(w)            / (n * h),
    f1 = -sum(u * w)        / (n * h^2),
    f2 =  sum((u^2 - 1) * w)/ (n * h^3))
}

## Mode of the KDE, and the Laplace variance at that mode.
##
## The grid max seeds the search because which.max() over the whole grid cannot be
## fooled by a local peak the way a bracketed optimize() can; optimize() then refines
## it to the true stationary point of the smooth f within +/-2 grid cells. Expanding
## about an actual stationary point is what the Laplace approximation assumes -- on a
## raw grid the "mode" and the point where f' vanishes need not coincide.
##
## (log f)'' uses the full expression rather than f"/f: after numerical refinement f'
## is near zero but not exactly zero, and the full form is correct either way.
##
## TWO DIFFERENT QUANTITIES, and this returns BOTH. They are not interchangeable.
##
##   variance = -1/(log f)''(m). The FLATNESS of the LFC density at its mode. As
##       h -> 0 it converges to -g(m)/g''(m), a property of the POPULATION of
##       LFCs rather than of how well the mode is pinned down, so it does NOT
##       shrink as D grows. On Jin2022 it tracks the bandwidth almost exactly --
##       sigma = 0.485, 0.808, 1.400 at h = 0.5x, 1x, 2x nrd0, i.e. sigma ~ 2h --
##       which is the sanity check above (all LFCs equal => sigma^2 = h^2) showing
##       up on real data. Read as an uncertainty it is a restatement of a
##       smoothing choice.
##
##   var_mode = the sampling variance of m-hat ITSELF: the Parzen (1962) kernel
##       mode result, sqrt(D h^3)(m-hat - m) asymptotically normal with variance
##       g(m) R(K') / g''(m)^2, where R(K') = int K'(u)^2 du. For the Gaussian
##       kernel R(K') = int phi'(u)^2 du = 1/(4 sqrt(pi)). This one shrinks with D
##       and is stable in h (0.122, 0.128, 0.150 over that same 4x range on
##       Jin2022). It is variance only -- the O(h^2) smoothing bias of the mode is
##       not in it -- and it treats the D LFCs as iid draws from a fixed density,
##       so it covers TAXON sampling and not the shared N subjects behind them.
##
## THE FIT USES `variance`. var_mode is returned for inspection only, and this
## is the measured reason why -- tried across all six scored datasets, three
## bandwidths each, 2026-08-29. Substituting it for `variance` narrows the
## laplace arm until that arm floods: on nrd0 it called 24 taxa on Zemb2020 for
## 1 TP / 23 FP (against 5 -> 1/4 under `variance`), 54 on Jin2022 for 1/53
## (against 4 -> 1/3), and 64 on Prasse2015/QSeq for 10/54 (against 0 -> 0/0).
## The width rule caught every one of those by switching to bootstrap, so the
## DELIVERED hybrid held up -- 6 TP / 10 FP against 4 TP / 7 FP at nrd0 -- but
## theta^perp coverage against the external measurements fell from 4/5 to 2/5.
##
## NEITHER quantity covers the term that dominates the error in theta^perp: the
## mode of the LFC density need not sit at the true load shift. That bias is
## 0.47 to 0.90 log2 on Jin2022/Zemb2020 -- larger than var_mode -- and it is
## not estimable from the LFCs, being exactly the direction the data does not
## identify. `variance` is kept because it is wide enough to absorb it.
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
  list(mode = m, variance = -1 / logf2, logf2 = logf2, h = h,
       var_mode = d[["f"]] / (4 * sqrt(pi) * length(v) * h^3 * d[["f2"]]^2))
}


################################################################################
## ssrv_arm_table: theta draws -> the posterior summary table an arm reports.
##
## Pulled out of sparse_auto() / sparse_auto_sens() so that the hybrid, the
## laplace-only and the bootstrap-only scale models are all summarised by ONE
## piece of code. Before this, only the SELECTED arm got a full table (interval,
## posterior mean, p, q, sig.q) while the other two got a bare logical vector,
## so the unselected arm could be counted but not inspected: you could see that
## laplace called 72 taxa and bootstrap 16, and not which 72, at what effect
## size, at what q.
##
## The arithmetic is unchanged from the block it replaces:
##   - two-sided posterior predictive p-value, divided by the arm's OWN number of
##     retained draws (the laplace arm can be narrower than the bootstrap arm
##     once failed draws are dropped),
##   - BH across taxa, significance at q < 0.05,
##   - 95% equal-tailed credible interval,
##   - posterior mean taken WITHOUT na.rm, so a taxon with no usable draw reports
##     NA instead of a mean over nothing.
##
## An arm with fewer than 2 usable draws returns an all-NA table of D rows rather
## than an empty one, so it reads as "withheld" rather than "found nothing".
################################################################################
ssrv_arm_table <- function(theta, ids, arm = NA_character_){

  if(is.null(dim(theta)) || ncol(theta) < 2){
    ## Built through the same cbind() path as the live branch so the column
    ## names match exactly ("2.5%", "97.5%", mean, p.val, q.val, sig.q) -- code
    ## downstream indexes this table by position as well as by name.
    ci   <- matrix(NA_real_, nrow = 2, ncol = length(ids),
                   dimnames = list(c("2.5%", "97.5%"), ids))
    out  <- as.data.frame(t(ci))
    mean <- rep(NA_real_, length(ids))
    out  <- cbind(out, mean)
    p.val <- rep(NA_real_, length(ids))
    q.val <- rep(NA_real_, length(ids))
    sig.q <- rep(NA, length(ids))
    out   <- cbind(out, p.val, q.val, sig.q)
    rownames(out) <- ids
    attr(out, "scale_model") <- arm
    attr(out, "n_draws")     <- if(is.null(dim(theta))) 0L else ncol(theta)
    return(out)
  }

  # calculate posterior predictive p-value
  p.val <- apply(theta, 1, FUN = function(x){min(sum(x < 0,na.rm = T), sum(x > 0,na.rm = T)) * 2 / (ncol(theta)) })
  q.val <- p.adjust(p.val, method = "BH")

  # get 95% credible interval and mean
  out <- apply(theta, 1, FUN = function(x){quantile(x, c(0.025, 0.975), na.rm = T)})
  mean <- apply(theta, 1, FUN = function(x){mean(x)})
  out <- cbind(as.data.frame(t(out)), mean)

  # decide significance based on multiple testing adjusted p-value
  sig.q <- (q.val < 0.05)
  out <- cbind(out, p.val, q.val, sig.q)

  rownames(out) <- ids
  attr(out, "scale_model") <- arm
  attr(out, "n_draws")     <- ncol(theta)
  out
}

################################################################################
## sparse_arms: the THREE results of one sparse_auto() / sparse_auto_sens() fit.
##
##   sparse_arms(out)$hybrid      the arm the selection rule picked
##   sparse_arms(out)$laplace     the laplace scale model on its own
##   sparse_arms(out)$bootstrap   the bootstrap scale model on its own
##
## All three are the same shape -- the data.frame sparse_auto() has always
## returned -- so anything that worked on the fit works on any arm:
##
##   arms <- sparse_arms(sparse_out)
##   rownames(arms$laplace)[which(arms$laplace$sig.q)]     # WHICH taxa, not how many
##   arms$bootstrap[arms$bootstrap$sig.q %in% TRUE, ]
##
## `hybrid` is the returned object itself with the arm bookkeeping stripped, so
## it always equals whichever of laplace / bootstrap was selected; a run where it
## matches neither means the tables came from different fits.
##
## A result produced before per-arm tables were recorded has no `arm_out`
## attribute. That comes back as hybrid + two NULLs with a warning, rather than
## an error, so an old saved fit can still be loaded and read.
################################################################################
sparse_arms <- function(out){
  arms <- attr(out, "arm_out")
  hyb  <- out
  for(a in c("arm_out", "sig_by_arm", "n_sig")) attr(hyb, a) <- NULL
  if(is.null(arms)){
    warning("this result carries no per-arm tables -- it predates them; re-run ",
            "the fit to get the laplace-only and bootstrap-only results",
            call. = FALSE)
    return(list(hybrid = hyb, laplace = NULL, bootstrap = NULL))
  }
  list(hybrid = hyb, laplace = arms$laplace, bootstrap = arms$bootstrap)
}

################################################################################
## theta^perp -- the SCALE the method estimates, reported rather than left
## inside the fit.
##
## theta_d = theta^||_d + theta^perp. theta^||_d is the log2 fold change of taxon
## d's RELATIVE abundance, which the sequencing data identifies; theta^perp is
## the common shift that turns it into an ABSOLUTE log2 fold change, which the
## data alone does not. It is a scalar -- one shift for every taxon -- and under
## the model it IS the log2 fold change of the total load between the two groups,
##
##     theta^perp  =  mean_{treat} log2(total)  -  mean_{cont} log2(total),
##
## which is exactly the quantity qPCR (Zemb), flow cytometry (Vandeputte) and the
## absolute-total sheet (Jin) measure. Until now it was computed on every draw
## and then thrown away, so the fit could only be judged by its CONSEQUENCES --
## which features came out significant. Keeping it makes the estimate itself
## checkable against those external measurements: bias, RMSE, interval coverage.
##
## The two scale models produce it differently, and the difference is the whole
## point of the selection rule:
##
##   bootstrap  theta^perp[s] = -mode of the KDE of the LFCs, computed on a
##              resample of the D taxa. Its posterior IS the spread of those
##              modes over the S draws; there is no closed form and none is used.
##
##   laplace    theta^perp[s] ~ N(-mode[s], sigma^2[s]), sigma^2 the Laplace
##              variance at the mode. Marginally over draws theta^perp follows the
##              MIXTURE (1/S) sum_s N(-mode[s], sigma^2[s]), so its mean, sd and
##              quantiles are available EXACTLY -- no Monte Carlo error beyond the
##              S modes. sparse_auto() draws ONE N(0, sigma^2) value per draw and
##              adds it to every taxon, so the draws behind the DA table and the
##              mixture summarised here are the same scalar theta^perp.
##
##   hybrid     whichever arm the selection rule picked -- i.e. the arm the
##              returned table came from, so the reported scale always belongs to
##              the reported effect sizes.
##
## Read off a finished fit; nothing here refits anything or consumes a random
## number:
##
##   sparse_perp(out)                 3-row table: hybrid / laplace / bootstrap
##   theta_perp_draws(out, "laplace") that arm's posterior draws
##
## Sign convention, worth stating because it is easy to get backwards: `x` is
## coded 1 = treatment, 0 = control, theta^|| is treatment minus control, and so
## theta^perp is too. A positive theta^perp means the total load is HIGHER in the
## treatment group.
################################################################################

## Quantiles of the Gaussian mixture (1/S) sum_s N(center[s], var[s]) -- the
## laplace arm's exact marginal for theta^perp.
##
## Worth the uniroot: the alternative is the empirical quantile of one draw per
## s, which throws away the closed form the arm already has. Falls back to that
## empirical quantile if the root cannot be bracketed (a degenerate variance,
## say) rather than returning NA, so a summary always comes back.
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

## One arm's theta^perp posterior -> one row of summary.
##
## `center`/`var` present (the laplace arm) switches on the exact mixture above;
## absent (the bootstrap arm) the draws themselves are the posterior and are
## summarised directly. Either way the columns are the same, so the two arms
## stack.
##
## mc_se is the Monte Carlo error of `mean` -- how much of the reported number is
## the sampler rather than the data. The Dirichlet draws are iid, so it is just
## sd/sqrt(n); for the mixture the mean is driven by the S modes, so sd(center)
## is the right numerator there and not the total sd.
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

## Build the hybrid / laplace / bootstrap summary from the pieces sparse_auto()
## keeps. Split out so sparse_auto(), sparse_auto_sens() and sparse_perp() all
## produce the table the same way -- the hybrid row is a COPY of the arm it
## names, never a re-summary, so the two can never disagree.
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

################################################################################
## sparse_perp: the estimated theta^perp of a finished fit, one row per arm.
##
##   sparse_perp(sparse_out)
##   #>        arm n_draws    mean     sd  median   lower   upper  width ...
##   #> 1   hybrid    5000  1.9421 0.4137  1.9388  1.1401  2.7622 1.6221
##   #> 2  laplace    5000  1.9421 0.4137  1.9388  1.1401  2.7622 1.6221
##   #> 3 bootstrap    5000  1.8887 0.2130  1.8901  1.4700  2.3050 0.8350
##
## `hybrid` is whichever arm the selection rule chose, so it always duplicates one
## of the other two rows -- that is the point: it says which scale the reported
## effect sizes are on.
##
## `level` recomputes the interval at a different credible level from the stored
## draws (and, for the laplace arm, from the stored mixture), so a coverage check
## does not need a refit to move off 95%.
##
## A fit produced before theta^perp was recorded returns NULL with a warning
## rather than erroring, so an old saved result still loads.
################################################################################
sparse_perp <- function(out, level = 0.95){
  perp <- attr(out, "theta_perp")
  if(is.null(perp)){
    warning("this result carries no theta^perp -- it predates the recording; ",
            "re-run the fit to get the estimated scale", call. = FALSE)
    return(NULL)
  }
  ## The stored summary is at 0.95. Anything else is rebuilt from the draws.
  if(isTRUE(all.equal(level, 0.95)) && !is.null(attr(out, "perp_summary")))
    return(attr(out, "perp_summary"))
  ssrv_perp_summary(perp, level = level)
}

## One arm's raw theta^perp draws -- what you need to score coverage yourself, or
## to plot the posterior against a measured value. "hybrid" resolves to the arm
## the selection rule picked.
theta_perp_draws <- function(out, arm = c("hybrid", "laplace", "bootstrap")){
  arm  <- match.arg(arm)
  perp <- attr(out, "theta_perp")
  if(is.null(perp)){
    warning("this result carries no theta^perp -- re-run the fit", call. = FALSE)
    return(NULL)
  }
  if(identical(arm, "hybrid")) arm <- perp$scale_model
  perp[[arm]]$draws
}

################################################################################
## sparse_auto: selects the scale model on the variance theta^perp CONTRIBUTES.
##
## The previous rule (sparse_auto_old) compared the mean 90% interval width of
## theta = theta^|| + theta^perp. Two problems:
##
##  1. It compared theta, not theta^perp. theta^|| is not the same in the two
##     branches -- the bootstrap branch resamples the N samples while the
##     laplace branch uses the plug-in LFC -- so the extra resampling variance
##     sat in the bootstrap total regardless of the scale model.
##  2. It averaged over taxa, so a heavy right tail (a few high-variance rare
##     taxa) could decide the comparison while the taxa near significance went
##     the other way.
##
## What is compared here instead, per taxon d, across the S draws:
##
##     v(d) = Var(theta_d) - Var(theta^||_d)
##          = Var(theta^perp) + 2Cov(theta^||_d, theta^perp)
##
## i.e. the variance theta^perp actually adds once its covariance with theta^||
## is accounted for. This is an exact identity, not an approximation, and it
## reduces to Var(theta^perp) whenever the two are independent (the laplace
## branch). They differ in the bootstrap branch because theta^perp is the mode
## of the same resampled LFCs that form theta^||, so the two partly cancel.
##
## scale_rule sets how v(d) is summarised across taxa. The default "majority"
## is a paired per-taxon comparison, so a minority of extreme taxa cannot decide
## the choice; "q10"/"median" are robust unpaired summaries, and "mean"
## reproduces the averaging behaviour of the old rule.
##
## Returns the same data.frame as before -- the SELECTED arm, i.e. the hybrid
## rule's answer -- so every existing caller is unaffected. All THREE results are
## now produced, though, and the two unselected ones ride along:
##
##   arms <- sparse_arms(out)
##   arms$hybrid       what the selection rule chose (== the returned object)
##   arms$laplace      the laplace scale model on its own
##   arms$bootstrap    the bootstrap scale model on its own
##
## each a full table (interval, posterior mean, p.val, q.val, sig.q), so the
## question "WHICH features does each scale model call?" can be answered from a
## finished fit without refitting. The counts alone were never enough for that.
##
## The decision and its inputs ride along as attributes: scale_model, scale_rule,
## scale_stat, v_contribution, n_laplace_fail, plus arm_out / sig_by_arm / n_sig
## for the per-arm results.
################################################################################
sparse_auto <- function(dat.array, Y, x, D, N, S, output_dir, filtered.ind = NULL,
                        scale_rule = c("width", "majority", "median", "q10", "mean"),
                        scale_draw = c("shared", "independent"),
                        var_max = 100){

  scale_rule <- match.arg(scale_rule)
  ## "shared" is the model: ONE theta^perp per draw. "independent" reproduces the
  ## pre-2026-08-28 runs, which gave each taxon its own jitter -- kept only so
  ## those numbers can be regenerated, not because it is a modelling choice.
  scale_draw <- match.arg(scale_draw)

  theta_perpB <- matrix(NA, nrow = D, ncol = S)
  theta_perpL <- matrix(NA, nrow = D, ncol = S)
  theta_pll <- matrix(NA, nrow = D, ncol = S)
  thetaB <- matrix(NA, nrow = D, ncol = S)
  thetaL <- matrix(NA, nrow = D, ncol = S)
  theta <- matrix(NA, nrow = D, ncol = S)

  # theta_pll is written twice per draw (bootstrapped for thetaB, then plug-in
  # for thetaL). Keep both so Var(theta^||) can be subtracted from each branch.
  theta_pllB_store <- matrix(NA, nrow = D, ncol = S)
  theta_pllL_store <- matrix(NA, nrow = D, ncol = S)

  log_wn <- matrix(NA, nrow = D, ncol = N)
  theta_pllB <- matrix(NA, nrow = D, ncol = S)

  # the laplace arm's theta^perp, draw by draw: the centre and the Laplace
  # variance around it. Recorded -- as sparse_auto_sens() has always done --
  # because together they ARE that arm's posterior for the scale, which the fit
  # otherwise discards. Two numeric vectors of length S, no extra computation.
  #
  # READ THE NAME CAREFULLY. mode_estL stores -mode, i.e. theta^perp ITSELF, not
  # the mode of the LFC density. theta^perp = -mode throughout: a community whose
  # relative log fold changes pile up at -1 has a total load 2x HIGHER in the
  # treatment arm, not lower. The name is kept as-is only because
  # sparse_auto_sens() has always spelled it this way and writes it out under that
  # name (as does `mode_est` for the bootstrap branch), so renaming it here would
  # split one quantity across two spellings.
  mode_estL <- rep(NA_real_, S)
  var_estL  <- rep(NA_real_, S)

  na.ind <- c()
  var.rej <- c()   # subset of na.ind rejected for an implausible laplace variance

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
    theta_pllB_store[,s] <- theta_pll[,s]      # theta^|| going into thetaB

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
    theta_pllL_store[,s] <- theta_pll[,s]      # theta^|| going into thetaL

    # get density
    dens.LFC <- density(theta_pll[,s], na.rm = T)

    # mode and Laplace variance from the KDE's exact derivatives -- see the note on
    # kde_mode_var() above for why this replaced numDeriv::hessian() on approx()
    mv <- kde_mode_var(theta_pll[,s], dens.LFC)
    mode <- mv$mode
    theta_perpL[,s] <- -mode
    variance <- mv$variance
    mode_estL[s] <- -mode
    var_estL[s]  <- variance

    # Retained as a safety net rather than a working guard. At a genuine local
    # maximum of the KDE, f" < 0 and (log f)" < 0, so variance is positive by
    # construction -- unlike the numDeriv path this arm should essentially never
    # fire, and n_var_reject > 0 now signals something odd about the density (e.g.
    # optimize() bracketing a saddle), not routine numerical noise.
    # var_max = 100 is sd = 10 on the log2 scale, a 2^10 fold change.

    if((!is.na(variance)) & variance > 0 & variance < var_max){
      ## rnorm(1), NOT rnorm(D). theta^perp is ONE number -- the log2 fold change
      ## of the total load -- shared by every taxon, so the draw is a scalar that
      ## recycles across the D rows. Drawing D independent values leaves each
      ## taxon's MARGINAL variance unchanged, which is why it went unnoticed, but
      ## it zeroes the cross-taxon correlation the scale model exists to induce
      ## (measured at D = 200: cor 0.996 shared vs 0.001 independent) and shrinks
      ## the community-wide shift the posterior represents by sqrt(D) (sd 0.801
      ## vs 0.075 at sigma = 0.8). The bootstrap branch above has always used one
      ## scalar per draw, and ssrv_perp_table() summarises this arm as a scalar
      ## mixture, so both of those now describe the model that is actually fitted.
      theta_perpL[,s] <- theta_perpL[,s] +
        rnorm(if (scale_draw == "shared") 1 else D, 0, sd = sqrt(variance))
      thetaL[,s] <- theta_pll[,s] + theta_perpL[,s]
    }else{
      thetaL[,s] <- NA
      na.ind <- c(na.ind, s)
      if(!is.na(variance) && variance >= var_max) var.rej <- c(var.rej, s)
    }
  }

  # drop failed draws from thetaL AND from its theta^|| store, so the two stay
  # column-aligned and the variances below are taken over the same draws.
  # drop = FALSE: without it a single surviving draw collapses the matrix to a
  # vector and the ncol() check below reads NULL instead of 1.
  if(length(na.ind)>0){
    thetaL <- thetaL[,-na.ind, drop = FALSE]
    theta_pllL_store <- theta_pllL_store[,-na.ind, drop = FALSE]
    # separate variable, NOT a reassignment of theta_perpB: the bootstrap branch
    # keeps all S draws, and only the SELECTION statistic below is restricted to
    # the draws the laplace branch survived. Overwriting theta_perpB would leave
    # a silently truncated matrix behind for anything downstream.
    theta_perpB_keep <- theta_perpB[,-na.ind, drop = FALSE]
  }else{
    theta_perpB_keep <- theta_perpB
  }
  if(is.null(dim(thetaL)) || ncol(thetaL) < 2){
    stop("laplace branch retained fewer than 2 draws; cannot compare scale models")
  }

  ### variance theta^perp contributes to theta, taxon by taxon
  ##
  ## COMMON theta^||. Both branches are evaluated against
  ## the SAME theta^|| -- the plug-in one -- so the statistic reflects only the
  ## scale model:
  ##
  ##     v_arm(d) = Var_s( theta^||_plug[d,] + theta^perp_arm[d,] )
  ##                - Var_s( theta^||_plug[d,] )
  ##              = Var(theta^perp_arm) + 2Cov(theta^||_plug, theta^perp_arm)
  ##
  ## Why this changed. The two branches are built on DIFFERENT theta^||: the
  ## bootstrap branch resamples subjects, the laplace branch uses the plug-in.
  ## Subtracting each branch's OWN Var(theta^||) therefore removed exactly the
  ## term in which they differ, so the comparison was not like for like. The
  ## bootstrap branch's extra conservativeness lives in its subject-resampling
  ## theta^||, and the old subtraction discarded it, while the laplace branch's
  ## independent rnorm survived intact (it is uncorrelated with theta^||). The
  ## result was a statistic that favoured laplace essentially always: measured
  ## over 1300 simulated replicates, v_L > v_B in 100% of them, even at cells
  ## where the bootstrap posterior was 20% WIDER and made fewer discoveries.
  ##
  ## Note this changes only the SELECTION. Each branch is still FITTED on its own
  ## theta^|| -- thetaB keeps the subject bootstrap, thetaL the plug-in -- so the
  ## posterior that is returned is unchanged for a given choice.
  varPll_plug <- apply(theta_pllL_store, 1, var, na.rm = TRUE)
  vB <- apply(theta_pllL_store + theta_perpB_keep, 1, var, na.rm = TRUE) - varPll_plug
  vL <- apply(thetaL, 1, var, na.rm = TRUE) - varPll_plug

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

  ##### choose the scale model contributing the larger variance
  ## `chosen` is what the summary below acts on -- theta / theta_perp are kept
  ## because they name the selected draws for anything debugging inside the
  ## function, not because the returned table is built from them.
  if(scale_stat[1] > scale_stat[2]){
    theta <- thetaB
    theta_perp <- theta_perpB
    chosen <- "bootstrap"
  }else{
    theta <- thetaL
    theta_perp <- theta_perpL
    chosen <- "laplace"
  }
  print(chosen)
  print(c(rule = scale_rule, round(scale_stat, 4)))
  # a non-trivial var.rej count means the laplace approximation is straining on
  # this dataset -- worth looking at before trusting the comparison
  if(length(var.rej) > 0){
    message(sprintf("note: %d of %d draws rejected for variance >= var_max (%g); %d rejected in total",
                    length(var.rej), S, var_max, length(na.ind)))
  }

  ## ALL THREE results, from one code path (see ssrv_arm_table() above).
  ##
  ## outB and outL are the bootstrap-only and laplace-only scale models, each
  ## summarised over its own retained draws. `out` -- the hybrid, i.e. what the
  ## selection rule picked -- is the SAME object as the arm it names rather than
  ## a re-computation, so the selected result can never disagree with that arm's
  ## own table. Both arms were already fitted above; this costs two apply()
  ## passes and no extra sampling.
  ids  <- rownames(Y)
  outB <- ssrv_arm_table(thetaB, ids, "bootstrap")
  outL <- ssrv_arm_table(thetaL, ids, "laplace")
  out  <- if(chosen == "bootstrap") outB else outL

  # in case a list of filtered taxa is given
  if(sum(filtered.ind) == 0){
    filter = "FALSE"
  }else{
    filter = "TRUE"
  }

  ## Every arm's FULL table travels with the result, so a finished run can be
  ## asked which features each scale model called -- with the effect size and q
  ## that go with them -- and not just how many. sparse_arms(out) returns all
  ## three; sig_by_arm / n_sig are derived from these same tables rather than
  ## recomputed, so the counts and the identities cannot drift apart.
  attr(out, "arm_out") <- list(bootstrap = outB, laplace = outL)
  arm_sig    <- function(tab) stats::setNames(tab$sig.q, rownames(tab))
  sig_by_arm <- list(bootstrap = arm_sig(outB), laplace = arm_sig(outL))
  n_arm <- function(v) if (all(is.na(v))) NA_integer_ else sum(v, na.rm = TRUE)
  attr(out, "sig_by_arm") <- sig_by_arm
  attr(out, "n_sig")      <- c(bootstrap = n_arm(sig_by_arm$bootstrap),
                               laplace   = n_arm(sig_by_arm$laplace))

  # decision inputs travel with the result; out is still a plain data.frame
  attr(out, "scale_model")    <- chosen
  attr(out, "scale_rule")     <- scale_rule
  attr(out, "scale_stat")     <- scale_stat
  attr(out, "v_contribution") <- list(bootstrap = vB, laplace = vL)
  attr(out, "n_laplace_fail") <- length(na.ind)
  attr(out, "n_var_reject")   <- length(var.rej)
  attr(out, "var_max")        <- var_max

  ## theta^perp, the estimated scale --------------------------------------------
  ## Both arms' scalar posterior for the shift that turns relative log fold
  ## changes into absolute ones, so the fit can be scored against an external
  ## total-abundance measurement (qPCR / flow cytometry / absolute counts) and not
  ## only by which features it calls. Definition and sign convention: see the
  ## block above sparse_perp().
  ##
  ## Row 1 of theta_perpB is the whole bootstrap story -- -mode is a scalar per
  ## draw, recycled down all D rows, so every row is identical.
  ##
  ## Row 1 of theta_perpL is ONE exact draw from that arm's marginal
  ## N(-mode[s], sigma^2[s]); the jitter is iid across taxa so no row is special,
  ## and it is taken from draws ALREADY made rather than by calling rnorm() again
  ## -- a new call would advance the random stream and silently change every
  ## result computed after this function. The exact mixture rides along in
  ## center/var, so the summary never rests on that single draw.
  ##
  ## Same retained draws as the arm's own table: the laplace columns dropped for a
  ## failed variance are dropped here too, the bootstrap arm keeps all S.
  perp_keep <- if(length(na.ind) > 0) setdiff(seq_len(S), na.ind) else seq_len(S)
  theta_perp_out <- list(
    scale_model = chosen,
    bootstrap   = list(draws = theta_perpB[1, ], center = NULL, var = NULL),
    laplace     = list(draws  = theta_perpL[1, perp_keep],
                       center = mode_estL[perp_keep],
                       var    = var_estL[perp_keep]))
  attr(out, "theta_perp")   <- theta_perp_out
  attr(out, "perp_summary") <- ssrv_perp_summary(theta_perp_out)

  return(out)
}

#sensitivity to bandwidth
check_sparse_bw_sens <- function(dat.array = resam.array, Y, X, D, N, S, bandwidth = NULL, file_path = NULL,
                                 filtered.ind = rep(FALSE, D),
                                 scale_draw = c("shared", "independent")){
  ## Threaded through to sparse_auto_sens() rather than left at its default, or
  ## a sweep asked for under one scale draw would silently fit the other.
  scale_draw <- match.arg(scale_draw)
  
  tax_sens_list <- list()
  if(is.null(bandwidth)){
    bandwidth = c("nrd0", "nrd", "SJ-ste")
  }
  
  for(i in 1:length(bandwidth)){
    bw = bandwidth[i]
    # These used to read the global `resam.array` and hardcode filtered.ind, so
    # BOTH arguments were silently ignored: passing a different dat.array or a real
    # filter had no effect, and the call failed outright if no global resam.array
    # existed. Now the function's own arguments are used.
    if(is.null(file_path)){
      out <- sparse_auto_sens(dat.array = dat.array, Y, X, D, N, S, bandwidth = bw, file_path = NULL,
                              filtered.ind = filtered.ind, scale_draw = scale_draw)
    }else{
      file_path_bw <- paste0(file_path, "sparse_mode_est_bw_", bw, ".RData")
      out <- sparse_auto_sens(dat.array = dat.array, Y, X, D, N, S, bandwidth = bw, file_path = file_path_bw,
                              filtered.ind = filtered.ind, scale_draw = scale_draw)
    }
    tax_ssrv_qval <- out[out$sig.q == "TRUE",]
    tax_ssrv_qval$category <- rownames(tax_ssrv_qval)
    colnames(tax_ssrv_qval)[4] <- "Sig"

    ## Subsetting a data.frame drops custom attributes, so the per-arm tables the
    ## fit produced would be lost here and the bandwidth sweep would again report
    ## the selected arm only. Re-attached, so sparse_arms(sparse_sens$nrd) says
    ## what each scale model called AT THAT BANDWIDTH. Attributes are invisible to
    ## rownames() and nrow(), so every existing reader is unaffected.
    ##
    ## One asymmetry to know about when reading it back: $laplace and $bootstrap
    ## are the FULL tables of D rows, while $hybrid is this already-subset frame
    ## of significant rows only (with `category` added and column 4 renamed).
    ## Compare like for like with tab[tab$sig.q %in% TRUE, ].
    attr(tax_ssrv_qval, "arm_out")     <- attr(out, "arm_out")
    attr(tax_ssrv_qval, "scale_model") <- attr(out, "scale_model")
    attr(tax_ssrv_qval, "sig_by_arm")  <- attr(out, "sig_by_arm")
    attr(tax_ssrv_qval, "n_sig")       <- attr(out, "n_sig")
    ## theta^perp too, for the same reason: sigma^2 is O(h^2), so the estimated
    ## scale is the quantity this sweep is most informative about. Without it,
    ## sparse_perp(sparse_sens$nrd) would report nothing at the one bandwidth the
    ## user actually varied.
    attr(tax_ssrv_qval, "theta_perp")   <- attr(out, "theta_perp")
    attr(tax_ssrv_qval, "perp_summary") <- attr(out, "perp_summary")

    tax_sens_list[[i]] <- tax_ssrv_qval
  }
  names(tax_sens_list) <- bandwidth
  return(tax_sens_list)
}

################################################################################
## sparse_auto_sens: sparse_auto() with the KDE bandwidth exposed.
##
## Kept deliberately in step with sparse_auto() -- same closed-form Laplace
## curvature, same variance-contribution selection rule, same guards. The ONLY
## intended difference is the `bandwidth` argument, which is threaded into every
## density() call so a sensitivity sweep varies exactly one thing.
##
## Why matching sparse_auto() matters more here than anywhere else: sigma^2 is
## O(h^2) -- exactly h^2 when every LFC coincides -- so the bandwidth IS the Laplace
## scale uncertainty. The previous numDeriv-on-approx() path reported a step-size
## artefact instead, which left this function, the bandwidth sensitivity check, the
## one least able to see the bandwidth it was varying. bw_numeric and var_estL are
## now recorded so the h^2 relationship can be read off directly.
################################################################################
## Kept in step with sparse_auto() deliberately, including the common-theta^||
## comparison -- this function exists to vary the bandwidth, and it would be
## useless for that if it also differed in the selection rule.
sparse_auto_sens <- function(dat.array, Y, x, D, N, S, bandwidth = NULL , file_path=NULL, filtered.ind = NULL,
                             scale_rule = c("width", "majority", "median", "q10", "mean"),
                             scale_draw = c("shared", "independent"),
                             var_max = 100){

  scale_rule <- match.arg(scale_rule)
  scale_draw <- match.arg(scale_draw)

  theta_perpB <- matrix(NA, nrow = D, ncol = S)
  theta_perpL <- matrix(NA, nrow = D, ncol = S)
  theta_pll <- matrix(NA, nrow = D, ncol = S)
  thetaB <- matrix(NA, nrow = D, ncol = S)
  thetaL <- matrix(NA, nrow = D, ncol = S)
  theta <- matrix(NA, nrow = D, ncol = S)

  # theta^|| for each branch, so Var(theta^||) can be subtracted from the matching
  # arm -- see sparse_auto() for the identity these feed
  theta_pllB_store <- matrix(NA, nrow = D, ncol = S)
  theta_pllL_store <- matrix(NA, nrow = D, ncol = S)

  log_wn <- matrix(NA, nrow = D, ncol = N)
  theta_pllB <- matrix(NA, nrow = D, ncol = S)
  # Both of these store -mode, i.e. theta^perp itself, NOT the mode -- see the
  # note at the same declaration in sparse_auto(). The names predate theta^perp
  # being reported and are kept so the saved sparse_mode_est_bw_*.RData files
  # stay readable by whatever already reads them.
  mode_est  <- c()   # bootstrap branch: -mode, as before
  mode_estL <- c()   # laplace branch:   -mode
  var_estL  <- c()   # laplace variance per draw -- the quantity that scales as h^2
  bw_numeric <- c()  # the bandwidth density() actually used, per draw

  na.ind <- c()
  var.rej <- c()

  if(is.null(bandwidth)){
    bandwidth = "nrd0"
  }

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
    theta_pllB_store[,s] <- theta_pll[,s]

    ### bootstrap over D
    d.index <- sample(1:D, replace = T)
    theta_pllB[,s] <- theta_pll[d.index,s]

    dens.LFC <- density(theta_pllB[,s], na.rm = T, bw = bandwidth)

    # find mode
    theta_perpB[,s] <- -dens.LFC$x[which.max(dens.LFC$y)]
    mode_est[s] <- -dens.LFC$x[which.max(dens.LFC$y)]

    # get theta using bootstrap method
    thetaB[,s] <- theta_pll[,s] + theta_perpB[,s]


    ##### laplace approximation
    theta_pll[!filtered.ind,s] = rowMeans(log_wn[,x==1]) - rowMeans(log_wn[,x==0])
    theta_pll[filtered.ind,s] = NA
    theta_pllL_store[,s] <- theta_pll[,s]

    # get density at the bandwidth under test; kde_mode_var() reads h from dens$bw,
    # so the closed-form curvature is evaluated at that same bandwidth
    dens.LFC <- density(theta_pll[,s], na.rm = T, bw = bandwidth)
    bw_numeric[s] <- dens.LFC$bw

    mv <- kde_mode_var(theta_pll[,s], dens.LFC)
    mode <- mv$mode
    theta_perpL[,s] <- -mode
    variance <- mv$variance
    mode_estL[s] <- -mode
    var_estL[s]  <- variance

    # safety net only -- at a genuine KDE maximum the closed form gives variance > 0
    if((!is.na(variance)) & variance > 0 & variance < var_max){
      ## rnorm(1), NOT rnorm(D). theta^perp is ONE number -- the log2 fold change
      ## of the total load -- shared by every taxon, so the draw is a scalar that
      ## recycles across the D rows. Drawing D independent values leaves each
      ## taxon's MARGINAL variance unchanged, which is why it went unnoticed, but
      ## it zeroes the cross-taxon correlation the scale model exists to induce
      ## (measured at D = 200: cor 0.996 shared vs 0.001 independent) and shrinks
      ## the community-wide shift the posterior represents by sqrt(D) (sd 0.801
      ## vs 0.075 at sigma = 0.8). The bootstrap branch above has always used one
      ## scalar per draw, and ssrv_perp_table() summarises this arm as a scalar
      ## mixture, so both of those now describe the model that is actually fitted.
      theta_perpL[,s] <- theta_perpL[,s] +
        rnorm(if (scale_draw == "shared") 1 else D, 0, sd = sqrt(variance))
      thetaL[,s] <- theta_pll[,s] + theta_perpL[,s]
    }else{
      thetaL[,s] <- NA
      na.ind <- c(na.ind, s)
      if(!is.na(variance) && variance >= var_max) var.rej <- c(var.rej, s)
    }
  }

  # drop failed draws from thetaL and its theta^|| store together, so the variances
  # below are taken over the same columns
  if(length(na.ind)>0){
    thetaL <- thetaL[,-na.ind, drop = FALSE]
    theta_pllL_store <- theta_pllL_store[,-na.ind, drop = FALSE]
    # separate variable, not a reassignment -- see the note in sparse_auto()
    theta_perpB_keep <- theta_perpB[,-na.ind, drop = FALSE]
  }else{
    theta_perpB_keep <- theta_perpB
  }
  if(is.null(dim(thetaL)) || ncol(thetaL) < 2){
    stop("laplace branch retained fewer than 2 draws at bandwidth '", bandwidth,
         "'; cannot compare scale models")
  }

  ### variance theta^perp contributes to theta, taxon by taxon (see sparse_auto)
  ## COMMON theta^|| -- both branches compared against the plug-in theta^||, so
  ## the statistic isolates the scale model rather than also charging each branch
  ## for the different theta^|| it is built on. Full rationale in sparse_auto().
  varPll_plug <- apply(theta_pllL_store, 1, var, na.rm = TRUE)
  vB <- apply(theta_pllL_store + theta_perpB_keep, 1, var, na.rm = TRUE) - varPll_plug
  vL <- apply(thetaL, 1, var, na.rm = TRUE) - varPll_plug

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
    majority = c(bootstrap = mean(vB >  vL, na.rm = TRUE),
                 laplace   = mean(vL >= vB, na.rm = TRUE)),
    median   = c(bootstrap = median(vB, na.rm = TRUE),
                 laplace   = median(vL, na.rm = TRUE)),
    q10      = c(bootstrap = unname(quantile(vB, 0.1, na.rm = TRUE)),
                 laplace   = unname(quantile(vL, 0.1, na.rm = TRUE))),
    mean     = c(bootstrap = mean(vB, na.rm = TRUE),
                 laplace   = mean(vL, na.rm = TRUE)))

  ##### choose the scale model contributing the larger variance
  if(scale_stat[1] > scale_stat[2]){
    theta <- thetaB
    theta_perp <- theta_perpB
    chosen <- "bootstrap"
  }else{
    theta <- thetaL
    theta_perp <- theta_perpL
    chosen <- "laplace"
  }

  if(!is.null(file_path)){
    # mode_est keeps its original meaning (bootstrap-branch modes) so existing
    # readers still work; the rest is what makes a bandwidth sweep interpretable.
    saveRDS(list(mode_est    = mode_est,
                 scale_model = chosen,
                 bandwidth   = bandwidth,
                 bw_numeric  = bw_numeric,
                 mode_estL   = mode_estL,
                 var_estL    = var_estL,
                 var_over_h2 = var_estL / bw_numeric^2,
                 scale_rule  = scale_rule,
                 scale_stat  = scale_stat,
                 n_laplace_fail = length(na.ind),
                 n_var_reject   = length(var.rej)),
            file = file_path)
  }

  ## ALL THREE results, exactly as in sparse_auto() -- the hybrid, the laplace
  ## scale model on its own and the bootstrap scale model on its own, all built
  ## by ssrv_arm_table() so a bandwidth sweep can compare arms at each bandwidth
  ## feature by feature and not only by count. `out` IS the selected arm's table.
  ids  <- rownames(Y)
  outB <- ssrv_arm_table(thetaB, ids, "bootstrap")
  outL <- ssrv_arm_table(thetaL, ids, "laplace")
  out  <- if(chosen == "bootstrap") outB else outL

  # in case a list of filtered taxa is given
  if(sum(filtered.ind) == 0){
    filter = "FALSE"
  }else{
    filter = "TRUE"
  }

  ## Every arm's FULL table travels with the result -- see sparse_auto(). Read
  ## them back with sparse_arms(out); sig_by_arm / n_sig are derived from these
  ## same tables, so counts and identities cannot drift apart.
  attr(out, "arm_out") <- list(bootstrap = outB, laplace = outL)
  arm_sig    <- function(tab) stats::setNames(tab$sig.q, rownames(tab))
  sig_by_arm <- list(bootstrap = arm_sig(outB), laplace = arm_sig(outL))
  n_arm <- function(v) if (all(is.na(v))) NA_integer_ else sum(v, na.rm = TRUE)
  attr(out, "sig_by_arm") <- sig_by_arm
  attr(out, "n_sig")      <- c(bootstrap = n_arm(sig_by_arm$bootstrap),
                               laplace   = n_arm(sig_by_arm$laplace))

  # same attributes sparse_auto() attaches, plus the bandwidth actually used
  attr(out, "scale_model")    <- chosen
  attr(out, "scale_rule")     <- scale_rule
  attr(out, "scale_stat")     <- scale_stat
  attr(out, "v_contribution") <- list(bootstrap = vB, laplace = vL)
  attr(out, "n_laplace_fail") <- length(na.ind)
  attr(out, "n_var_reject")   <- length(var.rej)
  attr(out, "var_max")        <- var_max
  attr(out, "bandwidth")      <- bandwidth
  attr(out, "bw_numeric")     <- median(bw_numeric)
  attr(out, "var_estL")       <- var_estL

  ## theta^perp at THIS bandwidth -- see sparse_auto() for what is stored and why
  ## row 1 is a valid draw. Carried here as well because sigma^2 is O(h^2): the
  ## bandwidth sweep is the sensitivity analysis for the scale interval itself,
  ## and it could not report on it while theta^perp went unrecorded.
  perp_keep <- if(length(na.ind) > 0) setdiff(seq_len(S), na.ind) else seq_len(S)
  theta_perp_out <- list(
    scale_model = chosen,
    bootstrap   = list(draws = theta_perpB[1, ], center = NULL, var = NULL),
    laplace     = list(draws  = theta_perpL[1, perp_keep],
                       center = mode_estL[perp_keep],
                       var    = var_estL[perp_keep]))
  attr(out, "theta_perp")   <- theta_perp_out
  attr(out, "perp_summary") <- ssrv_perp_summary(theta_perp_out)

  return(out)
}

sparse_auto_sens_multi <- function(dat.array, Y, x, D, N, S, bandwidth = c("nrd0", "nrd", "SJ-ste"), filtered.ind = NULL){
  
  theta_perpB1 <- matrix(NA, nrow = D, ncol = S)
  theta_perpL1 <- matrix(NA, nrow = D, ncol = S)
  theta_perpB2 <- matrix(NA, nrow = D, ncol = S)
  theta_perpL2 <- matrix(NA, nrow = D, ncol = S)
  theta_perpB3 <- matrix(NA, nrow = D, ncol = S)
  theta_perpL3 <- matrix(NA, nrow = D, ncol = S)
  theta_pll <- matrix(NA, nrow = D, ncol = S)
  thetaB1 <- matrix(NA, nrow = D, ncol = S)
  thetaL1 <- matrix(NA, nrow = D, ncol = S)
  thetaB2 <- matrix(NA, nrow = D, ncol = S)
  thetaL2 <- matrix(NA, nrow = D, ncol = S)
  thetaB3 <- matrix(NA, nrow = D, ncol = S)
  thetaL3 <- matrix(NA, nrow = D, ncol = S)
  theta1 <- matrix(NA, nrow = D, ncol = S)
  theta2 <- matrix(NA, nrow = D, ncol = S)
  theta3 <- matrix(NA, nrow = D, ncol = S)
  
  log_wn <- matrix(NA, nrow = D, ncol = N)
  theta_pllB <- matrix(NA, nrow = D, ncol = S)
  na.ind <- c()
  
  if(is.null(bandwidth)){
    bandwidth = "nrd0"
  }
  
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
    
    ### 
    dens.LFC <- density(theta_pllB[,s], na.rm = T, bw = bandwidth[1])
    theta_perpB1[,s] <- -dens.LFC$x[which.max(dens.LFC$y)]
    thetaB1[,s] <- theta_pll[,s] + theta_perpB1[,s]
    
    dens.LFC <- density(theta_pllB[,s], na.rm = T, bw = bandwidth[2])
    theta_perpB2[,s] <- -dens.LFC$x[which.max(dens.LFC$y)]
    thetaB2[,s] <- theta_pll[,s] + theta_perpB2[,s]
    
    dens.LFC <- density(theta_pllB[,s], na.rm = T, bw = bandwidth[3])
    theta_perpB3[,s] <- -dens.LFC$x[which.max(dens.LFC$y)]
    thetaB3[,s] <- theta_pll[,s] + theta_perpB3[,s]
    
    ##### laplace approximation
    theta_pll[!filtered.ind,s] = rowMeans(log_wn[,x==1]) - rowMeans(log_wn[,x==0])
    theta_pll[filtered.ind,s] = NA
    
    # get density
    dens.LFC <- density(theta_pll[,s], na.rm = T, bw = bandwidth[1])
    log_density <- function(x) {
      approx(dens.LFC$x, log(dens.LFC$y), xout = x, rule = 2)$y
    }
    mode <- optimize(f = log_density, interval = range(theta_pll[,s], na.rm = T), maximum = TRUE)$maximum
    theta_perpL1[,s] <- -mode
    hessian <- hessian(func = log_density, x = mode)
    variance <- -1 / hessian 
    if((!is.na(variance)) & variance > 0){
      theta_perpL1[,s] <- theta_perpL1[,s] + rnorm(D, 0, sd = sqrt(variance))
      thetaL1[,s] <- theta_pll[,s] + theta_perpL1[,s]
    }else{
      thetaL[,s] <- NA
      na.ind <- c(na.ind, s)
    }
    
    dens.LFC <- density(theta_pll[,s], na.rm = T, bw = bandwidth[2])
    log_density <- function(x) {
      approx(dens.LFC$x, log(dens.LFC$y), xout = x, rule = 2)$y
    }
    mode <- optimize(f = log_density, interval = range(theta_pll[,s], na.rm = T), maximum = TRUE)$maximum
    theta_perpL2[,s] <- -mode
    hessian <- hessian(func = log_density, x = mode)
    variance <- -1 / hessian 
    if((!is.na(variance)) & variance > 0){
      theta_perpL2[,s] <- theta_perpL2[,s] + rnorm(D, 0, sd = sqrt(variance))
      thetaL2[,s] <- theta_pll[,s] + theta_perpL2[,s]
    }else{
      thetaL2[,s] <- NA
      na.ind <- c(na.ind, s)
    }
    
    dens.LFC <- density(theta_pll[,s], na.rm = T, bw = bandwidth[3])
    log_density <- function(x) {
      approx(dens.LFC$x, log(dens.LFC$y), xout = x, rule = 2)$y
    }
    mode <- optimize(f = log_density, interval = range(theta_pll[,s], na.rm = T), maximum = TRUE)$maximum
    theta_perpL3[,s] <- -mode
    hessian <- hessian(func = log_density, x = mode)
    variance <- -1 / hessian 
    if((!is.na(variance)) & variance > 0){
      theta_perpL3[,s] <- theta_perpL3[,s] + rnorm(D, 0, sd = sqrt(variance))
      thetaL3[,s] <- theta_pll[,s] + theta_perpL3[,s]
    }else{
      thetaL3[,s] <- NA
      na.ind <- c(na.ind, s)
    }
  } 
  
  if(length(na.ind)>0){
    theta1 <- theta1[,-na.ind]
  }
  
  out_bws <- list()
  for(b in 1:3){
    if(b==1){
      thetaB = thetaB1
      theta_perpB = theta_perpB1
      thetaL = thetaL1
      theta_perpL = theta_perpL1
    }else if(b==2){
      thetaB = thetaB2
      theta_perpB = theta_perpB2
      thetaL = thetaL2
      theta_perpL = theta_perpL2
    }else if(b==3){
      thetaB = thetaB3
      theta_perpB = theta_perpB3
      thetaL = thetaL3
      theta_perpL = theta_perpL3
    }
    ### compare uncertainty added from Bootstrap method and Laplace approximation
    # get average standard deviation of each rows (sd of theta across samples)
    uncertaintyB_interval <- t(apply(thetaB, MARGIN = 1, function(x) quantile(x,c(0.05, 0.95), na.rm = TRUE)))
    uncertaintyB <- mean(apply(uncertaintyB_interval, MARGIN = 1, function(x) {x[2]-x[1]}), na.rm = TRUE)
    uncertaintyL_interval <- t(apply(thetaL, MARGIN = 1, function(x) quantile(x,c(0.05, 0.95), na.rm = TRUE)))
    uncertaintyL <- mean(apply(uncertaintyL_interval, MARGIN = 1, function(x) {x[2]-x[1]}), na.rm = TRUE)
    
    ##### choose the method with bigger uncertainty as final output
    if(uncertaintyB > uncertaintyL){
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
    out_bws[[b]] <- out
  }
  names(out_bws) <- bandwidth
  return(out_bws)
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

# # For simulation ; difference is just taking metadata_df as an input
# sim_sparse_auto <- function(dat.array, Y, x, D, N, S, output_dir, metadata_df = NULL, filtered.ind = NULL){
#   
#   theta_pll <- matrix(NA, nrow = D, ncol = S)
#   theta_perpB <- matrix(NA, nrow = D, ncol = S)
#   theta_perpL <- matrix(NA, nrow = D, ncol = S)
#   thetaB <- matrix(NA, nrow = D, ncol = S)
#   thetaL <- matrix(NA, nrow = D, ncol = S)
#   theta <- matrix(NA, nrow = D, ncol = S)
#   
#   log_wn <- matrix(NA, nrow = D, ncol = N)
#   theta_pllB <- matrix(NA, nrow = D, ncol = S)
#   
#   na.ind <- c()
#   
#   for(s in 1:S){
#     
#     log_wn <- dat.array[s, !filtered.ind, ]
#     
#     ## bootstrap over N
#     n.index <- sample(1:N, replace = T)
#     log_wn_bootstrap <- log_wn[,n.index]
#     treat.ind2 <- x[n.index] == 1
#     cont.ind2 <- x[n.index] == 0 
#     
#     treat_bootstrap <- if (is.matrix(log_wn_bootstrap[,treat.ind2])) rowMeans(log_wn_bootstrap[,treat.ind2]) else log_wn_bootstrap[,treat.ind2]
#     cont_bootstrap <- if (is.matrix(log_wn_bootstrap[,cont.ind2])) rowMeans(log_wn_bootstrap[,cont.ind2]) else log_wn_bootstrap[,cont.ind2]
#     theta_pll[!filtered.ind,s] = treat_bootstrap - cont_bootstrap
#     theta_pll[filtered.ind,s] = NA
#     
#     ## bootstrap over D
#     d.index <- sample(1:D, replace = T)
#     theta_pllB[,s] <- theta_pll[d.index,s]
#     
#     dens.LFC <- density(theta_pllB[,s], na.rm = T)
#     theta_perpB[,s] <- -dens.LFC$x[which.max(dens.LFC$y)]
#     
#     thetaB[,s] <- theta_pll[,s] + theta_perpB[,s]
#     
#     ## laplace approximation
#     
#     theta_pll[!filtered.ind,s] = rowMeans(log_wn[,x==1]) - rowMeans(log_wn[,x==0])
#     theta_pll[filtered.ind,s] = NA
#     
#     dens.LFC <- density(theta_pll[,s], na.rm = T)
#     
#     log_density <- function(x) {
#       approx(dens.LFC$x, log(dens.LFC$y), xout = x, rule = 2)$y
#     }
#     
#     mode <- optimize(f = log_density, interval = range(theta_pll[,s], na.rm = T), maximum = TRUE)$maximum
#     theta_perpL[,s] <- -mode
#     hessian <- hessian(func = log_density, x = mode)
#     variance <- -1 / hessian 
#     
#     
#     # in case hessian have numerical issue
#     if((!is.na(variance)) & variance > 0){
#       theta_perpL[,s] <- theta_perpL[,s] + rnorm(D, 0, sd = sqrt(variance))
#       thetaL[,s] <- theta_pll[,s] + theta_perpL[,s]
#     }else{
#       thetaL[,s] <- NA
#       na.ind <- c(na.ind, s)
#     }
#   }  
#   
#   if(length(na.ind)>0){
#     theta <- theta[,-na.ind]
#   }
#   
#   uncertaintyB <- mean(apply(thetaB, MARGIN = 1, function(x) sd(x, na.rm = T)), na.rm = T)
#   uncertaintyL <- mean(apply(thetaL, MARGIN = 1, function(x) sd(x, na.rm = T)), na.rm = T)
#   
#   if(uncertaintyB > uncertaintyL){
#     theta <- thetaB
#   }else{
#     theta <- thetaL
#   }
#   
#   # calculate posterior predictive p-value
#   p.val <- apply(theta, 1, FUN = function(x){min(sum(x < 0), sum(x > 0)) * 2 / ncol(theta) })
#   q.val <- p.adjust(p.val, method = "BH")
#   
#   out <- apply(theta, 1, FUN = function(x){quantile(x, c(0.025, 0.975), na.rm = T)})
#   mean <- apply(theta, 1, FUN = function(x){mean(x)})
#   out <- cbind(as.data.frame(t(out)), mean)
#   sig <- (out$`2.5%` * out$`97.5%` > 0)
#   sig.q <- (q.val < 0.05)
#   out <- cbind(out, sig, sig.q)
#   rownames(out) <- rownames(Y)
#   
#   if(!is.null(metadata_df)){
#     #spiked.ind <- as.numeric(str_remove_all(metadata_df$feature_spiked,"Feature"))
#     matched_indices <- match(as.character(metadata_df$feature_spiked), rownames(Y), nomatch=0)
#     #spiked.ind <- rownames(Y) %in% metadata_df$feature_spiked
#     effect_sizes <- rep(0, D)
#     effect_sizes[matched_indices] <- metadata_df$effect_size
#     
#     out <- cbind(out, effect_sizes)
#   }
#   
#   if(sum(filtered.ind) == 0){
#     filter = "FALSE"
#   }else{
#     filter = "TRUE"
#   }
#   
#   #return(out)
#   write.csv(out, paste0(output_dir, "/sparse_auto_filter",filter,".csv"))
#   
#   write.csv(outL, paste0(output_dir, "/sparse_auto_filter",filter,".csv"))
#   write.csv(outB, paste0(output_dir, "/sparse_auto_filter",filter,".csv"))
# }

sim_sparse_auto <- function(dat.array, Y, x, D, N, S, output_dir, metadata_df = NULL, filtered.ind = NULL){
  
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
    
    # mode <- optimize(f = log_density, interval = range(theta_pll[,s], na.rm = T), maximum = TRUE)$maximum
    # #mode <- dens.LFC$x[which.max(dens.LFC$y)]
    # theta_perpL[,s] <- -mode
    # hessian <- hessian(func = log_density, x = mode)
    # variance <- -1 / hessian 
    
    # different way of finding mode
    # first find mode using KDE
    mode <- dens.LFC$x[which.max(dens.LFC$y)]
    # but KDE can oversmooth density - can flush away dirichlet noise
    # use optimization - which is more sensitive, but can stuck at local mode, not global mode
    # hence, use optimization but using interval that include the mode found via KDE
    mode <- optimize(f = log_density, interval = range(c(mode-0.5,mode+0.5), na.rm = T), maximum = TRUE)$maximum
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
  
  uncertaintyB <- mean(apply(thetaB, MARGIN = 1, function(x) sd(x, na.rm = T)), na.rm = T)
  uncertaintyL <- mean(apply(thetaL, MARGIN = 1, function(x) sd(x, na.rm = T)), na.rm = T)
  
  if(uncertaintyB > uncertaintyL){
    print("bootstrap")
    theta <- thetaB
  }else{
    theta <- thetaL
  }
  
  # calculate posterior predictive p-value
  p.val <- apply(theta, 1, FUN = function(x){min(sum(x < 0), sum(x > 0)) * 2 / ncol(theta) })
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
  write.csv(out, paste0(output_dir, "/sparse_auto_filter",filter,".csv"))
  
  
  ################
  # calculate posterior predictive p-value
  p.val <- apply(thetaL, 1, FUN = function(x){min(sum(x < 0), sum(x > 0)) * 2 / ncol(thetaL) })
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
  write.csv(out, paste0(output_dir, "/sparse_L_filter",filter,".csv"))
  
  
  ################
  # calculate posterior predictive p-value
  p.val <- apply(thetaB, 1, FUN = function(x){min(sum(x < 0), sum(x > 0)) * 2 / ncol(thetaB) })
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
  #return(out)
  write.csv(out, paste0(output_dir, "/sparse_B_filter",filter,".csv"))
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

