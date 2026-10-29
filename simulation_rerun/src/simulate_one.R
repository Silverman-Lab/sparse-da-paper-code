################################################################################
## simulate_one.R -- the SparseDOSSA2 data-generating block, defined ONCE.
##
## Ported from simulation/src/05052025/sparse_analysis_trycatch.R:35-137. The
## Phase 0 pilot, Study A's worker and Study B's worker all call this, so they
## generate identically by construction rather than by inspection. The existing
## pipeline had three near-copies of this block (the worker, the oracle worker,
## and scale_choice/run_scale_choice.R) which had already drifted apart.
##
## Faithful to the original in every respect that touches the RNG:
##   * three set.seed(seed) calls at the same three points;
##   * SparseDOSSA2 called TWICE -- once at median_read_depth = 50000 only to
##     obtain taxa names, then again with the spike metadata at the real lambda;
##   * lfcs drawn between the two calls, in the same order and with the same
##     rbeta parameterisation (ratio1_4: 1/5 decreasing, 4/5 increasing).
## Reordering any of it changes every downstream draw.
################################################################################

## effect_size is a shift in the NATURAL-LOG mean of absolute abundance, so the
## true log2 fold change is effect_size / log(2). Magnitudes land in [1, 6]
## natural-log units, i.e. [1.44, 8.66] in log2.
##
## `conc` controls how SHARPLY PEAKED the effect-size distribution is, using
## Beta(conc, 3*conc). That family has mean 1/4 for every conc, so
##
##     |lfc| = 5 * Beta(conc, 3*conc) + 1   has mean 2.25 at every conc
##
## and only the SPREAD changes -- which is what isolates "sharper effect sizes"
## from "smaller effect sizes". Var(Beta(a,3a)) = (3/16)/(4a+1), so
##
##   conc      1     3     10    30    Inf
##   sd|lfc|   0.97  0.60  0.34  0.20  0    (natural-log units)
##
## conc = 1 is Beta(1,3): the original distribution, drawn with the same single
## rbeta() call per sign group, so it is bit-identical to the previous pipeline.
## conc = Inf is a point mass at the common mean.
##
## Why this might matter: the Laplace variance is read off the curvature of the
## log-KDE at the mode, and nrd0 = 0.9 * min(sd, IQR/1.34) * D^(-1/5). A wide
## effect distribution inflates `sd`. Whether that reaches the bandwidth depends
## on which branch of the min() binds -- at low pi the nulls dominate both
## quartiles, so IQR/1.34 is the smaller and the effect spread is discarded.
## This lever measures it instead of assuming.
##
## ---------------------------------------------------------------------------
## prop_pos -- the SIGN SPLIT of the spiked effects. Added 2026-09-16.
## ---------------------------------------------------------------------------
## 4/5 is the design of the main manuscript: 80% of the spiked taxa increase and
## 20% decrease. The supplementary scenario is the MIRROR, prop_pos = 1/5, and
## nothing else about the design changes -- same |lfc| distribution, same number
## of spiked taxa, same everything downstream.
##
## Why it is not a cosmetic switch: the majority sign is what the compositional
## methods see as the "background". With 80% increasing, the reference (the
## geometric mean / the CLR centre) is dragged UP by the signal, so the nulls
## look like they went DOWN and the sign of the induced compositional bias is
## negative; at 1/5 it is the other way round. Sparse SSRV's KDE mode, which
## tracks the null cluster, therefore sits on the opposite side of the true load
## shift. So the two scenarios are not related by any relabelling of the output.
##
## The counts are ceiling() on BOTH sides, as the original was, so k_neg + k_pos
## can exceed num.spiked by one at odd values -- that is the published behaviour
## and it is the SAME total under either split (verified over every (pi, D) in
## the pilot and main grids), so the two scenarios spike the same number of taxa.
##
## Bit-identity of the default: ceiling(num.spiked * (1 - 4/5)) and the original
## ceiling(num.spiked * 1/5) agree at every setting in either grid, even though
## 1 - 4/5 and 1/5 are not the same double (they differ by 1 ulp, and ceiling of
## an integer minus an ulp is that integer). tests/test_edits_additive.R asserts
## it, and the default path draws in the same order -- negatives first -- so a
## replicate run at prop_pos = 4/5 is bit-identical to one run before this
## argument existed.
make_lfcs <- function(num.spiked, conc = 1, prop_pos = 4/5) {
  k_neg <- ceiling(num.spiked * (1 - prop_pos))
  k_pos <- ceiling(num.spiked * prop_pos)
  draw <- function(k) {
    if (is.infinite(conc)) rep(0.25, k) else rbeta(k, shape1 = conc, shape2 = 3 * conc)
  }
  c(-draw(k_neg) * 5 - 1,
     draw(k_pos) * 5 + 1)
}

## The template surgery: draw per-feature parameters for d features from the
## Stool fit, then switch OFF zero inflation (pi0 = 0) so the zeros that appear
## are sampling zeros rather than structural ones.
##
## NOTE the seed. generate_featureParam() consumes RNG, so the template must be
## drawn inside the seeded stream or a replicate is not reproducible from its
## logged seed.
## sigma_scale is the INDEPENDENT (per-taxon) heterogeneity lever: it multiplies
## the template's per-feature sigma. sigma_scale = 1 is the original template and
## consumes no extra RNG, so every non-heterogeneity result is unchanged.
make_template <- function(d, sigma_scale = 1) {
  Stool <- Stool_mod <- SparseDOSSA2:::Stool
  feature_param_template <- cbind("pi0"   = Stool$EM_fit$fit$pi0,
                                  "mu"    = Stool$EM_fit$fit$mu,
                                  "sigma" = Stool$EM_fit$fit$sigma)
  feature_param <- SparseDOSSA2:::generate_featureParam(
    F_fit = Stool$F_fit,
    param_original = feature_param_template,
    n_feature = d)

  Stool_mod$EM_fit$fit$pi0   <- 0
  Stool_mod$EM_fit$fit$mu    <- feature_param[, "mu"]
  Stool_mod$EM_fit$fit$sigma <- feature_param[, "sigma"] * sigma_scale
  Stool_mod
}

## SHARED (compositional) heterogeneity: a per-sample evenness tilt. Raising a
## composition to a power changes its geometric mean directly, so this moves the
## component COMMON across taxa -- the one the KDE mode tracks and cancels, and
## therefore the one that drives Cov(theta^||, theta^perp) negative.
##
## Deliberately not a rank-one tilt with a ~ rnorm(D): mean(a) ~ 0 makes that
## mostly INDEPENDENT across taxa, which is what sigma_scale already varies.
##
## Sequencing depth is preserved exactly, so this axis cannot smuggle in a lambda
## effect. tau = 0 returns Y untouched and draws nothing.
##
## CAVEAT: the tilt modifies OBSERVED counts only; sim$simulated_matrices$a_spiked
## is unchanged, so theta_perp_true no longer corresponds exactly to the tilted
## data. That is fine for scale-model choice (which never uses it) but tau must
## NOT be carried into Study B, where theta_perp_true is the reference quantity.
tilt <- function(Y, tau, seed) {
  if (tau <= 0) return(Y)
  N <- ncol(Y); dep <- colSums(Y)
  set.seed(seed)
  u <- rnorm(N)
  p <- sweep(Y + 0.5, 2, colSums(Y + 0.5), "/")
  q <- sapply(seq_len(N), function(k) { w <- p[, k]^(1 + tau * u[k]); w / sum(w) })
  Ynew <- sapply(seq_len(N), function(k) rmultinom(1, dep[k], q[, k]))
  dimnames(Ynew) <- dimnames(Y)
  Ynew
}

## Pooled within-group SD of the per-sample log geometric mean. Isolates the
## SHARED component of between-sample variation, which is the statistic
## scale_model_report.Rmd uses to place the simulated cells against the real
## datasets (Jin 0.247, McMurrough 0.285, Stammler 0.658). Reported so the
## heterogeneity sweeps can be read on the same axis as the real studies rather
## than in units of an arbitrary multiplier.
sd_between_sample <- function(Y, X) {
  g  <- colMeans(log2(sweep(Y + 0.5, 2, colSums(Y + 0.5), "/")))
  n1 <- sum(X == 1); n0 <- sum(X == 0)
  sqrt(((n1 - 1) * sd(g[X == 1])^2 + (n0 - 1) * sd(g[X == 0])^2) / (n1 + n0 - 2))
}

## The prevalence filter, and the aggregated `other` row.
##
## other_row = "legacy" reproduces the HPC pipeline EXACTLY, including a real
## indexing bug at sparse_analysis_trycatch.R:125:
##
##     Y     <- sim$simulated_data[rownames(sim$simulated_data) %in% keep_taxa, ]
##     other <- colSums(sim$simulated_data[rownames(Y) %in% other_taxa, ])
##
## The second line indexes with rownames(Y) -- already the KEPT taxa -- so the
## logical is all FALSE, and being shorter than nrow(sim$simulated_data) it is
## recycled. Zero rows are selected and `other` is a vector of zeros. Every
## published number therefore has an all-zero `other` row: the filtered-out reads
## are discarded rather than aggregated, and SSRV's KDE mode estimation sees one
## extra taxon whose theta^|| is pure Dirichlet noise.
##
## other_row = "fixed" aggregates the filtered-out taxa as intended, matching
## apply_shared_filter() in simulation/src/oracle/oracle_arms.R:127.
##
## Left as a switch because it changes every number in the study. Whichever is
## used has to be stated in the write-up.
##
## ---------------------------------------------------------------------------
## filter_type -- WHICH rule decides the feature set. 2026-08-29, reverted
## 2026-09-01.
## ---------------------------------------------------------------------------
## "mean" (the DEFAULT since 2026-09-01): rowMeans(counts) >= min_mean,
##     min_mean = 5. This is what the real-data scripts do again -- Prasse,
##     Stammler, Zemb and selex all filter this way (Vandeputte and Jin keep
##     their source papers' rules) -- and the simulation exists to say something
##     about those analyses, so it filters the way they do.
##
##     It is an ABUNDANCE filter and implies NO prevalence floor: a taxon with all
##     of its reads in one sample passes if that count is large enough. That is
##     precisely why prv_cut / prev.cut must be 0 rather than 0.10 in
##     run_methods_rerun.R -- at 0.10 they would cut features the script kept.
##
## "prevalence": rowSums(counts > 0)/n >= min_prev, min_prev = 0.10.
##     The rule used between 2026-08-29 and 2026-09-01, kept selectable so those
##     results stay regenerable. Its attraction was that it is written the way
##     LinDA writes it (`rowSums(Y > 0)/n >= prev.cut`) and ANCOM-BC writes it
##     (`prevalence >= prv_cut`), making the script's rule and the packages' own
##     arguments literally one comparison. Those arguments are 0 either way, so
##     whichever rule is selected is the ONLY filter any method sees.
##
## THE DEFAULT MATTERS: test_edits_additive.R calls simulate_one() with
## min_mean = FILTER_MIN_MEAN and no filter_type, so it took the default while
## passing a mean threshold. Under a "prevalence" default that combination was
## silently incoherent.
##
## The two are NOT close. Measured over the Study A grid, 3 seeds per cell:
##
##     cell                D_kept mean>=5    D_kept prev>=0.10
##     baseline (5e5)          234                299.7
##     lambda = 5e3             74                184
##     lambda = 1e6            257                300
##     D = 1000                683                994
##
## At the baseline the prevalence rule keeps 299.7 of 300, so the study sits close
## to the min_mean = 0 arm that config.R's nofilter check was written to probe:
## rare taxa, whose theta^|| is hugely dispersed (rdirichlet on near-zero counts),
## stay in the LFC density and widen nrd0. Since the Laplace variance is read off
## the curvature of the log-KDE at the mode, expect the Laplace arm to come out
## MORE conservative than under min_mean = 5, most visibly at low pi.
##
## Both rules are depth-dependent, so D_kept still moves with lambda and that
## confound is still recorded rather than assumed -- see D_kept below.
apply_filter <- function(counts, min_mean = 5, other_row = c("legacy", "fixed", "drop"),
                         filter_type = c("mean", "prevalence"), min_prev = 0.10) {
  other_row   <- match.arg(other_row)
  filter_type <- match.arg(filter_type)

  ## One statistic and one threshold, so `keep` and `other` are exact complements
  ## by construction. Writing them as two independent which() calls -- which is
  ## what this used to do -- lets a taxon fall in both or neither if the
  ## comparison is ever changed on one line and not the other.
  stat <- if (filter_type == "prevalence") rowSums(counts > 0) / ncol(counts)
          else                             rowMeans(counts)
  cut  <- if (filter_type == "prevalence") min_prev else min_mean

  keep       <- stat >= cut
  keep_taxa  <- names(which(keep))
  other_taxa <- names(which(!keep))

  Y <- counts[rownames(counts) %in% keep_taxa, , drop = FALSE]

  other <- if (other_row == "legacy") {
    colSums(counts[rownames(Y) %in% other_taxa, , drop = FALSE])   # bug reproduced
  } else if (other_row == "drop") {
    NULL                                                            # no row at all
  } else {
    colSums(counts[rownames(counts) %in% other_taxa, , drop = FALSE])
  }

  ## No `other` row when nothing was filtered out (min_mean = 0). Appending one
  ## anyway would add an all-zero taxon, making D_used = D + 1 and handing the
  ## KDE a pure-noise feature -- which is precisely the artefact the unfiltered
  ## run is meant to avoid. scale_choice/run_scale_choice.R:201 keeps
  ## D_used == D for the same reason.
  ## other_row = "drop" (added 2026-09-02): the filtered-out taxa are DISCARDED
  ## rather than aggregated. This is neither of the two existing options -- the
  ## published pipeline appended an all-zero row (its recycled-logical bug), and
  ## "fixed" appends their true colSums -- so it is the third arm of that
  ## comparison, and the only one where the analysed table is exactly the
  ## surviving taxa. D_used = length(keep_taxa), not +1.
  if (length(other_taxa) > 0 && !is.null(other)) {
    Y <- rbind(Y, other)
    rownames(Y)[nrow(Y)] <- "other"
  }
  list(Y = Y, keep_taxa = keep_taxa, other_taxa = other_taxa)
}

## True between-group log2 load contrast -- the scalar separating the
## compositional estimand from the absolute one, and what SSRV's theta^perp
## estimates. Same definition as true_theta_perp() in oracle_arms.R:32.
true_theta_perp <- function(load, x) {
  mean(log2(load[x == 1])) - mean(log2(load[x == 0]))
}

## n, d, lambda, pi   the parameter set
## seed               replicate seed; everything downstream is deterministic in it
## other_row          see apply_filter()
##
## Returns the full sim object (so a_spiked is available to Study B), the
## filtered count matrix Y with its `other` row, the group labels, and the truth
## tables. Writes nothing -- callers decide what to persist.
## sigma_scale / tau default to the no-op values, so a call without them is
## bit-identical to the original pipeline.
## prop_pos defaults to the main manuscript's 4/5 (80% of spiked taxa up), so a
## call without it is bit-identical to the pipeline before the argument existed.
## 1/5 is the supplementary mirror scenario. See make_lfcs().
simulate_one <- function(n, d, lambda, pi, seed,
                         other_row = c("legacy", "fixed", "drop"), min_mean = 5,
                         filter_type = c("mean", "prevalence"), min_prev = 0.10,
                         sigma_scale = 1, tau = 0, effect_conc = 1,
                         prop_pos = 4/5) {

  filter_type <- match.arg(filter_type)
  other_row   <- match.arg(other_row)

  seed <- as.integer(seed)

  set.seed(seed)
  Stool_mod <- make_template(d, sigma_scale = sigma_scale)

  set.seed(seed)
  # first call: taxa names only, hence the fixed depth
  sim0 <- SparseDOSSA2(template = Stool_mod, n_sample = n, n_feature = d,
                       new_features = TRUE, median_read_depth = 50000)

  num.spiked <- pi * d
  lfcs <- make_lfcs(num.spiked, conc = effect_conc, prop_pos = prop_pos)
  subset_features <- sample(rownames(sim0$params$feature_param), length(lfcs))

  mat_metadata <- as.matrix(data.frame(disease = c(rep(0, n/2), rep(1, n/2))))
  metadata_df  <- data.frame(
    metadata_datum      = rep(1, length(subset_features)),
    feature_spiked      = subset_features,
    associated_property = rep("abundance", length(subset_features)),
    effect_size         = lfcs,
    stringsAsFactors    = FALSE)

  set.seed(seed)
  sim <- SparseDOSSA2(template = Stool_mod, n_sample = n, n_feature = d,
                      new_features = TRUE,
                      metadata_matrix = mat_metadata,
                      spike_metadata  = metadata_df,
                      median_read_depth = lambda)

  ## Shared-heterogeneity tilt, applied to the counts BEFORE filtering so the
  ## filter sees the data the analysis would see. Depth-preserving, and a no-op
  ## at tau = 0. The offset keeps the tilt's RNG stream clear of the data stream.
  counts <- tilt(sim$simulated_data, tau, seed + 500000L)

  flt <- apply_filter(counts, min_mean = min_mean, other_row = other_row,
                      filter_type = filter_type, min_prev = min_prev)
  Y   <- flt$Y
  x   <- mat_metadata

  filtered_metadata_df <- metadata_df[metadata_df$feature_spiked %in% rownames(Y), ]

  load <- colSums(sim$simulated_matrices$a_spiked)

  list(
    sim   = sim,
    Y     = Y,
    x     = x,
    load  = load,
    theta_perp_true = true_theta_perp(load, as.numeric(x)),
    metadata_df          = metadata_df,
    filtered_metadata_df = filtered_metadata_df,
    keep_taxa  = flt$keep_taxa,
    other_taxa = flt$other_taxa,
    # D_kept moves with lambda because the filter is depth-dependent -- record it
    # so that confound stays visible instead of silently inverting a sweep.
    D_kept = nrow(Y),
    sd_between_sample = sd_between_sample(Y, as.numeric(x)),
    ## spread of the realised |lfc| -- the thing effect_conc actually controls,
    ## recorded so the sweep can be read in natural-log units rather than in
    ## units of an arbitrary concentration parameter
    lfc_sd = sd(abs(lfcs)),
    ## REALISED sign split, not the requested one: prop_pos is a target that
    ## ceiling() rounds, so at small num.spiked the two differ (pi = 0.02, D = 300
    ## gives 5/7 = 0.714, not 0.8). Recorded per replicate so the scenario a
    ## results tree was produced under is readable from the tree itself.
    n_pos = sum(lfcs > 0),
    n_neg = sum(lfcs < 0),
    prop_pos_realised = mean(lfcs > 0),
    ## filter_type / min_prev recorded alongside min_mean: a result directory has
    ## to say which rule produced it, or the two filters' numbers are
    ## indistinguishable once they are sitting in the same results tree.
    params = list(n = n, d = d, lambda = lambda, pi = pi, seed = seed,
                  other_row = other_row, min_mean = min_mean,
                  filter_type = filter_type, min_prev = min_prev,
                  sigma_scale = sigma_scale, tau = tau,
                  effect_conc = effect_conc, prop_pos = prop_pos)
  )
}
