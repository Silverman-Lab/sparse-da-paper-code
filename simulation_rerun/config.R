################################################################################
## config.R -- single source of truth for simulation_rerun.
##
## Every script sources this FIRST. Nothing else defines a path, a grid, a
## replicate count or a core count. The original pipeline had the run tag
## hardcoded in four places that all had to agree (foo.R's RESULT_TAG, the
## src/<tag>/ script path, the MIMIX output tree, and each aggregator's string);
## that class of drift is not possible here.
################################################################################

options(scipen = 999)

## ---------------------------------------------------------------------------
## Paths. Located by walking up until sparse_SSRV.R appears, so every script runs
## from anywhere without a hardcoded absolute path -- the laptop equivalent of
## the /storage/... roots the HPC scripts assume.
## ---------------------------------------------------------------------------
find_proj <- function(start = getwd()) {
  p <- normalizePath(start, mustWork = FALSE)
  for (i in 1:6) {
    if (file.exists(file.path(p, "sparse_SSRV.R"))) return(p)
    p <- dirname(p)
  }
  stop("could not locate sparse_SSRV.R above ", start)
}

PROJ    <- find_proj()
HERE    <- file.path(PROJ, "simulation_rerun")
OLD_SIM <- file.path(PROJ, "simulation")          # read-only: sourced, never written
RES     <- file.path(HERE, "results")

RES_PILOT  <- file.path(RES, "pilot")
RES_MAIN   <- file.path(RES, "main")
RES_ORACLE <- file.path(RES, "oracle")
RES_FIG    <- file.path(RES, "figures")

for (d in c(RES, RES_PILOT, RES_MAIN, RES_ORACLE, RES_FIG,
            file.path(HERE, "params"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ---------------------------------------------------------------------------
## Machine
## ---------------------------------------------------------------------------
## 8 physical cores; one is left for the OS so an interactive session stays
## responsive during a multi-hour run. Each worker holds its own S x D x N
## resampling array, so cores are bounded by RAM, not by CPU: at S = 2000,
## D = 800, N = 20 that is 2000*800*20*8 bytes = 256 MB per worker.
N_CORES <- max(1L, parallel::detectCores() - 1L)

## ---------------------------------------------------------------------------
## Design constants shared by Phase 0, Study A and Study B
## ---------------------------------------------------------------------------
S_DRAWS <- 2000    # Dirichlet posterior draws inside ONE fit
PVAL    <- 0.05
## Which rule decides the feature set, and its threshold. See apply_filter() in
## src/simulate_one.R. The difference between the two rules is LARGE -- D_kept
## 234 under "mean" against ~300 under "prevalence" at the baseline -- so
## FILTER_TYPE must be stated in the write-up alongside OTHER_ROW, and any figure
## carried over from a run at the other setting is not comparable.
##
## SET BACK TO "mean" ON 2026-09-01, matching the real-data revert: Prasse,
## Stammler, Zemb and selex all decide their feature set with
## rowMeans(counts) >= 5 again (Vandeputte and Jin keep their papers' own rules).
## The simulation exists to say something about those analyses, so it filters the
## way they do.
##
## The argument for "prevalence" was that it is written the same way the four
## filtering packages write theirs, so the script's filter and their
## prv_cut / prev.cut arguments were one comparison instead of two that happen to
## agree. That convenience is gone, and it is replaced by something stronger:
## those arguments are all 0 (see the count-input-policy block in
## src/run_methods_rerun.R), so this remains the ONLY filter in the pipeline no
## matter which rule it implements. A mean rule implies no prevalence floor at
## all, which is exactly why they must stay 0 rather than 0.10.
FILTER_TYPE     <- "mean"
FILTER_MIN_MEAN <- 5      # rowMeans(counts) >= 5   <- IN USE
FILTER_MIN_PREV <- 0.10   # rowSums(counts > 0)/n >= 0.10, only when FILTER_TYPE = "prevalence"

## The threshold the optional 6th param_string field carries, read according to
## FILTER_TYPE. Defined once so config, the worker and the grids cannot disagree
## about which of the two constants above is in play.
FILTER_CUT <- if (FILTER_TYPE == "prevalence") FILTER_MIN_PREV else FILTER_MIN_MEAN

## How the aggregated `other` row is built. See simulate_one().
##   "legacy" reproduces the HPC pipeline exactly, INCLUDING its indexing bug,
##            which makes the row all zeros.
##   "fixed"  aggregates the filtered-out taxa as intended.
## Kept as a switch rather than a silent correction because it changes every
## number in the study; whichever is used must be stated in the write-up.
##
## Set to "fixed": Phase 0 measured both on the same seeds and the scale-model
## choice is untouched by it (max |dP(bootstrap)| = 0, max |d margin| = 0.011,
## d D_kept = 0), so the correct version costs nothing there. It also makes the
## main and oracle paths agree, which Study B needs -- it derives its shared
## filter from apply_shared_filter(), which never had the bug.
## Overridable from the environment so the two settings can be compared without
## editing this file between runs -- the same Sys.getenv idiom theta_perp_check.R
## and scale_draw_check.R already use. The DEFAULT is unchanged.
## One results tree is still one rule: probe_grid() below puts OTHER_ROW into the
## tag, so a "legacy" run cannot resume into a "fixed" tree.
OTHER_ROW <- Sys.getenv("SIM_OTHER_ROW", "fixed")
stopifnot(OTHER_ROW %in% c("fixed", "legacy", "drop"))

## Baseline history, all measured rather than assumed:
##   N=20  D=50  lam=1e4  pi=0.2   first pass         -> archive_pi0.2/
##   N=20  D=50  lam=1e4  pi=0.1   pi lowered         -> archive_D50_pi0.1/
##   N=20  D=100 lam=2e4  pi=0.1   power fix
##   N=40  D=200 lam=5e5  pi=0.1   effect-size probe
##   N=100 D=300 lam=7.5e4 pi=0.2  Phase 0c pass 1    -> ssrv_N100_D300_*/
##   N=50  D=300 lam=7.5e4 pi=0.2  Phase 0c pass 2    -> ssrv_N50_D300_*/
##   N=100 D=300 lam=5e5  pi=0.2   Study A, 09012026 tree
##   N=100 D=300 lam=3e6  pi=0.2   CURRENT -- Study A, 09022026 tree
##
## N, D and pi are the ORIGINAL published design. 60 taxa are spiked.
##
## LAMBDA RAISED 5e5 -> 3e6 ON 2026-09-02, to put the simulation on the depth of
## the template it is generated from. Measured:
##
##   SparseDOSSA2:::Stool  depth_fit = (mu_depth 14.8764, sigma_depth 0.8198)
##                         -> median depth 2,888,909 reads over 332 features
##                         -> 8,702 reads/taxon; fitted on 239 samples
##
## At 3e6 over D = 300 the simulation runs 10,000 reads/taxon, against the
## template's 8,702. The previous 5e5 was 1,667 and the AOAS published run 7.5e5
## was 2,500 -- both a quarter to a fifth of the template's own depth, i.e. the
## generator was being asked for data far shallower than the fit it came from.
##
## The real studies sit the other side: 24-689 reads/taxon (Prasse raw 24, Jin
## 50, Zemb 62, selex 152, Vandeputte 375, Stammler 689), so no single lambda is
## simultaneously template-faithful and study-faithful. The overLambda sweep is
## what covers that span; the BASELINE is pinned to the template.
##
## ACCEPTED CONSEQUENCE: rowMeans(counts) >= 5 stops binding at this depth.
## D_kept measured 72.5 at lam=5e3, 237.1 at 5e5, 262.8 at 1e6 -> ~299 of 300 at
## 3e6. The filter rule is kept identical to the real-data analysis anyway, and
## D_kept/D ~ 1.0 removes the depth-dependent-filter confound noted below.
BASE <- list(N = 100, D = 300, lambda = 3000000, pi = 0.2)

## lambda scales with D as round(READS_PER_TAXON * D), holding reads/taxon fixed.
## Derived from BASE rather than hardcoded so the D = 300 entry is EXACTLY
## BASE$lambda -- main_grid() dedups on paste(N, D, lambda, pi), and if the D
## sweep's lambda at D = 300 missed 3e6 by a rounding unit the baseline would be
## simulated twice and appear as two points on the overLambda curve.
READS_PER_TAXON <- BASE$lambda / BASE$D

## Concentration of the effect-size distribution: |lfc| = 5*Beta(a,3a) + 1, whose
## mean is 2.25 for every a, so only the SPREAD moves. See make_lfcs().
BASE_EFFECT_CONC <- 1

## ---------------------------------------------------------------------------
## Sign split of the spiked effects -- the MAIN / SUPPLEMENTARY scenario switch.
## Added 2026-09-16.
## ---------------------------------------------------------------------------
## 0.8: 80% of spiked taxa increase, 20% decrease. This is the design of the main
##      manuscript and the design of every existing results tree.
## 0.2: the mirror, 20% up and 80% down. The supplementary scenario. NOTHING else
##      about the design changes -- same grid, same |lfc| distribution, same
##      number of spiked taxa, same replicate count, same seeds.
##
## A STUDY-WIDE constant, so it is handled exactly like OTHER_ROW: read from the
## environment, and carried in the TAG rather than in the param_string. Two
## scenarios must not share a results tree -- the directory name is the resume
## key, so without it in the tag a supplementary replicate would look "already
## done" against the main tree and silently never run.
##
## Set it on the COMMAND LINE, not by editing this file:
##     SIM_PROP_POS=0.2 Rscript run_main.R --sweep overN
## run_main.R launches each replicate through system() from a PSOCK worker, and
## both inherit the environment, so the value reaches the worker unchanged.
EFFECT_PROP_POS <- as.numeric(Sys.getenv("SIM_PROP_POS", "0.8"))
stopifnot(is.finite(EFFECT_PROP_POS), EFFECT_PROP_POS > 0, EFFECT_PROP_POS < 1)

## The suffix a non-default scenario puts on its tag AND on every output name.
## "" for the main design, so nothing existing moves or is renamed.
SCENARIO_SFX <- if (isTRUE(all.equal(EFFECT_PROP_POS, 0.8))) "" else
  sprintf("_pos%02d", round(100 * EFFECT_PROP_POS))

## The inverse: the scenario suffix a TAG carries. Defined here so the aggregator
## and the plotter derive the same output names from the same tag, whether or not
## SIM_PROP_POS happens to be set in the shell they are run from -- aggregating a
## supplementary tree must never be able to overwrite the main figures just
## because the env var was forgotten.
scenario_sfx_of <- function(tag) {
  m <- regmatches(tag, regexpr("_pos[0-9]{2}$", tag))
  if (length(m)) m else ""
}

## ---------------------------------------------------------------------------
## Phase 0 grid
## ---------------------------------------------------------------------------
## Two D sweeps on purpose. The existing scale_choice study ties lambda = 150*D,
## so its D axis also moves depth and the two effects cannot be separated --
## scale_model_report.Rmd finds D up -> bootstrap and lambda up -> laplace, i.e.
## they push opposite ways. overD_fixedLambda isolates D; overD_fixedReads keeps
## reads/taxon constant and asks the "wider study, same effort per taxon"
## question. `corner` probes the large-D / low-lambda cell where both effects
## point at bootstrap.
PILOT_N_REP <- 50

pilot_grid <- function() {
  cell <- function(sweep, N, D, lambda, pi = BASE$pi, sigma_scale = 1, tau = 0,
                   effect_conc = BASE_EFFECT_CONC)
    data.frame(sweep = sweep, N = N, D = D, lambda = lambda, pi = pi,
               sigma_scale = sigma_scale, tau = tau, effect_conc = effect_conc,
               stringsAsFactors = FALSE)

  Dv <- c(20, 50, 100, 200, 400, 800)

  rbind(
    cell("overN",             c(10, 16, 20, 40, 60, 100, 200), BASE$D, BASE$lambda),
    ## PRIMARY D sweep: lambda rises with D, holding reads/taxon at 200. This is
    ## the "a wider study sequences proportionally deeper" question, and it keeps
    ## D_kept/D roughly constant (~0.45 measured) instead of letting the
    ## depth-dependent filter collapse D_kept at large D.
    cell("overD_fixedReads",  BASE$N, Dv, READS_PER_TAXON * Dv),
    ## Retained as the control: it is the only way to see a D effect that is not
    ## a depth effect, since the two push opposite ways on the scale rule.
    cell("overD_fixedLambda", BASE$N, Dv, BASE$lambda),
    ## Recentred for D = 200. The previous grid started at lambda = 1e3, which
    ## here is 5 reads/taxon and collapses D_kept to ~2 -- a degenerate cell, not
    ## a low-depth one. Aligned to main_grid()'s list on 2026-09-02 so the pilot
    ## sweep CONTAINS its own baseline (3e6); the previous list topped out at
    ## 2.5e6 and would not have. 33-33,333 reads/taxon; real is 24-689.
    cell("overLambda",        BASE$N, BASE$D,
         c(1e4, 1e5, 5e5, 1e6, 3e6, 1e7)),
    cell("overPi",            BASE$N, BASE$D, BASE$lambda,
         pi = c(0.02, 0.05, 0.10, 0.20, 0.30)),
    ## depth-starved corner: wide studies at low reads/taxon (25-250 here)
    do.call(rbind, lapply(c(2e4, 1e5), function(lam)
      cell("corner", BASE$N, c(400, 800), lam))),

    ## ---- heterogeneity probe -------------------------------------------
    ## The four axes above returned P(bootstrap) = 0 in all 37 cells
    ## (0 / 3700, 95% UB 0.08%), with margin ~ -1: the rule is decisively
    ## Laplace everywhere on N, D, lambda and pi. Heterogeneity is the only
    ## lever the earlier no-filter study ever saw move the choice, so it is
    ## the remaining place a bootstrap regime could live.
    ##
    ## INDEPENDENT (per-taxon) heterogeneity: scales the template's per-feature
    ## sigma. Pushed BELOW the earlier study's floor of 0.15, because
    ## scale_model_report.Rmd places bootstrap selection at the LOW end of the
    ## real sd_between_sample range -- Jin2022, the one real dataset that picks
    ## bootstrap, sits at 0.247, with McMurrough 0.285 and Stammler 0.658.
    cell("overSigma", BASE$N, BASE$D, BASE$lambda,
         sigma_scale = c(0.05, 0.10, 0.15, 0.30, 0.50, 0.75, 1.00, 1.50)),

    ## SHARED (compositional) heterogeneity, crossed with pi. These two move
    ## OPPOSITE sides of the comparison, which is why every one-dimensional
    ## sweep so far came back flat:
    ##     tau -> v_B   (mode cancellation)
    ##     pi  -> v_L   (signal density)
    do.call(rbind, lapply(c(0.02, BASE$pi), function(p)
      cell("overShared", BASE$N, BASE$D, BASE$lambda,
           pi = p, tau = c(0, 0.05, 0.10, 0.15, 0.25, 0.40)))),

    ## ---- effect-size concentration -------------------------------------
    ## |lfc| = 5*Beta(a,3a) + 1 has mean 2.25 at every a, so this varies the
    ## SPREAD of the effect sizes with their magnitude held fixed:
    ##   conc      1     3     10    30    Inf
    ##   sd|lfc|   0.97  0.60  0.34  0.20  0
    ## The question is whether a sharply peaked effect distribution changes which
    ## scale model wins. Mechanism: nrd0 = 0.9*min(sd, IQR/1.34)*D^(-1/5) and
    ## v_L ~ h^2, so a wide effect distribution inflates `sd` -- but only reaches
    ## the bandwidth if the sd branch of the min() binds. At pi = 0.1 the nulls
    ## occupy both quartiles, so IQR/1.34 is the smaller and the effect spread is
    ## discarded. Crossed with pi because that is what decides which branch binds.
    do.call(rbind, lapply(c(BASE$pi, 0.30), function(p)
      cell("overEffect", BASE$N, BASE$D, BASE$lambda, pi = p,
           effect_conc = c(1, 3, 10, 30, Inf))))
  )
}

## Seeds are namespaced BY SWEEP so running one sweep at a time cannot reuse
## another's seeds. The offset clears the existing scale_choice runs (1e6 block)
## so old and new results never share a replicate.
PILOT_SEED_OFF <- 2000000L
PILOT_SEED_BLOCK <- c(overN             = 0L,
                      overD_fixedLambda = 100000L,
                      overD_fixedReads  = 200000L,
                      overLambda        = 300000L,
                      overPi            = 400000L,
                      corner            = 500000L,
                      overSigma         = 600000L,
                      overShared        = 700000L,
                      overEffect        = 800000L)

## ---------------------------------------------------------------------------
## Study A -- set at the Phase 0 gate.
## ---------------------------------------------------------------------------
## Phase 0 found the scale-model choice flat across every axis (1 bootstrap in
## 2850; margin monotone in N and D but never crossing 0), so the grid is set for
## the METHOD COMPARISON and for real-data resemblance, not for scale choice.
##
## Each sweep spans the real-data range and extends past it into the asymptotic
## regime. lambda is tied to D in the D sweep, holding reads/taxon at 200.
## 50 as of 2026-09-02, down from 100. The deeper baseline retains ~299 of 300
## features against the previous 237, so every cell costs more per replicate and
## the replicate count is where that is bought back.
##
## The cost is precision: Monte Carlo SE scales as 1/sqrt(n_rep), so halving the
## replicates widens it by sqrt(2). On the previous pi sweep the Sparse SSRV TDR
## se at pi = 0.65 was 0.0165 at n = 100; expect about 0.023 at n = 50.
##
## This is RECOVERABLE without re-running: run_main.R resumes by replicate index,
## so `--reps 100` later tops the SAME tree up from 50 to 100 and the two halves
## are one homogeneous set -- PROVIDED no code changed in between. If anything in
## config.R or src/ changes first, that top-up is a code mix and needs a new tag.
MAIN_N_REP <- 50
RESULT_TAG <- format(Sys.Date(), "%m%d%Y")

## Study A runs one sweep at a time over several days, so its tag is FROZEN at the
## date the run started rather than tracking Sys.Date(). The tag is the resume key
## -- run_main.R decides a replicate is done by looking for its directory -- so a
## tag that rolls over at midnight re-simulates the baseline into a second tree
## and splits the result set between them. Pinned by hand, and it must be repinned
## by hand if BASE ever changes: unlike the ssrv/sigma/nofilter tags, a bare date
## does not encode the baseline, so a new design would resume into the old tree.
##
## REPINNED 2026-09-02 for the depth change (BASE$lambda 5e5 -> 3e6, and with it
## READS_PER_TAXON 1667 -> 10000, so the overD ladder moves too). A NEW tag is
## mandatory: every cell in the grid changes, so 09012026 and 09022026 share no
## comparable replicate except by coincidence at lambda = 5e5, which now appears
## only as an overLambda POINT rather than as the baseline.
##
## The 2026-09-01 repin note, kept because the same argument applies: it must be
## a NEW tag, not --reps on the previous one, for two independent reasons:
##
##   1. The code changed. FILTER_TYPE is back to "mean" (D_kept ~300 -> 234 at
##      the baseline), which changes the feature set every method is scored on.
##      Nothing produced under 08292026 is comparable to what this tag produces.
##      08292026 itself was pinned for the same class of reason -- it carried the
##      switch TO "prevalence", the shared-theta^perp draw in sim_sparse_auto(),
##      and lib_cut / lib.cut going from 1000 to 0.
##
##   2. 08172026 was ALREADY a mix of code states, which is a warning about how
##      easy this is to miss. Its D = 300 baseline cell ran before the ALDEx2
##      mc.samples 2000 -> 1000 edit and the rest of the sweep ran after, so
##      ALDEx2_def is 4.7x slower at D = 300 than at D = 500 -- a NON-MONOTONE
##      timing curve, which is how it was found. The same window covers the
##      holm -> BH and edgeR/limma raw-p -> FDR fixes.
##
## Nothing in a tag directory records which code produced it, so the tag is the
## only thing standing between two code states and one results tree.
## ---------------------------------------------------------------------------
## Sweep axis labels -- ONE definition, shared by every figure.
## ---------------------------------------------------------------------------
## The lever a sweep varies has a code name (overD) and a symbol (D), and neither
## belongs on a published axis: a reader should not need to know that D is the
## feature count or pi the proportion of true associations. These are the plain
## English names, defined here because plot_studyA_combined.R and
## aggregate_main.R both draw sweep figures and had drifted to different labels
## -- the first to "Number of taxa (D)", the second to the raw "overD".
##
## The ORDER here is also the panel order of the combined Study A figure.
SWEEP_LABEL <- c(overD      = "Number of Features",
                 overN      = "Sample Size",
                 overLambda = "Sequencing Depth",
                 overPi     = "Proportion of True Associated Features")

## The label for a sweep, falling back to the sweep's own name so an unlabelled
## new sweep still plots -- with a visibly raw axis rather than a blank one.
sweep_label <- function(sw) {
  sw <- as.character(sw)
  ifelse(sw %in% names(SWEEP_LABEL), unname(SWEEP_LABEL[sw]), sw)
}

## The scenario rides in the tag (see EFFECT_PROP_POS): "09022026" is the main
## 80/20 design and "09022026_pos20" the supplementary 20/80 mirror, same code
## state and same grid. Written as one base string plus the suffix so the two can
## never drift to different baselines.
MAIN_TAG <- paste0("09022026", SCENARIO_SFX)

## MIMIX is a Julia MCMC fit per replicate (Mamba.jl, 2000 iters, 1 chain).
## Measured under a minute at most settings, so it is on by default.
RUN_MIMIX  <- TRUE
JULIA_BIN  <- Sys.which("julia")
MIMIX_SRC  <- file.path(OLD_SIM, "src", "mimix")

## ONE sweep definition, used by both Phase 0c (SSRV only) and Study A (13
## methods). Keeping them identical is the point: Phase 0c is a dry run of the
## Study A pipeline, so any difference between their numbers would have to come
## from the method set, never from the design.
##
## Both sweeps stop short of the original study's top end: lambda at 5e5 rather
## than 1e6+, pi at 0.65 rather than 0.75. pi = 0.75 means three quarters of taxa
## are differentially abundant, which is far outside anything real and is a
## regime where the identification argument is already strained.
main_grid <- function() {
  cell <- function(sweep, N, D, lambda, pi)
    data.frame(sweep = sweep, N = N, D = D, lambda = lambda, pi = pi,
               stringsAsFactors = FALSE)

  # 300 is the baseline and sits between the original's 200 and 350, so both are
  # dropped rather than crowding the curve with near-duplicates.
  Dv <- c(20, 50, 100, 300, 500, 700, 1000)

  g <- rbind(
    # the original study's N sweep exactly
    cell("overN",      c(40, 70, 100, 200, 300, 400), BASE$D, BASE$lambda, BASE$pi),
    # lambda tied to D, holding reads/taxon at BASE$lambda/BASE$D = 1666.67
    cell("overD",      BASE$N, Dv, round(READS_PER_TAXON * Dv), BASE$pi),
    # 17 -> 3333 reads/taxon at D = 300; the real studies run 68-692
    ## 33 / 333 / 1,667 / 3,333 / 10,000 / 33,333 reads/taxon at D = 300.
    ## The first four bracket every real study (24-689 r/t); 5e5 is retained
    ## because it was the 09012026 baseline, so that tree stays comparable as a
    ## point INSIDE this grid; 3e6 is the new baseline and the template's own
    ## depth; 1e7 is the asymptote and still inside the template's 97.5th
    ## percentile (1.44e7).
    cell("overLambda", BASE$N, BASE$D, c(1e4, 1e5, 5e5, 1e6, 3e6, 1e7), BASE$pi),
    # Extended to 0.75 on 2026-09-05 at the user's request. The original stop at
    # 0.65 was deliberate and its reasoning still holds: pi = 0.75 means three
    # quarters of taxa differentially abundant, which is outside anything real
    # and past the point where the scale model is identifiable -- read 0.75 as a
    # stress point, not as a realistic setting. Appended, not inserted: the seed
    # is the replicate index (worker_main.R:54), so adding a level cannot
    # perturb any existing cell, and main_grid() dedups on parameters, so the
    # 09022026 tree is topped up rather than re-simulated.
    cell("overPi",     BASE$N, BASE$D, BASE$lambda,
         c(0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.65, 0.75))
  )

  ## The baseline appears in all four sweeps. Deduplicate on the PARAMETERS so it
  ## is simulated once; `sweep` is only a plotting label and the aggregator
  ## re-attaches each setting to every sweep it belongs to.
  key <- paste(g$N, g$D, g$lambda, g$pi)
  g[!duplicated(key), ]
}

## Which sweeps a given parameter set belongs to -- the inverse of the dedup
## above, used by the aggregator to rebuild each sweep's curve.
main_sweeps_of <- function(N, D, lambda, pi) {
  s <- character(0)
  if (D == BASE$D && lambda == BASE$lambda && pi == BASE$pi) s <- c(s, "overN")
  if (N == BASE$N && lambda == round(READS_PER_TAXON * D) && pi == BASE$pi) s <- c(s, "overD")
  if (N == BASE$N && D == BASE$D && pi == BASE$pi) s <- c(s, "overLambda")
  if (N == BASE$N && D == BASE$D && lambda == BASE$lambda) s <- c(s, "overPi")
  s
}

## Phase 0c uses the same grid; only the method set and the replicate count
## differ. Defined as an alias rather than a copy so they cannot drift.
ssrv_grid  <- main_grid
SSRV_N_REP <- 50

## ---------------------------------------------------------------------------
## sigma_scale sweep -- does the dispersion of the NULL taxa drive the Laplace
## arm's FDR?
## ---------------------------------------------------------------------------
## Measured chain (350 baseline replicates): larger Laplace variance genuinely
## lowers FDR (v_L vs FDR_L: r = -0.167, p = 0.008), but effect-size MAGNITUDE
## does not move v_L at all (r = +0.024, p = 0.71). That is what the nrd0
## structure predicts -- h = 0.9*min(sd, IQR/1.34)*D^(-1/5), and at low pi the
## nulls occupy both quartiles, so IQR/1.34 binds and the spiked effects sit in
## tails that never touch it. So the lever on v_L is the NULL cluster's width,
## which sigma_scale sets.
##
## Swept BOTH directions. sigma_scale < 1 is the direction that matters: measured
## sd_between_sample runs 0.195 (at 0.05) to 1.44 (at 1.5), while the real studies
## sit at 0.247-0.658 -- i.e. BELOW the sigma_scale = 1 baseline of ~1.3. If the
## simulation's nulls were too tight we would need to go up; the measurements say
## the simulation is already MORE dispersed than reality on this axis, so the
## sweep has to cover both or it cannot distinguish the two stories.
##
## Crossed with pi because the FDR problem only exists at low pi: at the baseline
## pi = 0.2 the Laplace FDR is already 0.025, with nothing to move.
SIGMA_N_REP <- 50

## ---------------------------------------------------------------------------
## No-prevalence-filter check -- ONE cell.
## ---------------------------------------------------------------------------
## The filtered counterpart is already measured (ssrv_N50_D300_*, pi = 0.05:
## FDR_L 0.081, FDR_B 0.015, D_kept ~155), so only the unfiltered arm is run and
## the comparison is against real numbers rather than a re-run.
##
## min_mean = 0 keeps every taxon, so D_used = D = 300 and the rare taxa -- whose
## theta^|| is hugely dispersed, rdirichlet(alpha ~ 0.5) on near-zero counts --
## stay in the LFC density. That should widen nrd0 and hence v_L ~ h^2, the one
## quantity measured to lower FDR (v_L vs FDR_L: r = -0.167, p = 0.008).
NOFILTER_N_REP <- 100

## filter_cut = 0 keeps every taxon under EITHER rule (prevalence >= 0 and
## rowMeans >= 0 are both vacuous), so this cell stays the no-filter arm without
## needing to know which FILTER_TYPE is configured.
##
## NOTE: this check earned its keep against min_mean = 5 at lambda = 5e5, where
## it moved D_kept from ~155 to 300. It is VACUOUS AGAIN as of 2026-09-02, for a
## new reason: at BASE$lambda = 3e6 the mean filter keeps ~299 of 300 features
## anyway, so the filtered and unfiltered arms coincide at the baseline. To make
## this contrast mean anything now it has to be run at a SHALLOW lambda -- 1e4
## keeps only ~30% -- rather than at BASE.
## ---------------------------------------------------------------------------
## pi = 0.65 depth probe -- added 2026-09-02.
## ---------------------------------------------------------------------------
## Question it answers: Sparse SSRV's TDR at pi = 0.65 is ~0.14, and every
## simulation_rerun tag (08172026, 08292026, 09012026) shows the same thing, so
## it is NOT caused by the filter revert or by the rnorm(1) scale draw -- those
## two differ across those tags and the number does not move. All three ran at
## lambda = 5e5, so depth is UNTESTED, and lambda 7.5e5 -> 5e5 is the one BASE
## parameter changed from the original published design.
##
## This crosses that against OTHER_ROW, the second untested difference from the
## published pipeline (Phase 0 checked it against the scale-model CHOICE, never
## against power). OTHER_ROW is a study-wide constant and rides in the TAG, not
## the param_string, so each setting gets its own tree -- run the probe twice,
## with SIM_OTHER_ROW=fixed and SIM_OTHER_ROW=legacy.
## 750000 is not a guess: the AOAS manuscript's Figure caption states "Baseline
## values are 300 features, 100 samples, 750,000 reads, and 20% true
## associations", and simulation/src/06072025/Results over Pi.R sets
## `lambda <- 2500 * d`, i.e. 750000 at D = 300. Current BASE$lambda is 5e5,
## which is 1667 reads/feature against the published 2500.
PROBE_LAMBDA <- c(500000, 750000)   # current baseline, and the PUBLISHED one
## Which pi the probe sits at. Env-overridable so the depth contrast can be run
## at more than one sparsity without editing this file; it rides in the TAG (see
## run_main.R) so two pi values cannot share a results tree.
PROBE_PI <- as.numeric(Sys.getenv("SIM_PROBE_PI", "0.65"))
stopifnot(is.finite(PROBE_PI), PROBE_PI > 0, PROBE_PI < 1)
probe_grid <- function()
  data.frame(sweep = "probe", N = BASE$N, D = BASE$D,
             lambda = PROBE_LAMBDA, pi = PROBE_PI, sigma_scale = 1,
             stringsAsFactors = FALSE)

nofilter_grid <- function()
  data.frame(sweep = "nofilter", N = BASE$N, D = BASE$D, lambda = BASE$lambda,
             pi = 0.05, sigma_scale = 1, filter_cut = 0, stringsAsFactors = FALSE)

sigma_grid <- function() {
  sig <- c(0.05, 0.10, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00)
  do.call(rbind, lapply(c(0.02, 0.20), function(p)
    data.frame(sweep = "overSigma", N = BASE$N, D = BASE$D, lambda = BASE$lambda,
               pi = p, sigma_scale = sig, stringsAsFactors = FALSE)))
}

## Memory, not CPU, is the binding constraint, and it binds EVERYWHERE -- not just
## on the big cells, which is what the first version of this function assumed.
##
## MEASURED 2026-08-17, live RSS of a baseline worker (N=100, D=300, D_kept~235):
## 0.8 GB early in the method block, rising to 2.5 GB at peak. The resampling
## array is only 0.38 GB of that; the rest is the R process itself with phyloseq,
## DESeq2, ALDEx2, ANCOMBC, baySeq and MRCKF loaded, plus each method's working
## set. Sizing the pool from the array alone (the old comment's 0.21 GB) was off
## by 12x and put 7 workers x ~2 GB = 15 GB on a 16 GB machine: the run went into
## swap and per-worker CPU collapsed to 12-30%, i.e. ~4x slower than the same work
## on a smaller pool. Throughput is maximised by staying OUT of swap, not by
## filling the cores -- 4 workers not swapping beat 5 that are by ~2.5x, measured.
##
## Take the 2.5 GB as a LOWER bound on peak: it grows as a worker advances through
## the method block, so any mid-run sample understates it. Size with margin.
##
## Budget ~10-11 GB usable: 16 GB less ~2 GB wired and ~3.5 GB of desktop apps
## (measured with the browser CLOSED -- Chrome alone was 1.6 GB and cost a worker).
##   baseline-ish   ~2.5 GB/worker -> 4   (10 GB, ~1 GB slack)
##   N >= 200       array grows with N: 1.5 GB at N=400, so ~3.65 GB/worker -> 3
##   D >= 700       ~3.5 GB/worker at D=1000 (D_kept ~685)                  -> 3
##
## VERIFIED 2026-08-30 on the overN sweep, N = 200/300/400 at 3 workers: peak
## 8.39 GB total RSS with all three workers in the method block, and mean worker
## CPU 100%. Comfortably inside the budget, so 3 STAYS. Do not raise it -- 8.39
## is the peak observed, not the peak possible, and D >= 700 has not been tested.
##
## How to read this if you re-check it, because the obvious measurement lies:
##
##   * `sysctl vm.swapusage` is NOT the signal. macOS reported 12.8 GB "used" of
##     14.3 GB while every worker sat at 99% CPU, and the swap file SHRANK to
##     9.2 GB over the next two minutes. It is compressed-memory bookkeeping.
##     `vm_stat`'s "Pages occupied by compressor" was 264 MB -- the real figure.
##
##   * Count JULIA, not just R. During the MiMIX fit the R worker blocks in
##     system(jcmd) at ~0% CPU while a julia process does the work, so a filter
##     matching only R processes reports "no active workers" on a perfectly
##     healthy run. All three replicates can be in that phase at once.
##
## The signal that actually means swapping is the one measured on 2026-08-17:
## per-worker CPU collapsing to 12-30%. Match `exec/R|julia` and look at CPU.
##
## (2 rather than 3 would cost ~9 h on the N sweep, which is why this was worth
## measuring rather than assuming.)
## (The HPC version requested 20 GB per core for exactly this reason.)
workers_for <- function(N, D, n_cores = N_CORES) {
  if (N >= 200 || D >= 700) max(2L, min(n_cores, 3L)) else max(2L, min(n_cores, 4L))
}


## ---------------------------------------------------------------------------
## Study B -- oracle. Baseline plus the lambda sweep.
## ---------------------------------------------------------------------------
ORACLE_N_REP <- 50

## FDR convention. A method that flags nothing has FDR = 0/0. Reported as 0, per
## the original aggregators -- but ALWAYS print n_sig beside it, because
## "discovered nothing" and "discovered only true positives" both give 0 and are
## not the same result. The per-replicate metrics_summary.csv keeps the raw NaN,
## so the conditional-on-discovery version is recoverable without re-running.
FDR_WHEN_NO_DISCOVERIES <- 0
