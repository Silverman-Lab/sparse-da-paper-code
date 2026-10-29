################################################################################
## worker_main.R -- ONE replicate of Study A.
##
## Invoked as a subprocess, exactly like the original pipeline:
##   Rscript src/worker_main.R "<N> <D> <lambda> <pi>" <output_dir>
## The replicate index is the trailing number of <output_dir> and IS the seed,
## matching simulation/src/05052025/sparse_analysis_trycatch.R:40.
##
## Run as a subprocess rather than in-process so a method that segfaults or
## exhausts memory takes down one replicate, not the whole worker pool.
##
## Differences from the original worker, all deliberate:
##   * data generation is simulate_one(), shared with the pilot and Study B, so
##     the three cannot drift apart;
##   * the method block is run_all_methods() rather than an inline copy;
##   * metrics come from tidy_metrics.R -- ONE classification rule for every
##     method (see the README note on convention);
##   * the three SSRV arms are stashed as three methods (edit (c));
##   * theta_perp_true is recorded (free, from a_spiked) so the bias columns are
##     well defined; Study A does not plot them, goal 3 owns that;
##   * `n` is not clobbered by the Dirichlet loop.
################################################################################

.bootstrap_config <- function(start) {
  p <- normalizePath(start, mustWork = FALSE)
  for (i in 1:6) {
    cfg <- file.path(p, "simulation_rerun", "config.R")
    if (file.exists(cfg)) return(cfg)
    p <- dirname(p)
  }
  stop("could not locate simulation_rerun/config.R above ", start)
}
source(.bootstrap_config(getwd()))
source(file.path(HERE, "src", "shared.R"))
shared_setup("methods")

args <- commandArgs(trailingOnly = TRUE)
param_string <- args[1]
output_dir   <- args[2]
## Optional 3rd argument: comma-separated method set. Defaults to all 13.
## Phase 0c passes "sparse" to run the method under test alone -- which still
## yields all three scale-model arms, since one sim_sparse_auto() call writes
## the hybrid, Laplace-only and bootstrap-only tables (edits (a) and (c)).
ALL_METHODS <- c("sparse", "aldex", "deseq2", "ancombc", "ancom",
                 "linda", "edger", "limma", "bayseq", "ancombc2", "mrckf")
methods_sel <- if (length(args) >= 3 && nzchar(args[3])) {
  strsplit(args[3], ",")[[1]]
} else ALL_METHODS
unknown <- setdiff(methods_sel, ALL_METHODS)
if (length(unknown)) stop("unknown method(s): ", paste(unknown, collapse = ", "))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed.num <- as.integer(tail(strsplit(output_dir, "_")[[1]], 1))
pl <- as.numeric(strsplit(param_string, " ")[[1]])
n <- pl[1]; d <- pl[2]; lambda <- pl[3]; pi <- pl[4]
## Optional 5th field: sigma_scale, the per-taxon heterogeneity multiplier.
## Absent (4-field param_string) means 1, which is a no-op and consumes no extra
## RNG, so every existing result stays reproducible from its param_string.
sigma_scale <- if (length(pl) >= 5 && is.finite(pl[5])) pl[5] else 1
## Optional 6th field: the FILTER THRESHOLD, read according to config's
## FILTER_TYPE -- a prevalence fraction when that is "prevalence", a mean count
## when it is "mean". Absent means config's FILTER_CUT. 0 keeps every taxon under
## either rule and, via apply_filter(), appends no `other` row -- so
## D_used == D exactly.
##
## The field is the THRESHOLD, not the rule: FILTER_TYPE is a study-wide constant
## like OTHER_ROW, so one results tree is one rule. Putting the rule in the
## param_string would let two rules share a tree and make D_kept unreadable.
filter_cut <- if (length(pl) >= 6 && is.finite(pl[6])) pl[6] else FILTER_CUT

##### data ####################################################################
## prop_pos comes from config's EFFECT_PROP_POS, i.e. from the environment, NOT
## from the param_string: it is a study-wide constant like OTHER_ROW and it rides
## in the tag, so one results tree is one scenario. 0.8 is the main manuscript's
## design and a no-op relative to every existing tree.
gen <- simulate_one(n = n, d = d, lambda = lambda, pi = pi, seed = seed.num,
                    other_row = OTHER_ROW, filter_type = FILTER_TYPE,
                    min_prev = if (FILTER_TYPE == "prevalence") filter_cut else FILTER_MIN_PREV,
                    min_mean = if (FILTER_TYPE == "mean")       filter_cut else FILTER_MIN_MEAN,
                    sigma_scale = sigma_scale,
                    prop_pos = EFFECT_PROP_POS)
Y <- gen$Y
x <- gen$x

## MIMIX handoff. run_main.R deletes these after the Julia fit and keeps taxa.csv.
write.csv(t(Y), file.path(output_dir, "Y.csv"), row.names = FALSE)
write.csv(x,    file.path(output_dir, "X.csv"), row.names = FALSE)
write.csv(seq_along(x), file.path(output_dir, "Z.csv"), row.names = FALSE)
write.csv(gen$metadata_df,          file.path(output_dir, "metadata_df.csv"))
write.csv(gen$filtered_metadata_df, file.path(output_dir, "filtered_metadata_df.csv"))
write.table(gen$keep_taxa,          file.path(output_dir, "keep_taxa.csv"))

## Recorded per replicate so the depth-dependence of the filter, and the true
## scale, are both auditable from the results rather than recomputed later.
write.csv(data.frame(N = n, D = d, lambda = lambda, pi = pi, seed = seed.num,
                     sigma_scale = sigma_scale,
                     ## Both the rule and its threshold, so replicate_info.csv
                     ## alone says how the feature set was chosen. min_mean is
                     ## kept under its old name so the aggregator and any script
                     ## reading older result trees still find the column.
                     filter_type = FILTER_TYPE, filter_cut = filter_cut,
                     min_mean = if (FILTER_TYPE == "mean") filter_cut else NA_real_,
                     D_kept = gen$D_kept,
                     median_depth = median(colSums(Y)),
                     prop_zero = mean(Y == 0),
                     sd_between_sample = gen$sd_between_sample,
                     theta_perp_true = gen$theta_perp_true,
                     n_spiked_kept = nrow(gen$filtered_metadata_df),
                     ## The sign split, requested and realised. ceiling() rounds
                     ## both arms, so the two differ at small num.spiked; both are
                     ## recorded because the first says which scenario was asked
                     ## for and the second says what the replicate actually got.
                     prop_pos = EFFECT_PROP_POS,
                     prop_pos_realised = gen$prop_pos_realised,
                     n_pos = gen$n_pos, n_neg = gen$n_neg,
                     other_row = OTHER_ROW,
                     methods = paste(methods_sel, collapse = "+"),
                     # WHICH rdirichlet is in scope. dirmult is attached first,
                     # but LaplacesDemon (pulled in by one of the DA packages)
                     # attaches later and masks it, so the method path uses
                     # LaplacesDemon's while the SSRV-only pilot uses dirmult's.
                     # Same distribution, different draws -- recorded so a load
                     # order change cannot silently alter results unnoticed.
                     rdirichlet_pkg = environmentName(environment(rdirichlet))),
          file.path(output_dir, "replicate_info.csv"), row.names = FALSE)

##### methods #################################################################
## "mrckf" is NOT in run_all_methods()'s default `methods` vector, so it has to
## be named explicitly or CKF silently never runs.
## MRCKF stays LAST: knockoff construction and cv.glmnet both consume RNG, so
## moving it would shift every other method's draws.
res <- run_all_methods(
  Y = Y, x = x, output_dir = output_dir,
  filtered_metadata_df = gen$filtered_metadata_df,
  n = n, d = d, S = S_DRAWS, pval = PVAL,
  methods = methods_sel,
  scale_model = "hybrid")

##### metrics #################################################################
tidy <- tidy_replicate(res, rep_id = seed.num, arm = "sample",
                       params = list(n = n, d_kept = gen$D_kept,
                                     lambda = lambda, pi = pi),
                       theta_perp_true = gen$theta_perp_true)
summ <- summarise_replicate(tidy)

write.csv(tidy, file.path(output_dir, "metrics_long.csv"),    row.names = FALSE)
write.csv(summ, file.path(output_dir, "metrics_summary.csv"), row.names = FALSE)

cat(sprintf("replicate %d done: D_kept=%d, %d methods, %d failed\n",
            seed.num, gen$D_kept, length(res$results), length(res$failures)))
