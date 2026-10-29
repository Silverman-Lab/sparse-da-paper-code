################################################################################
## run_main.R -- Study A driver. Replaces foo.R + foo_R.sh + master_driver.sh.
##
## Iterates settings INSIDE R over a PSOCK pool instead of submitting one SLURM
## job per line of params.txt. Everything else is kept from foo.R:
##   * one Rscript SUBPROCESS per replicate, so a method that segfaults or
##     exhausts memory loses that replicate and not the worker;
##   * the replicate index is the trailing number of the directory name AND the
##     seed, so a replicate is reproducible standalone;
##   * the Julia MIMIX fit;
##   * the disk reduction -- beta.tsv (iters x taxa, tens of MB) collapsed to
##     beta_summary.csv, and the Y/X/Z handoff files deleted once Julia has run.
##
## Usage:
##   Rscript run_main.R                  # full grid, MAIN_N_REP reps, resumable
##   Rscript run_main.R --sweep overN    # ONE sweep to completion (overN, overD,
##                                       # overLambda, overPi); the baseline is
##                                       # shared, so whichever sweep runs first
##                                       # pays for it and the rest resume past it
##   Rscript run_main.R --benchmark      # 3 reps at baseline + 3 at the worst
##                                       # cell, with a projected total
##   Rscript run_main.R --reps 10        # smaller run
##   Rscript run_main.R --no-mimix       # skip the Julia fits -- NOTE this also
##                                       # deletes Y/X/Z, so MIMIX can never be
##                                       # added to those replicates afterwards
################################################################################

suppressWarnings(suppressMessages({
  library(foreach)
  library(doParallel)
}))

.bootstrap_config <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) {
    cfg <- file.path(p, "simulation_rerun", "config.R")
    if (file.exists(cfg)) return(cfg)
    p <- dirname(p)
  }
  stop("could not locate simulation_rerun/config.R above ", getwd())
}
source(.bootstrap_config())

args      <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) && length(args) > i[1]) args[i[1] + 1] else default
}
benchmark <- "--benchmark" %in% args
grid_sel  <- opt("--grid", "main")
sweep_sel <- opt("--sweep", "")              # "" => the whole grid
methods_s <- opt("--methods", "")            # "" => all 13, passed through to the worker

## Phase 0c is SSRV-only, so its tag and output tree are separate from Study A's
## and MIMIX is off by default -- MIMIX is a competitor method, not part of SSRV.
ssrv_only <- identical(methods_s, "sparse")
## The tag carries the BASELINE, not just the date: these runs are re-done at
## different baselines and a date-only tag would collide (or worse, silently
## resume into a run with a different design, since resumability keys on the
## directory name).
##
## The effect-size SIGN SPLIT rides here too, via config's SCENARIO_SFX -- the
## main 80/20 design adds nothing, the supplementary 20/80 mirror adds "_pos20".
## It is appended to EVERY grid's tag, not just the main one: resumability keys
## on the directory name, so a supplementary replicate that landed in a main tree
## would be read as "already done" and never run.
RESULT_TAG <- switch(grid_sel,
  ssrv  = sprintf("ssrv_N%d_D%d_%s",  BASE$N, BASE$D, RESULT_TAG),
  sigma = sprintf("sigma_N%d_D%d_%s", BASE$N, BASE$D, RESULT_TAG),
  nofilter = sprintf("nofilter_N%d_D%d_%s", BASE$N, BASE$D, RESULT_TAG),
  ## OTHER_ROW is in the tag because it is a study-wide constant, not a
  ## param_string field: two settings must not share a results tree or D_kept
  ## and the `other` row become unreadable.
  ## rule2 is a FRESH tag, not a top-up of probe_pi065_*: sparse_SSRV_rerun.R
  ## gained the width90mean recording between the two, and the project rule is
  ## that a results tree is one code state.
  probe    = sprintf("probe_pi%s_%s_v3_%s", sub("[.]","",format(PROBE_PI)), OTHER_ROW, RESULT_TAG),
  MAIN_TAG)
## MAIN_TAG already carries it (config.R); the others are built from a bare date.
if (grid_sel != "main") RESULT_TAG <- paste0(RESULT_TAG, SCENARIO_SFX)

use_mimix <- RUN_MIMIX && !("--no-mimix" %in% args) && !ssrv_only
## Per-grid replicate counts. A lookup rather than a chain of ifs, so adding a
## grid without adding its N_REP here fails loudly at the switch instead of
## silently inheriting MAIN_N_REP (which is what happened to the sigma sweep --
## it ran at 100 rather than SIGMA_N_REP = 50).
REPS_FOR <- list(main = MAIN_N_REP, ssrv = SSRV_N_REP, sigma = SIGMA_N_REP,
                 nofilter = NOFILTER_N_REP)
n_rep     <- if ("--reps" %in% args) as.integer(opt("--reps")) else {
  stopifnot(grid_sel %in% names(REPS_FOR)); REPS_FOR[[grid_sel]]
}

if (use_mimix && !nzchar(JULIA_BIN)) {
  warning("julia not found on PATH; continuing without MIMIX")
  use_mimix <- FALSE
}

RES_TAG_DIR <- file.path(RES_MAIN, RESULT_TAG)
dir.create(RES_TAG_DIR, recursive = TRUE, showWarnings = FALSE)
message(sprintf("grid=%s  methods=%s  reps=%d  mimix=%s  -> %s",
                grid_sel, if (nzchar(methods_s)) methods_s else "all",
                n_rep, use_mimix, RES_TAG_DIR))
## Printed unconditionally, including for the default: the scenario is invisible
## in the param_string and in every per-cell directory name, so the run log is
## where it has to be legible.
message(sprintf("effect-size sign split: %.0f%% positive / %.0f%% negative  (%s)",
                100 * EFFECT_PROP_POS, 100 * (1 - EFFECT_PROP_POS),
                if (nzchar(SCENARIO_SFX)) "SUPPLEMENTARY scenario" else "MAIN manuscript design"))

## ---------------------------------------------------------------------------
## Work list
## ---------------------------------------------------------------------------
grid <- switch(grid_sel,
               ssrv     = ssrv_grid(),
               sigma    = sigma_grid(),
               nofilter = nofilter_grid(),
               probe    = probe_grid(),
               main_grid())
## Grids may carry extra generative columns beyond N/D/lambda/pi. They ride along
## in the param_string as additional fields and in the directory name, so cells
## that differ only in them cannot collide or silently resume into each other.
EXTRA <- intersect(c("sigma_scale", "filter_cut"), names(grid))  # order is the field order

## ---------------------------------------------------------------------------
## --sweep: run ONE sweep to completion, so its curve can be read while the rest
## of the study is still unrun. Without it the whole grid is interleaved and no
## sweep finishes until nearly all of them do.
## ---------------------------------------------------------------------------
if (nzchar(sweep_sel)) {
  if (grid_sel != "main") stop("--sweep applies to the main grid only")
  ## Membership via main_sweeps_of(), NOT the grid's `sweep` column. main_grid()
  ## dedups on the PARAMETERS and keeps only the first label, so the baseline is
  ## labelled "overN" alone; filtering on the column would leave overD,
  ## overLambda and overPi each missing their D=300 / lambda=5e5 / pi=0.2 point
  ## -- a hole in the middle of every curve but the first.
  keep <- vapply(seq_len(nrow(grid)), function(i)
    sweep_sel %in% main_sweeps_of(grid$N[i], grid$D[i], grid$lambda[i], grid$pi[i]),
    logical(1))
  if (!any(keep)) stop("no settings belong to sweep ", sweep_sel)
  grid <- grid[keep, , drop = FALSE]
  message("sweep ", sweep_sel, ": ", nrow(grid), " settings")
}

## Run order: baseline first, then the cheap cells, then the expensive ones. A
## partial run then still yields complete curves for most sweeps rather than a
## uniformly thin slice of all of them. Sorted on D AND N: within a single sweep
## D is often constant (every overN cell sits at D = 300), and sorting on D alone
## would leave the N axis in grid order, running N = 400 before N = 40.
is_base <- with(grid, N == BASE$N & D == BASE$D & lambda == BASE$lambda & pi == BASE$pi)
rest <- grid[!is_base, , drop = FALSE]
grid <- rbind(grid[is_base, , drop = FALSE], rest[order(rest$D, rest$N), ])

if (benchmark) {
  worst <- grid[which.max(grid$D), , drop = FALSE]
  grid  <- rbind(grid[1, , drop = FALSE], worst)
  n_rep <- 3
  message("BENCHMARK: baseline + worst cell (D = ", worst$D, "), ", n_rep, " reps each")
}

safe <- function(...) gsub("[^a-zA-Z0-9]", "_", paste(...))

work <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  vals <- c(grid$N[i], grid$D[i], grid$lambda[i], grid$pi[i],
            if (length(EXTRA)) unlist(grid[i, EXTRA]) else NULL)
  data.frame(N = grid$N[i], D = grid$D[i], lambda = grid$lambda[i], pi = grid$pi[i],
             rep = seq_len(n_rep),
             param_string = paste(vals, collapse = " "),
             safe = do.call(safe, as.list(vals)),
             stringsAsFactors = FALSE)
}))
work$dir <- file.path(RES_TAG_DIR, sprintf("results_%s_%d", work$safe, work$rep))

## Resumability: a replicate is done when its metrics_summary.csv exists.
done <- file.exists(file.path(work$dir, "metrics_summary.csv"))
message(sprintf("%d settings x %d reps = %d replicates; %d already done, %d to run",
                nrow(grid), n_rep, nrow(work), sum(done), sum(!done)))
work <- work[!done, ]
if (nrow(work) == 0) { message("nothing to do"); quit(status = 0) }

## ---------------------------------------------------------------------------
## One replicate
## ---------------------------------------------------------------------------
run_replicate <- function(w, here, mimix_src, julia_bin, use_mimix, tag,
                          methods_s = "") {
  t0 <- proc.time()[["elapsed"]]
  dir.create(w$dir, recursive = TRUE, showWarnings = FALSE)

  worker <- file.path(here, "src", "worker_main.R")
  ## OMP/VECLIB pinned to 1: each worker is already one process in an outer pool,
  ## so letting a method's BLAS spawn its own threads oversubscribes the machine.
  cmd <- sprintf("cd %s && OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 %s %s %s %s %s",
                 shQuote(here), shQuote(file.path(R.home("bin"), "Rscript")),
                 shQuote(worker), shQuote(w$param_string), shQuote(w$dir),
                 shQuote(methods_s))
  status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ok <- file.exists(file.path(w$dir, "metrics_summary.csv"))
  t_methods <- proc.time()[["elapsed"]] - t0

  ## ----- MIMIX -------------------------------------------------------------
  t_mimix <- NA_real_
  if (ok && use_mimix) {
    mo <- file.path(mimix_src, "SparseDOSSA-analysis-results", tag, basename(w$dir))
    dir.create(mo, recursive = TRUE, showWarnings = FALSE)
    jcmd <- paste(shQuote(julia_bin), shQuote(file.path(mimix_src, "scripts/fit-mcmc.jl")),
                  "--hyper",   shQuote(file.path(mimix_src, "nutnet-analysis/configs/hyper.yml")),
                  "--monitor", shQuote(file.path(mimix_src, "nutnet-analysis/configs/monitor-mimix-no-factors.yml")),
                  "--inits",   shQuote(file.path(mimix_src, "nutnet-analysis/configs/inits.yml")),
                  "--factors 0 --iters 2000 --chains 1",
                  shQuote(w$dir), shQuote(mo))
    tj <- proc.time()[["elapsed"]]
    system(jcmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    t_mimix <- proc.time()[["elapsed"]] - tj

    ## fit-mcmc.jl writes the full <iters> x <taxa> sample matrix; the aggregators
    ## only ever take quantiles, so collapse it and drop the raw samples. This was
    ## the dominant disk cost in the original study (~24 GB across a sweep).
    bp <- file.path(mo, "beta.tsv")
    if (file.exists(bp)) {
      beta <- as.matrix(read.delim(bp, check.names = FALSE))
      bs <- t(apply(beta, 2, function(v)
        c(quantile(v, c(0.025, 0.5, 0.975), na.rm = TRUE), mean = mean(v, na.rm = TRUE))))
      write.csv(bs, file.path(mo, "beta_summary.csv"))
      unlink(bp)
    }
    unlink(file.path(mo, "beta_var.tsv"))
  }

  ## Y/X/Z existed only to hand data to Julia. Keep the taxa names, drop the rest.
  yp <- file.path(w$dir, "Y.csv")
  if (file.exists(yp)) {
    writeLines(scan(yp, what = character(), nlines = 1, sep = ",", quiet = TRUE),
               file.path(w$dir, "taxa.csv"))
  }
  unlink(c(yp, file.path(w$dir, "X.csv"), file.path(w$dir, "Z.csv")))

  data.frame(dir = basename(w$dir), N = w$N, D = w$D, lambda = w$lambda, pi = w$pi,
             rep = w$rep, ok = ok, status = status,
             secs_methods = t_methods, secs_mimix = t_mimix,
             secs_total = proc.time()[["elapsed"]] - t0,
             stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## Drive
## ---------------------------------------------------------------------------
## Run in POOL-SIZE BATCHES, not one pool for everything. The resampling array is
## S x D_kept x N doubles, so the N>=200 / D>=700 cells need ~0.7-0.85 GB per
## worker; at 7 workers that is 5-6 GB and the machine has 16 GB total. Those
## cells run on a smaller pool. Sizing the whole run for the worst cell would
## waste most of the wall clock, since most cells are small.
work$pool <- mapply(workers_for, work$N, work$D)

run_batch <- function(w, pool) {
  message(sprintf("  batch: %d replicates on %d workers (N up to %d, D up to %d)",
                  nrow(w), pool, max(w$N), max(w$D)))
  cl <- parallel::makeCluster(pool, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  ## "pass", not "remove": with "remove" a systematically broken driver silently
  ## returns NULL for every replicate and the run looks like it merely produced
  ## nothing. Passing the condition through means the message is visible.
  res <- foreach(i = seq_len(nrow(w)),
                 .export = c("run_replicate", "HERE", "MIMIX_SRC", "JULIA_BIN",
                             "use_mimix", "RESULT_TAG", "methods_s"),
                 .errorhandling = "pass") %dopar% {
    run_replicate(w[i, ], HERE, MIMIX_SRC, JULIA_BIN, use_mimix, RESULT_TAG, methods_s)
  }
  bad <- !vapply(res, is.data.frame, logical(1))
  if (any(bad)) {
    msg <- unique(vapply(res[bad], function(e)
      if (inherits(e, "condition")) conditionMessage(e) else "non-data.frame result",
      character(1)))
    warning(sum(bad), " of ", length(res), " replicates errored in the driver: ",
            paste(head(msg, 3), collapse = " | "), call. = FALSE)
  }
  if (all(bad)) return(NULL)
  do.call(rbind, res[!bad])
}

t_start <- proc.time()[["elapsed"]]
# largest pool first: the cheap cells finish quickly and give an early read on
# whether anything is systematically broken before the expensive ones start
log <- do.call(rbind, lapply(sort(unique(work$pool), decreasing = TRUE),
                             function(p) run_batch(work[work$pool == p, ], p)))
elapsed <- proc.time()[["elapsed"]] - t_start

log_path <- file.path(RES_TAG_DIR, sprintf("run_log_%s.csv",
                                           format(Sys.time(), "%Y%m%d_%H%M%S")))
write.csv(log, log_path, row.names = FALSE)

message(sprintf("\n%d/%d replicates ok in %.1f min (%.1f s/replicate wall, %d workers)",
                sum(log$ok), nrow(log), elapsed / 60,
                elapsed / max(1, nrow(log)), N_CORES))
message("log -> ", log_path)

if (any(!log$ok)) {
  message("FAILED replicates:")
  print(log[!log$ok, c("dir", "status")])
}

## ---------------------------------------------------------------------------
## Benchmark projection
## ---------------------------------------------------------------------------
if (benchmark) {
  ## na.action = na.pass: with MIMIX off, secs_mimix is NA for every replicate and
  ## the formula interface would drop every row ("no rows to aggregate").
  per <- aggregate(cbind(secs_methods, secs_mimix, secs_total) ~ D, data = log,
                   FUN = function(z) mean(z, na.rm = TRUE), na.action = na.pass)
  cat("\n=== per-replicate seconds ===\n"); print(per, row.names = FALSE, digits = 3)

  ## Cost scales with D, so project each setting from the two measured points by
  ## linear interpolation in D rather than assuming the baseline holds everywhere
  ## -- the large-D cells are what actually determine the total.
  full <- if (grid_sel == "ssrv") ssrv_grid() else main_grid()
  ## reps_full is the count the REAL run would use, not the 3 the benchmark used
  reps_full <- if (grid_sel == "ssrv") SSRV_N_REP else MAIN_N_REP
  if (nrow(per) >= 2) {
    est <- approx(per$D, per$secs_total, xout = full$D, rule = 2)$y
  } else {
    est <- rep(per$secs_total[1], nrow(full))
  }
  total <- sum(est) * reps_full
  cat(sprintf("\nprojected: %d settings x %d reps\n", nrow(full), reps_full))
  cat(sprintf("  serial CPU time : %.1f h\n", total / 3600))
  cat(sprintf("  wall on %d cores: %.1f h\n", N_CORES, total / 3600 / N_CORES))
  cat("  (D-interpolated from the two measured cells; MIMIX ",
      if (use_mimix) "included" else "EXCLUDED", ")\n", sep = "")
}
