# simulation_rerun

Laptop re-run of the sparse_DA simulation study against the **current** scale model
rule. Separate from `simulation/`, which is a mirror of the Roar Collab (SLURM)
workspace and is left untouched — this folder sources from it read-only and never
writes to it.

Five goals:

1. Move the parameter baseline into the real-data regime, and record which scale model
   the rule selects at every setting.
2. Report Laplace-only / bootstrap-only / hybrid as three curves on FDR, TPR and F0.5.
3. Oracle study: every method in the ideal setting (true absolute abundance known)
   versus the practical setting (relative counts only). All scale-adjustment
   reporting lives here.
4. **theta^perp coverage** (added 2026-08-29). The simulation is the only place the
   true load contrast is known exactly, so it is the only place the scale model's
   interval can actually be scored. Built with `ssrv_perp_summary()` — the SAME
   function the real-data analysis uses — so simulated and measured coverage are one
   quantity. Written per replicate to `theta_perp*.csv`, aggregated to
   `theta_perp_coverage.csv`. Reported with bias, RMSE and width, never alone.
5. **The nine Sparse SSRV variants** (added 2026-08-29), supplementary: 3 arms x 3
   bandwidths (`nrd0`, `nrd`, `SJ-ste`), the same grid as the real data's
   `Task 3 sparse SSRV variant tables.R`. The three bandwidths reuse ONE Dirichlet
   array, so the grid costs ~+2.8% per replicate. Aggregated to `ssrv_variants.csv`;
   deliberately NOT stashed as methods, so nothing supplementary can reach a main
   figure.

## Status

| stage | script | state |
|---|---|---|
| Phase 0 — scale-choice pilot | `run_pilot.R`, `aggregate_pilot.R` | **built, run** |
| Study A — main sweep | `run_main.R`, `aggregate_main.R` | **built, benchmarked, running** |
| Study B — oracle | `run_oracle.R`, `aggregate_oracle.R` | not yet written |

Study B's grid is deliberately **not** in `config.R` yet; a placeholder would be the
thing that silently gets run.

### Study A — how it is being run

Grid set from the Phase 0 result: the original published baseline (`N = 100, D = 300,
lambda = 5e5, pi = 0.2`), 23 unique settings over four sweeps. Results live under the
**pinned** tag `results/main/09022026/` — `MAIN_TAG` in `config.R`, frozen rather than
tracking `Sys.Date()` because the run spans days and the tag is the resume key.

**2026-09-16: the supplementary sign-split scenario.** The manuscript's main
simulation spikes **80% of the associated taxa UP and 20% DOWN**. The
supplementary section is the **mirror — 20% up, 80% down** — and nothing else
about the design changes: same 24-setting grid, same 50 replicates, same seeds,
same `|lfc|` distribution, same method set, same filter, same depth, same code.

It is run from `run_supplement.sh`, which is the four sweeps in sequence:

```bash
nohup ./run_supplement.sh > results/main/supplement_pos20.log 2>&1 &
```

The lever is `SIM_PROP_POS`, read by `config.R` into `EFFECT_PROP_POS` and passed
to `simulate_one(prop_pos = )`. Like `SIM_OTHER_ROW` it is a **study-wide
constant, so it rides in the TAG and not in the `param_string`**:

| scenario | `SIM_PROP_POS` | tag | outputs |
|---|---|---|---|
| main manuscript | unset (0.8) | `09022026` | `studyA_combined_fdr_tpr_f0.5.pdf`, … |
| supplement | `0.2` | `09022026_pos20` | `studyA_combined_fdr_tpr_f0.5_pos20.pdf`, … |

Both halves of that matter. Resumability keys on the directory name, so without
the scenario in the tag a supplementary replicate would be read as "already done"
against the main tree and silently never run. And `results/figures/` is shared by
every tag, so every table and figure carries a `_pos20` suffix, derived from the
TAG by `scenario_sfx_of()` rather than from the environment — aggregating the
supplementary tree cannot overwrite the main figures even if `SIM_PROP_POS` is
forgotten in the shell the aggregation is launched from.

**The default path is unchanged, and that is tested rather than asserted.**
`prop_pos` defaults to `4/5`, the draws stay in the same order (negatives first),
and `ceiling(num.spiked * (1 - 4/5))` equals the original `ceiling(num.spiked *
1/5)` at every setting in both grids — `1 - 4/5` and `1/5` differ by one ulp, and
the ceiling of an integer minus an ulp is that integer. `tests/test_edits_additive.R`
asserts all of it, and replicate 1 of `results_100_300_10000_0_2` re-run under the
patched code is **byte-identical** to the copy stored in the `09022026` tree,
`metrics_summary.csv`, `scale_model.csv`, `theta_perp.csv` and all three SSRV arm
tables included.

**Why this is not a relabelling of the main result.** The majority sign is what
the compositional methods take as their background. With 80% of the signal going
up, the reference — the geometric mean, the CLR centre — is dragged up with it, so
the nulls appear to go down and the induced compositional bias is negative; at
20/80 it is the other way round. Sparse SSRV's KDE mode tracks that null cluster,
so it sits on the opposite side of the true load shift. No sign flip of the
output recovers one scenario from the other, which is why it is simulated rather
than argued.

`run_pilot.R` is deliberately left on the default 4/5: the pilot tags do not
encode the scenario, so a pilot run under `SIM_PROP_POS=0.2` would mix two designs
into one tree.

**2026-09-02: the depth baseline moved, `lambda` 5e5 -> 3e6.** The simulation now
runs at the depth of the template it is generated from: `SparseDOSSA2:::Stool` has
`depth_fit = (mu_depth 14.8764, sigma_depth 0.8198)`, i.e. a median of 2,888,909
reads over 332 features (8,702 reads/taxon), fitted on 239 samples. At 3e6 over
D = 300 the simulation runs 10,000 reads/taxon; the previous 5e5 was 1,667 and the
AOAS published 7.5e5 was 2,500 — both a quarter to a fifth of the template's own
depth. `READS_PER_TAXON` is derived, so the D sweep's multiplier moves with it,
from `c = 1,667` to **`c = 10,000`** (D = 20..1000 now spans lambda 2e5..1e7).

Two consequences, both accepted:

* **The `rowMeans >= 5` filter stops binding at the baseline.** `D_kept` is 237 at
  5e5 and ~299 of 300 at 3e6, so `OTHER_ROW` and `FILTER_TYPE` are inert there and
  the `nofilter` check cell is vacuous unless run at a shallow lambda.
* **`D_kept/D` is now ~1.0 across the whole D sweep** instead of ~0.45, which
  removes the depth-dependent-filter confound the D sweep previously carried.

`09012026` is the last 5e5 tree. Do **not** `--reps` onto it: every cell differs.
lambda = 5e5 survives as an `overLambda` POINT in the new grid, so that tree stays
comparable there.

**2026-09-01: `FILTER_TYPE` is back to `"mean"`** (`rowMeans(counts) >= 5`), matching
the real-data revert — Prasse, Stammler, Zemb and selex all filter that way again;
Vandeputte and Jin keep their source papers' rules. `MAIN_TAG` was repinned to
`09012026` for it. Do **not** `--reps` onto `08292026`: that tree was produced under the
10% prevalence rule and `D_kept` at the baseline differs by ~65 taxa, so the two are not
one result set. `08172026` is older still and was itself a mix of code states.

Run one sweep at a time so each curve can be read while the rest is still running:

```bash
nohup Rscript run_main.R --sweep overN --reps 50 > results/main/overN.log 2>&1 &
Rscript aggregate_main.R
```

then `overD`, `overLambda`, `overPi`. Staged at 50 replicates; `MAIN_N_REP = 100` is the
final target and the top-up pass is the same command with `--reps 50` dropped, since the
rep index is the seed and reps 51-100 are simply added. `aggregate_main.R` skips any
sweep with fewer than two settings so far, so partial runs do not emit single-point
figures.

**Benchmark, 2026-08-17** (`results/main/benchmark.log`, 6/6 replicates in 46.1 min,
all 13 methods, MIMIX on): **433 s per replicate at the baseline, 2036 s at D = 1000.**
The script's own projection reads 51.3 h wall on 7 cores, but it interpolates in `D`
only and so charges every `overN` cell the N = 100 price; correcting for N puts the full
study nearer 3-3.5 days at 100 reps. At D = 1000 the cost is two methods — ANCOM 812 s
and MRCKF 619 s of the 1858 s method block (MRCKF is 3.5 s at D = 300, a 179× jump).

**MIMIX is on**, and it is not optional after the fact: `run_main.R` deletes `Y/X/Z`
after every replicate whether or not the Julia fit ran, so a replicate produced with
`--no-mimix` can never have MIMIX added without being re-run.

## Running

```bash
cd simulation_rerun
Rscript run_pilot.R                     # all sweeps, PILOT_N_REP reps, config's OTHER_ROW
Rscript run_pilot.R overLambda          # one sweep
Rscript run_pilot.R corner 10           # sweep + n_rep override
Rscript run_pilot.R all 50 fixed        # sweep + n_rep + other-row variant
Rscript aggregate_pilot.R               # summary table + figures
```

Everything locates the project root by walking up for `sparse_SSRV.R`, so scripts run
from any working directory. No `.libPaths()`, no `/storage/...`, no SLURM.

## Layout

```
config.R              single source of truth: paths, baseline, grids, N_REP, N_CORES
run_pilot.R           Phase 0 driver
aggregate_pilot.R     Phase 0 summary + figures
src/
  sparse_SSRV_rerun.R COPY of simulation/src/sparse_SSRV_03122025.R + edits (a),(b),(d)
  simulate_one.R      the SparseDOSSA2 generation block, defined ONCE
  shared.R            loads packages; sources the unchanged files from ../simulation
results/
  pilot/              pilot_raw_*.csv, seeds_*.csv, pilot_summary.csv, run_*.log
  figures/            pilot_*.pdf
```

`simulate_one.R` exists because the original pipeline had three near-copies of the
generation block (the worker, the oracle worker, and `scale_choice/run_scale_choice.R`)
which had already drifted apart. Phase 0, Study A and Study B all call this one
function, so they generate identically by construction rather than by inspection.

## Changes to the method file

`src/sparse_SSRV_rerun.R` is a copy of `simulation/src/sparse_SSRV_03122025.R` with
three edits. The original is untouched, so `simulation/` stays reproducible.

**(a) Emit the bootstrap-only table.** The original builds the bootstrap-arm table and
then has `write.csv(...)` sitting *after* `return(out)`, commented out — so
`sparse_B_filter*.csv` never appeared. It is now written. Purely additive: no draw is
taken and no RNG consumed, so a replicate is bit-identical to one run with the
original. This is what makes goal 2's three curves come from **one** fit rather than
three: the hybrid, Laplace-only and bootstrap-only tables are all products of the same
`sim_sparse_auto()` call.

**(b) `scale_model = c("hybrid", "laplace", "bootstrap")`.** Forces the arm instead of
selecting it. Both `thetaB` and `thetaL` are computed unconditionally either way, so
forcing costs nothing. `scale_stat` and the hybrid decision are still computed and
recorded under a forced value, so `scale_model.csv` always reports what the rule *would*
have chosen (`chosen_hybrid`) alongside what was used (`chosen`). A forced `laplace`
run on a collapsed Laplace branch falls back to bootstrap and sets `forced_ok = FALSE`
rather than mislabelling itself.

**(d) Guard the subject bootstrap against losing an arm.** See below.

## Two pre-existing bugs found while porting

### 1. The subject bootstrap can lose an entire arm — and it is fatal at small N

`sample(1:N, replace = TRUE)` can return indices from one group only. Then
`log_wn_bootstrap[, treat.ind2]` is a `D x 0` matrix, `rowMeans()` of it returns `NaN`
rather than erroring, `theta_pll[, s]` becomes all-`NaN`, and the failure surfaces two
lines later as

```
density.default: need at least 2 points to select a bandwidth automatically
```

which takes down the **entire replicate** — `sim_sparse_auto()` is called outside the
per-method `tryCatch`, so every other method in that replicate dies with it.

Probability a replicate hits this is `1 - (1 - 2*(1/2)^N)^S`:

| N | S = 2000 |
|---|---|
| 10 | **~98%** |
| 16 | ~6% |
| 20 | ~0.4% |
| 40+ | negligible |

The HPC sweep ran `N >= 40` and never saw it. This re-run deliberately moves to the
real-data regime, where four of the six studies have `N = 10-16` — so without a guard
the method simply cannot be run there.

`simulation/src/oracle/sparse_oracle_scale.R:92-99` already documents this as the guard
`sim_sparse_auto()` is missing and implements it. Edit (d) is that same guard, so the
two paths now agree: redraw until both arms are present. Verified — 21/21 replicates
succeed at `N = 10-200` where the first smoke test lost one at `N = 20`, and the guard
fires ~4 times per 2000 draws at `N = 10` (predicted 3.9) and zero at `N >= 16`.

The redraw consumes extra `sample()` calls, so the RNG stream shifts from that point —
but only in draws that would otherwise have crashed. At `N >= 40` no redraw ever fires
and results are unchanged. `n_resample_reject` is recorded per replicate so the
frequency is auditable rather than assumed.

### 2. The aggregated `other` row is all zeros

`simulation/src/05052025/sparse_analysis_trycatch.R:125`:

```r
Y     <- sim$simulated_data[rownames(sim$simulated_data) %in% keep_taxa, ]
other <- colSums(sim$simulated_data[rownames(Y) %in% other_taxa, ])
```

The second line indexes with `rownames(Y)` — already the **kept** taxa — so the logical
is all `FALSE`, and being shorter than `nrow(sim$simulated_data)` it is recycled. Zero
rows are selected and `other` is a vector of zeros.

Consequences in every published number: the filtered-out reads are **discarded** rather
than aggregated, and SSRV's KDE mode estimation sees one extra taxon whose `theta^||` is
pure Dirichlet noise. `apply_shared_filter()` in `oracle_arms.R:127` does it correctly,
so the main and oracle paths currently disagree.

`simulate_one()` takes `other_row = c("legacy", "fixed")`. `legacy` reproduces the bug
exactly (comparability with the HPC results); `fixed` aggregates as intended. Phase 0 is
run **both ways on the same seeds** so the choice is settled by measurement rather than
argument. Whichever is used for Studies A and B must be stated in the write-up.

## Phase 0 design

Scale-model choice only — none of the 13 DA methods run. Three things separate it from
`simulation/scale_choice/run_scale_choice.R`:

1. **The prevalence filter is applied.** `scale_choice` deliberately drops it to keep D
   orthogonal to lambda, so it measures `sparse_auto()` on a table no analysis ever
   sees. Here `sim_sparse_auto()` runs on filtered data — the Study A path, and the
   path every real study uses. `D_kept` is recorded so the filter's depth-dependence
   stays visible.

   **The rule has changed twice, and as of 2026-09-02 it barely bites.** On
   2026-08-29 it went from `rowMeans >= 5` to `rowSums(counts > 0)/n >= 0.10`; on
   **2026-09-01 it went back to `rowMeans >= 5`**,
   following the real-data analysis both times (`config.R`: `FILTER_TYPE`, now
   `"mean"`). It is not a small change — `D_kept` at the baseline is 234 under the mean
   rule against ~300 of 300 under the prevalence rule — so every figure carries the
   rule that produced it. Results tagged `08292026` are the prevalence rule;
   `09012026` and everything before 2026-08-29 are the mean rule. `filter_type =
   "prevalence"` regenerates the middle period.

2. **Both D sweeps are run.** `overD_fixedReads` is primary — `lambda = 200*D`, holding
   reads per taxon constant, i.e. "a wider study sequences proportionally deeper". It
   keeps `D_kept/D` roughly constant (~0.45 measured) instead of letting the
   depth-dependent filter collapse `D_kept` at large D. `overD_fixedLambda` is retained
   as the control: with lambda fixed it is the only way to see a D effect that is not a
   depth effect, and the two push *opposite* ways on the scale rule
   (`scale_model_report.Rmd`: D up → bootstrap, lambda up → Laplace), so a sweep that
   moves both cannot identify either.

3. **pi is swept.** It never has been for scale choice.

### Baseline

`N = 20, D = 50, lambda = 10000, pi = 0.1` — 200 reads/taxon, `D_kept ≈ 33`. The real
studies run `D = 27-237` post-filter at 68-692 reads/taxon, so this sits inside their
range on both axes.

`pi = 0.1` rather than the original study's 0.2: five of the six real datasets have
0.005-0.079 true associations per taxon (Vandeputte 0.79 is the outlier), so 0.2 sat
3-40x above anything real. 0.1 moves toward that without reaching the very sparse end,
where a 50-taxon simulation would spike only 2-3 features and every detection metric
would be dominated by its own denominator.

It also uses `dirmult::rdirichlet`, matching Study A's worker. `scale_choice` uses
MCMCpack's, which is a different RNG stream — one reason its numbers are not directly
comparable to the main sweep.

### A structural constraint the pilot exposes

NOTE: this section describes `FILTER_TYPE = "mean"` — the rule in force until
2026-08-29 and **again from 2026-09-01**, so the bound below is live for the `09012026`
run. It does NOT apply to the `08292026` (prevalence) tree: there a taxon survives on
being NON-ZERO in 10% of samples, which no depth budget constrains the same way, so
large D at low depth is reachable and the measured `D_kept` at `lambda = 5e3, D = 300`
is 184 rather than 74.

A taxon survives on `rowMeans >= 5`, and the row means sum to the mean depth, so

```
D_kept  <=  mean_depth / 5
```

**by construction**. Large D at low depth is therefore unreachable: at `D = 800,
lambda = 1000` only 36-48 taxa survive (`D_kept/D ~ 0.05`). Under the filter, `D_kept`
is governed by depth, not by D — which is exactly why the D and lambda axes cannot be
separated by a sweep that moves both. `pilot_dkept.pdf` plots `D_kept` against this
ceiling.

### Reading the output

`p_bootstrap` and `margin` must be read together. Under `scale_rule = "majority"`,
`statB + statL = 1`, so `margin = statB - statL` lives on `[-1, 1]`:

- `margin = -1` — Laplace contributes more variance for **every** taxon; a decisive call
- `margin = 0` — a dead tie; the two models are *indistinguishable* at that setting,
  which is a different finding from a real preference

A 50/50 split at `margin ~ 0` must not be read as a genuine preference. The bar chart
alone cannot separate the two.

`p_bootstrap_old` is the superseded mean-sd rule evaluated on the **same fits**. It is
the cleanest available evidence on whether the rule change, rather than the baseline, is
what moved the simulation into line with the real-data analysis.

## Phase 0b RESULT — the effect-size hypothesis, and why nothing moves the choice

Baseline `N = 40, D = 200, lambda = 500000, pi = 0.1`. 66 cells x 50 reps =
**3300 replicates, zero failures**, 49 min.

**Bootstrap chosen 0 / 3300.** `margin` median **−1.00** — for the median replicate not a
single taxon favours the bootstrap arm. The superseded rule chose bootstrap **89.2%** on
the same fits.

### `overEffect`: sharpening the effect-size distribution does nothing

`|lfc| = 5*Beta(a,3a) + 1` holds the mean at 2.25 for every `a`, so only the spread moves.
Swept at two sparsities:

| pi | `effect_conc` | realised sd(\|lfc\|) | P(bootstrap) | margin | med_vB | med_vL |
|---|---|---|---|---|---|---|
| 0.1 | 1 (original) | 0.946 | 0 | −0.994 | −0.161 | 0.577 |
| 0.1 | 3 | 0.606 | 0 | −0.995 | −0.148 | 0.608 |
| 0.1 | 10 | 0.342 | 0 | −0.994 | −0.190 | 0.595 |
| 0.1 | 30 | 0.196 | 0 | −0.993 | −0.178 | 0.578 |
| 0.1 | Inf (point mass) | 0.000 | 0 | −0.996 | −0.181 | 0.659 |
| 0.3 | 1 → Inf | 0.964 → 0.000 | 0 | −0.999 → −1.000 | ~−0.15 | ~1.0 |

Margin varies by **0.006** across the entire range — noise. Even a literal point mass at a
single effect size does not move it.

### The mechanism, and why it is structural

The two arms build `theta^perp` differently, and the difference is not about the data:

| | bootstrap | laplace |
|---|---|---|
| `theta^perp` | `-mode(KDE)` — deterministic scalar | `-mode` **+ `rnorm(D, 0, sqrt(variance))`** |
| independent noise added | **none** | the Laplace variance |
| `Cov(theta^||, theta^perp)` | strongly negative — the mode *is* the common component | same negative term |

`v(d) = Var(theta^perp) + 2Cov(theta^||_d, theta^perp)`. The Laplace arm gets an additive
independent `variance` term with **no covariance offset**; the bootstrap arm has no
counterpart for it and is left with the negative mode-cancellation covariance — the very
cancellation the method exists to perform.

Measured consequence, over 10850 replicates spanning every baseline tried:

| | |
|---|---|
| `med_vL > 0` | **100%** |
| `med_vB < 0` | 41.6% at `N=20, D=50, lam=1e4`; **93.1%** at `N=40, D=200, lam=5e5` |
| `med_vL > med_vB` | **99.99%** (median ratio 13x) |
| `statB` (fraction of taxa favouring bootstrap), median | **0** |
| replicates with `statB > 0.5` | **1 of 10850** |

`v_B` goes negative more often at the deeper, larger baseline — consistent with
mode-cancellation dominating once Dirichlet noise is small.

**So under `scale_rule = "majority"` the rule is close to predetermined: it compares an arm
that injects independent noise against one that does not.** In practice the hybrid is
"always Laplace, unless the Laplace branch collapses." That is a defensible conservative
choice — Section 4.3 motivates it for Type I control — but it should be reported as such
rather than as a data-adaptive selection. No data-generating lever tested (N, D, lambda,
pi, per-taxon heterogeneity, shared heterogeneity, effect-size spread) changes it.

A like-for-like comparison would either add the Laplace `variance` term to the bootstrap
arm as well, or compare the two on `Var(theta^perp)` alone without the covariance term.
That is a method decision, not a simulation one, and is left open.

## Phase 0 RESULT

Baseline `N = 20, D = 50, lambda = 10000, pi = 0.1`. 57 cells x 50 reps =
**2850 replicates, zero failures**, 17.2 min on 7 cores.

| | |
|---|---|
| bootstrap chosen, current rule | **1 / 2850** |
| bootstrap chosen, superseded mean-sd rule, same fits | **87.6%** |
| `margin` range across cells | −0.85 to −1.00 |

### The proportion is the wrong statistic; the margin is the right one

`P(bootstrap)` is 0 in 56 of 57 cells, so it cannot show a gradient. `margin =
statB − statL` can: under `scale_rule = "majority"` the two statistics sum to 1, so
margin lives on [−1, 1], and −1 means Laplace contributes more variance for *every*
taxon. It moves smoothly and monotonically:

| lever | margin across the range | reads as |
|---|---|---|
| **D**, fixed lambda (20 → 800) | −0.853 → −0.998 | bootstrap least competitive at large D |
| **D**, lambda = 200*D (20 → 800) | −0.849 → −0.997 | same, so it is D and not depth |
| **N** (10 → 200) | −0.942 → −0.993 | same direction |
| lambda (1e3 → 7.5e5) | −0.906 → −0.961 → −0.925 | shallow, non-monotone |
| pi (0.02 → 0.30) | −0.947 → −0.979 | mild |
| `sigma_scale` (0.05 → 1.5) | −0.95 to −0.99 | flat |
| `tau` (0 → 0.40) | −0.94 to −0.975 | flat |

The one bootstrap selection is at the smallest cell, `D = 20, lambda = 1e4`
(`D_kept ~ 18`): P = 0.02, Wilson CI [0.004, 0.105]. Exactly where the gradient says the
call is closest to a tie.

### The direction, and why it is the one to expect

**Bootstrap is most competitive at SMALL N and SMALL D, and becomes uniformly less
competitive as either dimension grows** — the opposite of "bootstrap appears as we enter
the asymptotic regime", and the two dimension levers are the only ones that move it at
all.

This follows from the method's own construction. Section 4.3 of the paper assigns the
bootstrap scale model the finite-N and finite-D components of scale uncertainty, and the
Laplace model the dispersion of `g^||` near its mode. The bootstrap's contributed
variance therefore shrinks as N and D grow, while the Laplace variance is a population
quantity that does not. Asymptotically Laplace must win; the sweeps measure how fast.
The `sigma_scale` and `tau` axes being flat says the same thing from the other side —
heterogeneity changes the *width* of `g^||` but not the *rate* at which the bootstrap
term vanishes.

### It agrees with the real data

`scale_model_choice/scale_model_comparison.csv` (2026-08-12) lists Jin2022 as bootstrap,
but it predates the closed-form curvature fix of 2026-08-14.
`sparse_ssrv_method_changes.Rmd` records that fix flipping Jin2022 and Prasse2015 from
bootstrap to laplace at 5/5 seeds each, every flip in that direction, because the old
`numDeriv::hessian()`-of-`approx()` path understated `vL` and biased the rule toward
bootstrap. Under the current method the real datasets are **6/6 laplace**.

So simulation and real data now agree, and the original premise — that the simulation
picked bootstrap while the real analysis picked laplace — was an artefact of the
superseded rule, not of the parameter baseline. **Do not use that stale CSV as the
real-data reference.**

### Consequence for the study

The hybrid rule is effectively Laplace throughout the practical regime, so
`Sparse_SSRV` == `Sparse_SSRV_L` in Study A on all but a handful of replicates. That
makes goal 2's three-curve decomposition the load-bearing result: the bootstrap-only
curve is what shows readers what the variance-inflation rule is protecting against,
since the hybrid almost never selects it. Report the proportion with its Wilson interval
and the margin beside it — a bare 0 with no interval overstates what 50 replicates can
establish, and hides the gradient entirely.

### First pass, at pi = 0.2

Archived under `results/pilot/archive_pi0.2/` and `results/figures/archive_pi0.2/`.
57 cells, **4700 replicates, zero failures**.

| | |
|---|---|
| bootstrap chosen, current rule | **0 / 4700** (95% upper bound 0.064%) |
| bootstrap chosen, superseded mean-sd rule, same fits | 3292 / 4700 = **70.0%** |
| `margin` min / median / max | −1.00 / **−1.00** / −0.059 |

The rule selects Laplace at **every** setting, and does so decisively: at the median
replicate, Laplace contributes more variance than bootstrap for *every taxon*. Nothing
moves it — not N (10–200), D (20–800), lambda (1e3–7.5e5), pi (0.02–0.30), the
large-D/low-lambda corner, per-taxon heterogeneity (`sigma_scale` 0.05–1.5), or shared
heterogeneity (`tau` 0–0.40).

The heterogeneity probe reached `sd_between_sample` **0.219 to 2.168**, which brackets
every real dataset — Jin2022 0.247, McMurrough2014 0.285, Stammler2016 0.658. So the
simulation does span the real regime; it is not that the sweep missed it.

**This agrees with the real data.** `scale_model_choice/scale_model_comparison.csv`
(2026-08-12) lists Jin2022 as bootstrap, but it predates the closed-form curvature fix of
2026-08-14. `sparse_ssrv_method_changes.Rmd` records that fix flipping Jin2022 and
Prasse2015 from bootstrap to laplace at 5/5 seeds each, every flip in that direction,
because the old `numDeriv::hessian()`-of-`approx()` path understated `vL` and biased the
rule toward bootstrap. Under the current method the real datasets are **6/6 laplace**.

So simulation (0/4700) and real data (0/6) now agree completely, and the original premise
— that the simulation picked bootstrap while the real analysis picked laplace — was an
artefact of the superseded rule, not of the parameter baseline. Do not use that stale CSV
as the real-data reference.

**Consequence for the study.** The hybrid rule is, empirically, Laplace throughout the
practical regime, so `Sparse_SSRV` == `Sparse_SSRV_L` in Study A. That makes goal 2's
three-curve decomposition the load-bearing result: the bootstrap-only curve is what shows
readers what the variance-inflation rule is protecting against, since the hybrid never
selects it. Report the proportion with its Wilson interval rather than as a bare 0 — 50
replicates per cell bound it at 7%, and it is the pooled 4700 that bounds it at 0.06%.
