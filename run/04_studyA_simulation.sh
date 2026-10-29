#!/bin/zsh
# Stage 4 -- Study A (simulation_rerun), both sign-split scenarios.
#
#   main manuscript  80% of spiked taxa up   tag 09022026         about 3 days
#   supplement       20% up / 80% down       tag 09022026_pos20   about 3 days
#
# 24 settings x 50 replicates each, all 13 methods, MiMIX on. The main scenario
# follows the commands in simulation_rerun/README.md (one sweep at a time); the
# supplement is simulation_rerun/run_supplement.sh unchanged apart from its cd.
# Both are resumable: the replicate index is the seed, and a replicate counts as
# done once its metrics_summary.csv exists, so re-running this script picks up
# where it stopped. See run_supplement.sh for the interrupted-MiMIX preflight.
#
# Outputs: simulation_rerun/results/main/<tag>/ (replicates and aggregated CSVs),
# simulation_rerun/results/figures/ (main_*, decomp_*, scale_choice_*,
# studyA_*, runtime_*; the supplement's carry a _pos20 suffix), and MiMIX fits in
# simulation/src/mimix/SparseDOSSA-analysis-results/<tag>/.

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export JULIA_PROJECT="$ROOT/simulation/src/mimix/julia_env"   # pinned MiMIX environment
cd "$ROOT/simulation_rerun" || exit 1
mkdir -p results/main

# ---- main scenario: SIM_PROP_POS unset (0.8) -> tag 09022026 -----------------
unset SIM_PROP_POS
TAG=09022026
for SW in overN overD overLambda overPi; do
  echo "== $TAG / $SW: $(date) =="
  Rscript run_main.R --sweep "$SW" --reps 50 > "results/main/${SW}_${TAG}.log" 2>&1 \
    || echo "RUN FAILED: $SW"
done
Rscript aggregate_main.R "$TAG"       > "results/main/aggregate_${TAG}.log" 2>&1 || echo "AGGREGATE FAILED"
Rscript plot_studyA_combined.R "$TAG" || echo "STUDYA PLOT FAILED"
Rscript plot_runtime.R "$TAG"         || echo "TIMING PLOT FAILED"

# ---- supplementary scenario: SIM_PROP_POS=0.2 -> tag 09022026_pos20 ----------
zsh ./run_supplement.sh > results/main/supplement_pos20.log 2>&1 || echo "SUPPLEMENT FAILED"
