#!/bin/zsh
################################################################################
## run_supplement.sh -- Study A under the SUPPLEMENTARY sign split.
##
## The manuscript's main simulation spikes 80% of the associated taxa UP and 20%
## DOWN. The supplementary section is the MIRROR: 20% up, 80% down. Nothing else
## differs -- same grid (24 settings), same 50 replicates, same seeds, same
## method set, same filter, same depth, same code.
##
## The scenario is a study-wide constant, so it rides in the TAG exactly like
## OTHER_ROW does:
##     SIM_PROP_POS=0.2  ->  config's EFFECT_PROP_POS  ->  SCENARIO_SFX "_pos20"
##     ->  results/main/09022026_pos20/
## That is what stops a supplementary replicate from being read as "already done"
## against the main 09022026 tree, since resumability keys on the directory name.
##
## Runs the four sweeps ONE AT A TIME, in the README's order. Not concurrently:
## each sweep sizes its own worker pool for its largest cell (config's
## workers_for()), and two pools at once puts the machine into swap, which the
## 2026-08-17 measurement shows costs more than it buys.
##
## Resumable. A replicate is done when its metrics_summary.csv exists, so
## re-running this script after an interruption picks up where it stopped --
## PROVIDED SIM_PROP_POS is still 0.2, which is why it is set here and not left
## to the shell.
##
## Usage:  nohup ./run_supplement.sh > results/main/supplement_pos20.log 2>&1 &
################################################################################
set -u
cd "$(dirname "$0")"   # manuscript_code: this script's own folder

export SIM_PROP_POS=0.2
TAG=09022026_pos20
REPS=50

echo "== supplementary Study A (20% positive / 80% negative) =="
echo "tag  : $TAG"
echo "reps : $REPS"
echo "start: $(date)"

################################################################################
## PREFLIGHT -- the silent crash mode.
################################################################################
## run_main.R decides a replicate is done by the existence of metrics_summary.csv,
## which the WORKER writes. The Julia MIMIX fit runs AFTER that, back in the
## driver. So a replicate interrupted between the two -- a crash, a power loss, a
## kill -- is left looking complete while carrying no MiMIX row, and a resume
## skips it forever. The cell then reports 50 replicates for thirteen methods and
## fewer for MiMIX, which no log anywhere would say out loud.
##
## This happened: the 2026-09-18 power loss left four such replicates in the
## pi = 0.02 cell. They are cheap to redo (the seed is the replicate index, so a
## re-run regenerates them identically) and impossible to repair in place, since
## the driver deletes the Y/X/Z handoff files Julia needs.
##
## The 10% ceiling is the safety catch. If MIMIX is systematically broken -- julia
## off PATH, a config moved -- then EVERY replicate looks like this one, and a
## blind cleanup would delete the entire study. Above the ceiling this reports and
## refuses rather than deleting.
MIMDIR="../simulation/src/mimix/SparseDOSSA-analysis-results/$TAG"
BROKEN=()
for d in results/main/$TAG/results_*(N); do
  if [[ -f "$d/metrics_summary.csv" && ! -f "$MIMDIR/$(basename $d)/beta_summary.csv" ]]; then
    BROKEN+=("$(basename $d)")
  fi
done
NDONE=$(find results/main/$TAG -name metrics_summary.csv 2>/dev/null | wc -l | tr -d ' ')
echo "preflight: $NDONE completed replicates, ${#BROKEN[@]} complete-but-missing-MIMIX"

if (( ${#BROKEN[@]} > 0 )); then
  if (( NDONE > 0 && ${#BROKEN[@]} * 10 > NDONE )); then
    echo "ABORT: ${#BROKEN[@]} of $NDONE replicates have no MIMIX output -- more than 10%."
    echo "       That is a broken MIMIX setup, not an interrupted run. Check that"
    echo "       julia is on PATH and that $MIMDIR is writable. Nothing deleted."
    exit 1
  fi
  echo "removing ${#BROKEN[@]} interrupted replicate(s) so they are re-run with MIMIX:"
  for b in "${BROKEN[@]}"; do
    echo "  $b"
    rm -rf "results/main/$TAG/$b" "$MIMDIR/$b"
  done
fi

for SW in overN overD overLambda overPi; do
  echo
  echo "=== sweep $SW  ($(date)) ==="
  Rscript run_main.R --sweep "$SW" --reps "$REPS" \
    > "results/main/${SW}_pos20.log" 2>&1
  echo "sweep $SW exited $? at $(date)"
  tail -5 "results/main/${SW}_pos20.log"

  ## Aggregate after EVERY sweep rather than only at the end, so each curve is
  ## readable while the rest is still running. aggregate_main.R skips any sweep
  ## with fewer than two settings, so a partial study emits no single-point
  ## figures. Outputs are suffixed "_pos20" (from the tag, via scenario_sfx_of),
  ## so they cannot overwrite the main scenario's figures.
  Rscript aggregate_main.R "$TAG" > "results/main/aggregate_${SW}_pos20.log" 2>&1 \
    || echo "AGGREGATE FAILED after $SW"
done

echo
echo "=== final figures ($(date)) ==="
Rscript plot_studyA_combined.R "$TAG" || echo "STUDYA PLOT FAILED"
Rscript plot_runtime.R "$TAG"         || echo "TIMING PLOT FAILED"

echo
echo "=== outputs ==="
find "results/main/$TAG" -name metrics_summary.csv | wc -l | \
  awk '{print "completed replicates: "$1" / 1200"}'
ls -la "results/main/$TAG"/*.csv 2>/dev/null
echo "done at $(date)"
