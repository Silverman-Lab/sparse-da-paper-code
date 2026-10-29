#!/bin/bash
# Stage 2 -- rebuild manuscript_tables/ (except T3/S4, see stage 3) and
# manuscript_figures/ from the stage 1 fits.
#
# Adapted from run_logs_09012026/run_reporting.sh in the original project.
#
# ORDER IS NOT ARBITRARY:
#   `real data results.R` writes the four cross-dataset artifacts
#   (real_data_scale_model / _results / _arm_features / _truth_sets .RData) that
#   task_common.R stops on if they are missing, so it goes first.
#
#   theta_perp_check.R refits Sparse SSRV for Zemb, Vandeputte, Jin and Prasse
#   (about 3 min) and writes theta_perp_check/, which T5 and Task 6 read, so it
#   must run before them.
#
#   manuscript_tables.R is the umbrella for the tables: it sources
#   jin_DA_biological_interpretation.R, task_common.R and Tasks 1-3. Tasks 4-6
#   are the figures and stand alone.
#
#   WITH_BOOTSTRAP=1 adds the original driver's second pass, which forces the
#   bootstrap arm and writes *_bootstrap.tex beside the hybrid tables. Off by
#   default because manuscript_tables/ holds only the hybrid tables.

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
LOG="$ROOT/logs"
mkdir -p "$LOG"
STATUS="$LOG/report_status.txt"
: > "$STATUS"

run() {   # run <label> <rscript args...>
  local label="$1"; shift
  local t0 t1 rc
  t0=$(date +%s)
  echo "[$(date '+%H:%M:%S')] START  $label" >> "$STATUS"
  Rscript "$@" > "$LOG/report_$label.log" 2>&1
  rc=$?
  t1=$(date +%s)
  echo "[$(date '+%H:%M:%S')] END    $label  exit=$rc  $((t1 - t0))s" >> "$STATUS"
  return $rc
}

if ! run real_data_results "real data results.R"; then
  echo "[$(date '+%H:%M:%S')] ABORT -- real data results.R failed; nothing downstream can run" >> "$STATUS"
  exit 1
fi
if ! run theta_perp_check "theta_perp_check.R"; then
  echo "[$(date '+%H:%M:%S')] ABORT -- theta_perp_check.R failed; T5 and Task 6 would be stale" >> "$STATUS"
  exit 1
fi

failures=0
run tables_hybrid    "manuscript_tables.R" || failures=1
if [ "${WITH_BOOTSTRAP:-0}" = 1 ]; then
  run tables_bootstrap -e 'ARM <- "bootstrap"; source("manuscript_tables.R")' || failures=1
fi

run fig4_stammler_pairs    "Task 4 stammler all pairs figure.R" || failures=1
run fig5_stammler_variants "Task 5 stammler all pairs variant figure.R" || failures=1
run fig6_lfc_diagnostic    "Task 6 LFC distribution diagnostic.R" || failures=1

if [ "$failures" -eq 0 ]; then
  echo "[$(date '+%H:%M:%S')] REPORTING DONE" >> "$STATUS"
else
  echo "[$(date '+%H:%M:%S')] REPORTING DONE WITH FAILURES" >> "$STATUS"
fi
exit "$failures"
