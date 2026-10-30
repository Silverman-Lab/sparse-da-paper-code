#!/bin/bash
# Stage 1 -- the six real-data fits of the 09012026 round.
#
# Adapted from run_logs_09012026/run_all_da.sh in the original project: same
# scripts and the same scheduling, minus the three Prasse runs that no
# manuscript table or figure reads (pooled mean5, subestuary_family, nested
# contrasts), plus Vandeputte, which was run on its own just before that batch.
#
# Measured on the 2026-09-01 run (16 GB laptop): Prasse 10 min, Stammler 23 min,
# Zemb 5 min, Vandeputte 3 min, Jin 30 min, selex 5 h. selex and Jin are the
# memory hogs (3.6 GB and 1.2 GB per subset for the resample array), so they run
# alone after the small studies finish.
#
# Exit status of every script goes to logs/fit_status.txt rather than aborting
# the batch: one study failing should not cost the others. The script still
# exits with status 1 at the end if any study failed.

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
export JULIA_PROJECT="$ROOT/simulation/src/mimix/julia_env"   # pinned MiMIX environment
LOG="$ROOT/logs"
mkdir -p "$LOG"
STATUS="$LOG/fit_status.txt"
: > "$STATUS"

run() {   # run <label> <script path>
  local label="$1" script="$2" t0 t1
  t0=$(date +%s)
  echo "[$(date '+%H:%M:%S')] START  $label" >> "$STATUS"
  Rscript "$script" > "$LOG/$label.log" 2>&1
  local rc=$?
  t1=$(date +%s)
  echo "[$(date '+%H:%M:%S')] END    $label  exit=$rc  $((t1 - t0))s" >> "$STATUS"
  return "$rc"
}

# ---- lane 1: Prasse, then Vandeputte ----------------------------------------
# Prasse stages MiMIX input in its own Data/ folder, so it must never overlap
# another Prasse script (none other is in this folder).
(
  lane_rc=0
  run prasse_rawVSqseq_mean5 "2015_prasse_appliedEnvMicrobiology_wetlands/sparse_DA_09012026_rawVSqseq_mean5.R" || lane_rc=1
  run vandeputte             "Vandeputte/vandeputte_09012026.R" || lane_rc=1
  exit "$lane_rc"
) &
LANE1=$!

# ---- lane 2: Stammler then Zemb ---------------------------------------------
(
  lane_rc=0
  run stammler "2016_stammler_microbiome_micehuman/sparse_DA_09012026.R" || lane_rc=1
  run zemb     "2020_zemb_microOpen_spike/sparse_DA_09012026.R" || lane_rc=1
  exit "$lane_rc"
) &
LANE2=$!

failures=0
wait "$LANE1" || failures=1
wait "$LANE2" || failures=1

# ---- then the two memory-hungry studies, alone ------------------------------
run selex "selex/selex 09012026.R" || failures=1
run jin   "2022_Jin_natureComm_technicalReplicates/sparse_DA_09012026.R" || failures=1

if [ "$failures" -eq 0 ]; then
  echo "[$(date '+%H:%M:%S')] ALL DONE" >> "$STATUS"
else
  echo "[$(date '+%H:%M:%S')] ALL DONE WITH FAILURES" >> "$STATUS"
fi
exit "$failures"
