#!/bin/bash
# Stage 3 -- the CKF vs Sparse SSRV sign-balance simulation behind T3 and S4.
#
# Two designs, each written under its own tag so neither overwrites the other:
#   pi = 0.20 (the default)          -> T3_ckf_sign_balance.tex,      S4_ckf_sign_balance_arms.tex
#   pi = 0.05 (20 of 400 spiked taxa) -> T3_ckf_sign_balance_pi05.tex, S4_ckf_sign_balance_arms_pi05.tex
# ckf_simulation_09012026.R sources ckf_tables_09012026.R at the end, so each
# run writes its two tables into manuscript_tables/ itself. The 2026-09-10
# pi = 0.05 run took about 7 h on 5 cores.
#
# Outputs land in "variable selection on CLR/results with sparsity 0.20/" (the
# folder name is historical; both designs write there under their own tag).

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
LOG="$ROOT/logs"
mkdir -p "$LOG"

failures=0

Rscript ckf_simulation_09012026.R > "$LOG/ckf_simulation.log" 2>&1
rc=$?
echo "pi = 0.20: exit=$rc"
[ "$rc" -eq 0 ] || failures=1

CKF_PI=0.05 CKF_TAG=pi05_09012026 Rscript ckf_simulation_09012026.R > "$LOG/ckf_simulation_pi05.log" 2>&1
rc=$?
echo "pi = 0.05: exit=$rc"
[ "$rc" -eq 0 ] || failures=1

exit "$failures"   # 1 if either run failed
