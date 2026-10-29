# manuscript_code

Code and raw inputs that recreate the latest manuscript results, the September 2026
rebuild, from scratch:

- `manuscript_tables/`: T1–T5, S1, S3, S4 and `manuscript_tables_09012026.tex`
- `manuscript_figures/`: the nine PDFs
- `simulation_rerun/results/`: Study A, tags `09022026` (main) and `09022026_pos20` (supplement)

Every Sparse SSRV fit here uses the scalar Laplace scale draw (`scale_draw = "shared"`: one
N(0, τ²) value per posterior draw, added to every taxon). That is the default in
`sparse_SSRV.R` and in `simulation_rerun/src/sparse_SSRV_rerun.R`, and nothing in this folder
overrides it.

The layout mirrors the original `sparse_DA` project, so every relative path in the scripts
still works. All outputs land inside this folder.

## Quick start

```bash
cd manuscript_code
bash run/01_fit_real_data.sh          # six real-data fits             ~6 h (selex alone ~5 h)
bash run/02_build_tables_figures.sh   # tables + figures from stage 1  ~5 min
bash run/03_ckf_simulation.sh         # T3 / S4 (pi = 0.20 and 0.05)   ~14 h
zsh  run/04_studyA_simulation.sh      # Study A, both scenarios        ~6 days, resumable
```

Logs go to `logs/`. Times are from the September runs on a 16 GB laptop.

## What makes what

| Output | Script | Stage |
|---|---|---|
| T1, T2, T4_1–T4_5, T5 (three files), S1, S3, `manuscript_tables_09012026.tex` | `manuscript_tables.R`, which sources `task_common.R`, Tasks 1–3, `2022_Jin…/jin_DA_biological_interpretation.R`, `jin_scale_uncertainty.R` | 2 |
| T3 and S4 (`_pi05` versions too) | `ckf_simulation_09012026.R`, which ends by sourcing `ckf_tables_09012026.R` | 3 |
| `lfc_diagnose_*.pdf` (7) | `Task 6 LFC distribution diagnostic.R` | 2 |
| `stammler_all_pairs_methods.pdf` | `Task 4 stammler all pairs figure.R` | 2 |
| `stammler_all_pairs_ssrv_variants.pdf` | `Task 5 stammler all pairs variant figure.R` | 2 |
| `simulation_rerun/results/main/<tag>/` (replicates, `aggregated_metrics.csv`, `scale_choice_summary.csv`, `ssrv_variants.csv`, `theta_perp_coverage.csv`, `runtime_*.csv`) | `run_main.R` (runs `src/worker_main.R` per replicate, then MiMIX), then `aggregate_main.R` | 4 |
| `simulation_rerun/results/figures/`: `main_*`, `decomp_*`, `scale_choice_*` | `aggregate_main.R` | 4 |
| `simulation_rerun/results/figures/`: `studyA_*`, `runtime_*` | `plot_studyA_combined.R`, `plot_runtime.R` | 4 |

Stage 1 writes each study's fits under the `09012026` stamp, plus the record files
(`scale_model_record.rds`, `sparse_arm_sig.rds`, `theta_perp_record.*`, `ckf_record.*`) that
stage 2 reads. The stage 1 scripts, one per study:

| Study | Script |
|---|---|
| Prasse et al. (2015) | `2015_prasse_appliedEnvMicrobiology_wetlands/sparse_DA_09012026_rawVSqseq_mean5.R` |
| Stämmler et al. (2016) | `2016_stammler_microbiome_micehuman/sparse_DA_09012026.R` |
| Zemb et al. (2020) | `2020_zemb_microOpen_spike/sparse_DA_09012026.R` |
| McMurrough et al. (2014), selex | `selex/selex 09012026.R` |
| Jin et al. (2022) | `2022_Jin_natureComm_technicalReplicates/sparse_DA_09012026.R` |
| Vandeputte et al. (2017) | `Vandeputte/vandeputte_09012026.R` |

In stage 2, `real data results.R` runs first and writes the four `real_data_*.RData`
artifacts. `theta_perp_check.R` refits Sparse SSRV for the θ⊥ recovery table (T5) and Task 6.

## Requirements

**R 4.4.3.** These are the versions installed when the September results were made.

- **GitHub:** LinDA 0.2.0 (zhouhj1994/LinDA), MRCKF 0.1.0 (aasrinivasan/MRCKF),
  QSeq 0.0.0.9000 (djeppschmidt/QSeq), SparseDOSSA2 0.99.2 (biobakery/SparseDOSSA2),
  driver 0.1.1 (jsilve24/driver)
- **Bioconductor:** ALDEx2 1.36.0, ANCOMBC 2.6.0, DESeq2 1.44.0, S4Vectors 0.42.1,
  baySeq 2.38.0, edgeR 4.2.2, limma 3.60.6, phyloseq 1.48.0
- **CRAN:** compositions 2.0.8, cowplot 1.1.3, data.table 1.16.4, directlabels 2024.1.21,
  dirmult 0.1.3.5, doParallel 1.0.17, dplyr 1.1.4, foreach 1.5.2, gghighlight 0.4.1,
  ggnewscale 0.5.2, ggpattern 1.1.4, ggplot2 3.5.1, ggrepel 0.9.6, knockoff 0.3.6,
  LaplacesDemon 16.1.6, latex2exp 0.9.6, magrittr 2.0.3, modeest 2.4.0, numDeriv 2016.8.1.1,
  patchwork 1.3.0, plyr 1.8.9, purrr 1.0.4, readr 2.1.5, readxl 1.4.3, stringi 1.8.4,
  stringr 1.5.1, tibble 3.2.1, tidyverse 2.0.0, vegan 2.6.10, viridis 0.6.5

`fido` is referenced only by an unused function in `simulation/src/01a_main_functions.R`
and is not needed.

**Julia 1.10, for MiMIX.** `simulation/src/mimix/julia_env/` holds the exact Julia
environment the September runs used: its `Project.toml` and `Manifest.toml`, with
`MicrobiomeMixedModels` pointed at the copy in `simulation/src/mimix/MicrobiomeMixedModels.jl`.
The drivers export `JULIA_PROJECT` to it. On a new machine, set it up once:

```bash
JULIA_PROJECT=simulation/src/mimix/julia_env julia -e 'using Pkg; Pkg.instantiate(); Pkg.build("RCall")'
```

## How the scripts find this folder

Every R script outside `simulation_rerun/` finds the project root by walking up from the
working directory to the marker file `.manuscript_code_root`. So run them from inside this
folder, as the drivers do. A script run from outside stops with an error rather than falling
back to `~/Research/sparse_DA`, so a run here cannot overwrite the original project.

`simulation_rerun/` already found its root by walking up to `sparse_SSRV.R` and is unchanged.

The only changes from the original files are those root-finding lines. A few related path
fixes come with them:

- the Zemb and Stämmler data-folder lookups;
- two hardcoded home paths at the end of the Jin script;
- the `cd` line in `simulation_rerun/run_supplement.sh`;
- the MiMIX path in the Julia `Manifest.toml`.

**The Jin study folder is renamed.** Here it is `2022_Jin_natureComm_technicalReplicates`,
because Jin et al. was published in Nature Communications in 2022. The original project
still calls it `2024_Jin_natureComm_technicalReplicates`. Every path and comment naming the
folder was updated to match: eight R files plus the stage 1 driver. Nothing else about the
study changed. Its dataset label in the tables was already `Jin2022`.

That is 24 edited files in all.

## Verified on 2026-10-29

- **Stage 2.** Rebuilt in a scratch copy from the September fit outputs.
  - All 17 `.tex` files in `manuscript_tables/` came out identical apart from their
    generation-timestamp lines. That includes T3/S4, rebuilt from the saved CKF outputs.
  - All 9 PDFs in `manuscript_figures/` are identical once PDF creation dates are stripped.
  - `theta_perp_check.R` reproduced all 81 current rows of `theta_perp_check.csv` exactly.
- **Stage 1.** The Zemb fit, rerun from the raw files, reproduced `res_tpfp_09012026.RData`
  exactly (every method's TP/FP, MiMIX and MRCKF included). It also reproduced
  `results_09012026.RData`, `sparse_sens_09012026.RData`, `zemb_mrckf_out.RData` and its
  scale-model record. MiMIX ran through the pinned Julia environment.
- **Stage 4.** One replicate, `results_100_300_10000_0_2_1` (the λ = 10⁴ point of
  overLambda), was rerun with the current worker. It is byte-identical to the stored Sep 5
  copy in all 33 result files: `metrics_summary.csv`, every method's output, the three SSRV
  arm tables and `theta_perp.csv`. The only exceptions:
  - `execution time.csv`, which holds timings;
  - `replicate_info.csv`, which gained four scenario columns on Sep 16
    (`prop_pos = 0.8`, …) with every shared column equal;
  - the MiMIX handoff files, which `run_main.R` handles, not the worker.

  So the Sep 16 edits leave the main scenario unchanged.
- **Not rerun:** the other five real-data fits, the CKF simulation itself, and the rest of
  Study A. These take hours to days.

## Provenance notes

1. **Files edited after the run they reproduce.** The copies are the current versions.
   - `2022_Jin…/sparse_DA_09012026.R` was edited on 2026-09-29; its fits are from 2026-09-01.
     Its sibling `sparse_DA_09082026.R` is a mean-filter variant made on 2026-09-08. Apart
     from the variant's own edits, the two differ only in a 20-line block at the very end:
     it re-reads `jin_sig_summary_flat_09012026.csv` and writes `jin_summary_report.csv`
     after every result is saved. That block is most likely the later edit.
   - `ckf_simulation_09012026.R` and `ckf_tables_09012026.R` were edited 2026-09-10 10:02.
     The π = 0.20 simulation finished 2026-09-09 23:46 under the earlier version. The
     π = 0.05 run started at 10:07, after the edit.
   - In `simulation_rerun/`, these files were edited on 2026-09-16 to add the sign-split
     scenario: `config.R`, `run_main.R`, `src/worker_main.R`, `src/simulate_one.R`,
     `aggregate_main.R`, `plot_studyA_combined.R`, `plot_runtime.R` and `README.md`.
     That was after the `09022026` replicates (Sep 2–5) and figures (Sep 9), and before
     the `_pos20` run (Sep 16–19). See the Stage 4 check above.
2. **Left out on purpose.** These ran in the same round, but no manuscript table or figure
   reads them:
   - Prasse `sparse_DA_09012026_mean5.R`, `sparse_DA_09012026_subestuary_family_mean5.R`
     and `prasse_DA_nested_contrasts.R`.
   - Jin `sparse_DA_09082026.R`. It writes under separate `_mean5` record keys, so it never
     touches what the tables read.
   - The Study A pilot (`run_pilot.R`, `aggregate_pilot.R`, `aggregate_sigma.R`).
   - `finish_overN.sh`, which targets the older 08292026 tag.
   - `fix_aldex2_sig.R` and `redecide_scale_rule.R`.
   - `tests/test_edits_additive.R`, which needs the original `sparse_SSRV_03122025.R`.
3. **`2016_stammler…/stammler_mrckf_out_filter5.RDS`** is a 2026-08-15 file. The Stämmler
   script reads it only in its last block, for the legacy plot
   `stammler_all_pairs_results2.pdf`, which is not a manuscript figure. It is included so
   the script runs to the end.
4. **Records accumulate.** Each record file and `theta_perp_check.csv` replaces the rows for
   its own (dataset, pair) key and keeps the rest. A fresh run here therefore holds only
   this round's rows, while the original project's copies also hold older keys.
5. **Seeds.** Every fit script sets its seed, and Study A seeds each replicate by its index.
   On this machine, with these package versions, the reruns above reproduce exactly.

## Raw inputs

| Dataset | Files |
|---|---|
| Prasse et al. (2015) | `2015_prasse…/Data/TFW_comdat.RDS`, `TFW_metadata.csv` |
| Stämmler et al. (2016) | `2016_stammler…/40168_2016_175_MOESM{1,7,9,11}_ESM.*` |
| Zemb et al. (2020) | `2020_zemb…/zemb_counts.csv`, `zemb_metadata.csv`, `zemb_qPCR.csv` |
| McMurrough et al. (2014) | `data("selex")` from ALDEx2; `selex/selex-truth.txt` |
| Jin et al. (2022) | `2022_Jin…/Supplementary Data {3,4,5,7,8,9}.xlsx` |
| Vandeputte et al. (2017) | `Vandeputte/data/41586_2017_BFnature24460_MOESM10_ESM.csv`, `OTU_nochim.csv`, `taxa_assignments.csv` |
| Study A, CKF simulation | generated by SparseDOSSA2 (Stool template) |
