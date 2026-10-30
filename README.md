# sparse-da-paper-code

Research code accompanying **Sparse Bayesian Partially Identified Models for
Microbiome Data**, by Won Gu, Francesca Chiaromonte, and Justin D. Silverman.
The code was written by Won Gu.

This repository contains the September 2026 analysis scripts and input data for
six empirical studies and two simulation analyses. Fitted results, tables, and
figures are generated locally and are not included. The main Sparse SSRV
implementation is in `sparse_SSRV.R`; Study A has a separate implementation in
`simulation_rerun/src/sparse_SSRV_rerun.R`.

## Setup

Use Linux or macOS with Bash, zsh, R, and Julia on `PATH`. The original runs used
R 4.4.3 and Julia 1.10. Recorded R package versions are listed in
[dependencies.tsv](dependencies.tsv). Install the dependencies before running
analyses; these records are historical versions, not a resolved R lockfile.

The recorded installation sources are:

- **GitHub:** [LinDA](https://github.com/zhouhj1994/LinDA),
  [MRCKF](https://github.com/aasrinivasan/MRCKF),
  [QSeq](https://github.com/djeppschmidt/QSeq),
  [SparseDOSSA2](https://github.com/biobakery/SparseDOSSA2), and
  [driver](https://github.com/jsilve24/driver).
- **Bioconductor:** ALDEx2, ANCOMBC, DESeq2, S4Vectors, baySeq, edgeR, limma,
  and phyloseq.
- **CRAN:** the remaining packages in `dependencies.tsv`.

The fido reference in `simulation/src/01a_main_functions.R` belongs to an unused
function and does not require installing fido for these analyses.

Run commands from the repository root. Keep `.manuscript_code_root` when
copying the repository: R scripts use it to locate the project. Study A locates
its root using `sparse_SSRV.R`.

```sh
# Parse the R scripts without running analyses or loading analysis packages:
Rscript --vanilla run/check_environment.R --syntax-only

# Check installed R packages and required executables:
Rscript --vanilla run/check_environment.R

# Instantiate the supplied Julia environment and build its R bridge:
JULIA_PROJECT=simulation/src/mimix/julia_env julia -e 'using Pkg; Pkg.instantiate(); Pkg.build("RCall")'
```

The Julia environment includes a manifest and a local copy of MiMIX. The run
scripts select this environment through `JULIA_PROJECT`.

## Run the analyses

```sh
bash run/01_fit_real_data.sh          # approximately 6 hours
bash run/02_build_tables_figures.sh   # approximately 5 minutes, after stage 1
bash run/03_ckf_simulation.sh         # approximately 14 hours
zsh run/04_studyA_simulation.sh       # approximately 6 days, resumable
```

These durations were reported for a 16 GB laptop; SELEX accounts for about
five hours of stage 1. Simulations use multiple cores. Review
`simulation_rerun/config.R` before running on shared hardware.

Stage 1 fits the six empirical studies and writes results with the `09012026`
stamp, together with shared records used by stage 2:

| Study | Analysis script |
|---|---|
| Prasse et al. (2015) | `2015_prasse_appliedEnvMicrobiology_wetlands/sparse_DA_09012026_rawVSqseq_mean5.R` |
| Stämmler et al. (2016) | `2016_stammler_microbiome_micehuman/sparse_DA_09012026.R` |
| Zemb et al. (2020) | `2020_zemb_microOpen_spike/sparse_DA_09012026.R` |
| McMurrough et al. (2014), SELEX | `selex/selex 09012026.R` |
| Jin et al. (2022) | `2022_Jin_natureComm_technicalReplicates/sparse_DA_09012026.R` |
| Vandeputte et al. (2017) | `Vandeputte/vandeputte_09012026.R` |

Stage 2 first runs `real data results.R` to build the four cross-study
`real_data_*.RData` artifacts. It then runs `theta_perp_check.R`, which refits
Sparse SSRV for the scale-recovery tables and diagnostics, before generating
tables and figures. Stage 3 runs the CKF comparison at signal proportions
0.20 and 0.05 and generates its tables. Stage 4 runs both Study A scenarios: 80% of spiked
taxa increase in the main scenario (`09022026`), and 20% increase in the
supplementary scenario (`09022026_pos20`). See
[the Study A instructions](simulation_rerun/README.md) for individual sweeps and
configuration details.

Logs for stages 1–3 are in `logs/`; Study A logs are under
`simulation_rerun/results/main/`. Stages 1–3 return a failing exit status if an
analysis fails. Stage 4 can print failures and still exit successfully, so
inspect its logs and method completeness before accepting aggregates. Its
supplementary runner can remove incomplete replicate directories before
rerunning them; review its preflight before resuming existing results.

## Find the outputs

| Output | Generator | Stage |
|---|---|---|
| `manuscript_tables/`: T1, T2, T4_1–T4_5, T5 (three files), S1, S3, and `manuscript_tables_09012026.tex` | `manuscript_tables.R`, which sources the study summaries and Tasks 1–3 | 2 |
| `manuscript_tables/`: T3 and S4, including `_pi05` versions | `ckf_simulation_09012026.R` and `ckf_tables_09012026.R` | 3 |
| `manuscript_figures/lfc_diagnose_*.pdf` (seven figures) | `Task 6 LFC distribution diagnostic.R` | 2 |
| `manuscript_figures/stammler_all_pairs_methods.pdf` | `Task 4 stammler all pairs figure.R` | 2 |
| `manuscript_figures/stammler_all_pairs_ssrv_variants.pdf` | `Task 5 stammler all pairs variant figure.R` | 2 |
| `simulation_rerun/results/main/<tag>/`: replicates and aggregate CSVs | `run_main.R` and `aggregate_main.R` in `simulation_rerun/` | 4 |
| `simulation_rerun/results/figures/`: `main_*`, `decomp_*`, `scale_choice_*` | `simulation_rerun/aggregate_main.R` | 4 |
| `simulation_rerun/results/figures/`: `studyA_*`, `runtime_*` | `plot_studyA_combined.R` and `plot_runtime.R` in `simulation_rerun/` | 4 |

Generated results are ignored by Git. Input files and their sources are listed
in [THIRD_PARTY.md](THIRD_PARTY.md). Simulation inputs are generated with the
SparseDOSSA2 Stool template. The bundled legacy fit
`2016_stammler_microbiome_micehuman/stammler_mrckf_out_filter5.RDS` is read by a
final plotting block in the Stämmler script; that plot is not a manuscript figure.

## Rerunning and interpreting results

Shared record files replace entries for the current study/comparison key while
retaining other entries; `theta_perp_check.csv` is rebuilt from the stored
scale-recovery records. A partially rerun directory
can therefore contain records from several runs. Start from a clean results
state when a complete, consistent rebuild is needed. Study A uses the replicate
index as its seed and resumes existing replicates; fit scripts also set seeds.
Exact numerical agreement depends on the software environment.

Sparse SSRV uses `scale_draw = "shared"`: each posterior draw adds one scalar
Gaussian draw to every taxon. The default hybrid selection rule,
`scale_rule = "width"`, chooses the arm with the larger median taxon-wise 95%
credible-interval width of the delivered posterior.

The supplied reproduction notes reported matching September outputs for the
reporting stage, the Zemb fit, and one Study A worker replicate. Those checks
have not been independently repeated for this release. Release preparation
checked R and shell syntax and exercised the stage 1–3 wrappers with simulated
successes and failures; it did not rerun the full analyses.

The following manuscript/code differences remain unresolved:

- The manuscript describes hybrid selection by scale variance, while the code
  defaults to posterior interval width. Supplementary text describes independent
  Gaussian draws across taxa, while the code defaults to a shared scalar draw.
- The CKF simulation uses a sequencing depth of 4,000,000; the supplementary text
  gives 20,000,000.
- CKF tables handle zero-discovery replicates differently across methods:
  Sparse SSRV averages omit undefined false-discovery proportions, while CKF
  replaces NaNs with zero. Failed fits and valid zero-discovery replicates need
  to be distinguished before interpreting or changing these summaries.
- A generator for the manuscript figure
  `real/sensitivity_3combined 09012026 mean5.pdf` was not found in the supplied
  code. Generated table names also differ from manuscript input names, and
  Study A supplementary figures use `_pos20` rather than the manuscript's
  `_supp` suffix. Check contents before copying or renaming outputs.

## Licensing and citation


The authors’ original code and the included fido-derived excerpt are available
under the [MIT License](LICENSE). Justin D. Silverman, the author of that
excerpt, has authorized its inclusion under MIT. Bundled MiMIX code retains its
[upstream MIT License](simulation/src/mimix/LICENSE). Other third-party material
and datasets are excluded from the root license and retain their own terms.
Sources, applicable terms, and remaining provenance questions are documented in
[THIRD_PARTY.md](THIRD_PARTY.md).

If you use this code in your research, please cite the accompanying paper:

> Won Gu, Francesca Chiaromonte, and Justin D. Silverman.
> *Sparse Bayesian Partially Identified Models for Microbiome Data*.

See [CITATION.cff](CITATION.cff) for citation metadata. The paper DOI and
publication details will be added when available. Citation is requested as an
academic courtesy and is not an additional condition of the MIT License.
