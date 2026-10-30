# sparse-da-paper-code

Research code accompanying **Sparse Bayesian Partially Identified Models for
Microbiome Data**, by Won Gu, Francesca Chiaromonte, and Justin D. Silverman.

The code was written by Won Gu.

This repository contains the supplied September 2026 analysis scripts and input
data. Fitted results, manuscript tables, and figures must be generated locally.
The analysis code has been preserved while preparing this repository for release.
See [the release review](docs/RELEASE_REVIEW.md) for unresolved manuscript/code
discrepancies and publication requirements.

## Reproduce the analyses

Use Linux or macOS with Bash, zsh, R, and Julia on `PATH`. The original runs used
R 4.4.3 and Julia 1.10. R package versions and installation sources are recorded
in [Won’s supplied instructions](docs/STUDENT_README.md#requirements) and
[dependencies.tsv](dependencies.tsv). These are historical version records,
not a resolved R lockfile. Install those dependencies before running analyses.
The Julia environment includes a manifest and a local copy of MiMIX.

From the repository root:

```sh
# Parse scripts without running any analysis or loading analysis packages:
Rscript --vanilla run/check_environment.R --syntax-only

# Check installed R dependencies and available executables:
Rscript --vanilla run/check_environment.R

# Instantiate Julia dependencies and build the R bridge:
JULIA_PROJECT=simulation/src/mimix/julia_env julia -e 'using Pkg; Pkg.instantiate(); Pkg.build("RCall")'

bash run/01_fit_real_data.sh          # approximately 6 hours
bash run/02_build_tables_figures.sh   # approximately 5 minutes, after stage 1
bash run/03_ckf_simulation.sh         # approximately 14 hours
zsh run/04_studyA_simulation.sh       # approximately 6 days, resumable
```

Durations were reported for a 16 GB laptop. Simulations use multiple cores;
review `simulation_rerun/config.R` before running on shared hardware. Logs go
to `logs/` and simulation result directories. Stages 1–3 return a failing exit
status if any analysis fails. Stage 4 still requires inspection of its logs;
its inherited scripts can print failures and exit successfully.

The `.manuscript_code_root` marker file is required for path discovery. Keep it
when copying the repository. Generated results are ignored by Git; the supplied
raw data and the legacy Stämmler RDS input are included explicitly.

## Contents and provenance

- `sparse_SSRV.R`: core Sparse SSRV implementation.
- Study directories: real-data analyses and supplied inputs.
- `manuscript_tables.R`, `Task *.R`: reporting scripts.
- `simulation_rerun/`: Study A drivers and aggregation.
- `ckf_simulation_09012026.R`: sign-balance simulations.
- `simulation/src/mimix/`: Julia MiMIX code and environment.
- [Original instructions](docs/STUDENT_README.md): detailed output mapping,
  data inventory, and earlier reproduction claims supplied with the code.
- [Source manifest](docs/source-manifest.tsv): SHA-256 hashes of the supplied
  files before release preparation, excluding `.DS_Store`.

Earlier successful reruns described in the original instructions have not been
independently repeated as part of this repository preparation. Full analyses
take hours to days. Current validation and limitations are recorded in the
[release review](docs/RELEASE_REVIEW.md).

## Licensing and citation

A license for the authors’ code has not yet been selected. Included third-party
code and datasets retain their own terms; their redistribution terms need to be
documented before public release. See [THIRD_PARTY.md](THIRD_PARTY.md).
Use [CITATION.cff](CITATION.cff) for the manuscript title and authors; add the
paper DOI and release identifier when available.
