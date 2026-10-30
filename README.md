# sparse-da-paper-code

Code for the paper **Sparse Bayesian Partially Identified Models for Microbiome
Data** by Won Gu, Francesca Chiaromonte, and Justin D. Silverman. The code was
written by Won Gu.

The repository holds the analysis scripts and input data for six real-data
studies and two simulation studies, as run in September 2026. It does not
include results, tables, or figures; you generate those by running the scripts.
The main Sparse SSRV implementation is `sparse_SSRV.R`. Study A uses its own
copy, `simulation_rerun/src/sparse_SSRV_rerun.R`.

## Setup

You need Linux or macOS with Bash, zsh, R, and Julia. The original runs used
R 4.4.3 and Julia 1.10. The R package versions used are listed in
[dependencies.tsv](dependencies.tsv). Other versions may give slightly
different results.

Where to install the R packages from:

- **GitHub:** [LinDA](https://github.com/zhouhj1994/LinDA),
  [MRCKF](https://github.com/aasrinivasan/MRCKF),
  [QSeq](https://github.com/djeppschmidt/QSeq),
  [SparseDOSSA2](https://github.com/biobakery/SparseDOSSA2), and
  [driver](https://github.com/jsilve24/driver).
- **Bioconductor:** ALDEx2, ANCOMBC, DESeq2, S4Vectors, baySeq, edgeR, limma,
  and phyloseq.
- **CRAN:** everything else in `dependencies.tsv`.

You do not need to install fido. It is only used by a function in
`simulation/src/01a_main_functions.R` that the analyses never call.

Run everything from the repository root. The R scripts find the root by looking
for the file `.manuscript_code_root`, so keep that file if you copy the
repository. (Study A looks for `sparse_SSRV.R` instead.)

```sh
# Check that every R script parses (does not run anything):
Rscript --vanilla run/check_environment.R --syntax-only

# Check that the R packages, bash, zsh, and julia are installed:
Rscript --vanilla run/check_environment.R

# Set up the Julia environment used by MiMIX (once per machine):
JULIA_PROJECT=simulation/src/mimix/julia_env julia -e 'using Pkg; Pkg.instantiate(); Pkg.build("RCall")'
```

The Julia environment in `simulation/src/mimix/julia_env/` pins the exact
package versions and points to the copy of MiMIX in this repository. The run
scripts set `JULIA_PROJECT` for you.

## Running the analyses

There are four stages. Run them in order:

```sh
bash run/01_fit_real_data.sh          # about 6 hours
bash run/02_build_tables_figures.sh   # about 5 minutes
bash run/03_ckf_simulation.sh         # about 14 hours
zsh run/04_studyA_simulation.sh       # about 6 days; can be restarted
```

Times are from the original runs on a 16 GB laptop. SELEX alone takes about
5 hours of stage 1. The simulations use several cores; check
`simulation_rerun/config.R` before running them on a shared machine.

**Stage 1** fits the six real-data studies. Results are saved with the date
stamp `09012026`, along with record files that stage 2 reads.

| Study | Script |
|---|---|
| Prasse et al. (2015) | `2015_prasse_appliedEnvMicrobiology_wetlands/sparse_DA_09012026_rawVSqseq_mean5.R` |
| Stämmler et al. (2016) | `2016_stammler_microbiome_micehuman/sparse_DA_09012026.R` |
| Zemb et al. (2020) | `2020_zemb_microOpen_spike/sparse_DA_09012026.R` |
| McMurrough et al. (2014), SELEX | `selex/selex 09012026.R` |
| Jin et al. (2022) | `2022_Jin_natureComm_technicalReplicates/sparse_DA_09012026.R` |
| Vandeputte et al. (2017) | `Vandeputte/vandeputte_09012026.R` |

**Stage 2** builds the tables and figures. It first runs `real data results.R`,
which combines the stage 1 results into four `real_data_*.RData` files. Next,
`theta_perp_check.R` refits Sparse SSRV to check how well it recovers scale
(used by table T5 and Task 6). Then it makes the tables and figures.

**Stage 3** runs the CKF simulation twice, with 20% and with 5% of taxa
changed, and writes tables T3 and S4.

**Stage 4** runs Study A twice: once with 80% of the changed taxa increasing
(the main scenario, tag `09022026`) and once with 20% increasing (the
supplement, tag `09022026_pos20`). See the
[Study A README](simulation_rerun/README.md) to run one part at a time.

Logs for stages 1–3 go to `logs/`. Study A logs go to
`simulation_rerun/results/main/`.

If a script in stages 1–3 fails, the stage runs the remaining scripts and then
exits with an error. (Stage 2 stops right away if `real data results.R` or
`theta_perp_check.R` fails, because everything after them depends on their
output.) Stage 4 is different: it prints a message
when something fails but still exits normally, so read its logs before trusting
the results. When it is restarted, the supplement script deletes any replicate
that finished without its MiMIX results and runs it again (see the top of
`simulation_rerun/run_supplement.sh`).

## Outputs

| Output | Made by | Stage |
|---|---|---|
| `manuscript_tables/`: T1, T2, T4_1–T4_5, T5 (three files), S1, S3, and `manuscript_tables_09012026.tex` | `manuscript_tables.R`, which runs Tasks 1–3 and the study summaries | 2 |
| `manuscript_tables/`: T3 and S4, plus `_pi05` versions | `ckf_simulation_09012026.R`, which runs `ckf_tables_09012026.R` at the end | 3 |
| `manuscript_figures/lfc_diagnose_*.pdf` (seven figures) | `Task 6 LFC distribution diagnostic.R` | 2 |
| `manuscript_figures/stammler_all_pairs_methods.pdf` | `Task 4 stammler all pairs figure.R` | 2 |
| `manuscript_figures/stammler_all_pairs_ssrv_variants.pdf` | `Task 5 stammler all pairs variant figure.R` | 2 |
| `simulation_rerun/results/main/<tag>/`: replicates and summary CSVs | `simulation_rerun/run_main.R`, then `aggregate_main.R` | 4 |
| `simulation_rerun/results/figures/`: `main_*`, `decomp_*`, `scale_choice_*` | `simulation_rerun/aggregate_main.R` | 4 |
| `simulation_rerun/results/figures/`: `studyA_*`, `runtime_*` | `simulation_rerun/plot_studyA_combined.R` and `plot_runtime.R` | 4 |

`manuscript_tables/` and `manuscript_figures/` start out empty. Git ignores all
generated files.

The input data and where each file came from are listed in
[THIRD_PARTY.md](THIRD_PARTY.md). The simulations generate their own data from
the SparseDOSSA2 Stool template.

`2016_stammler_microbiome_micehuman/stammler_mrckf_out_filter5.RDS` is an older
saved fit. Only the last block of the Stämmler script reads it, to draw a plot
that is not in the paper. It is included so the script runs to the end.

## Notes on rerunning

**Record files keep old rows.** Each time a study runs, it replaces its own rows
in the shared record files and leaves other rows alone. If you rerun only some
studies, the records will mix results from different runs. For a clean rebuild,
delete old results first. (`theta_perp_check.csv` is rebuilt from scratch each
time.)

**Seeds.** Every fit script sets a random seed, and Study A uses the replicate
number as its seed. Even so, different package versions can change the numbers.

**Scale draw and scale rule.** By default, Sparse SSRV uses
`scale_draw = "shared"`: for each posterior draw, one random number is added to
every taxon. It picks between its two versions (arms) with
`scale_rule = "width"`, which chooses the arm whose 95% credible intervals have
the larger median width.

**What has been checked.** Before release, the reporting stage, the Zemb fit,
and one Study A replicate were rerun and matched the September results. For this
release, we only checked that the scripts parse and that the stage 1–3 scripts
report failures correctly. We did not rerun the full analyses.

## License and citation

The code is released under the [MIT License](LICENSE). This includes a short
section adapted from fido, which its author, Justin D. Silverman, has allowed
us to release under MIT. The bundled copy of MiMIX keeps its own
[MIT License](simulation/src/mimix/LICENSE). The datasets are not covered by
our license; each keeps its original terms. See [THIRD_PARTY.md](THIRD_PARTY.md)
for details.

If you use this code, please cite the paper:

> Won Gu, Francesca Chiaromonte, and Justin D. Silverman.
> *Sparse Bayesian Partially Identified Models for Microbiome Data*.

Citation details are also in [CITATION.cff](CITATION.cff). We will add the DOI
once the paper is published.
