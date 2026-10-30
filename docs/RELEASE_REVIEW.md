# Release preparation review — 2026-10-29

This is a source and packaging review, not a full statistical audit or numerical
reproduction. No inference code was changed. The manuscript examined was the
local `overleaf-aoas-resubmission` working tree; its HEAD was
`daab13506db92437919d351774251c07998e5ad2`. Working-tree content may differ from
that commit. Manuscript paths below are relative to that repository.

## Issues requiring author decisions

1. **Hybrid rule.** `main.tex:407` describes choosing the larger variance of
   theta-perp. `sparse_SSRV.R:560` defaults to `scale_rule="width"` and lines
   734–778 choose the larger median taxon-wise 95% credible interval width of
   the delivered theta posterior. These are different selection rules.
2. **Laplace draw.** `supplementary_body.tex:1712–1713` describes independent
   Gaussian noise across taxa. `sparse_SSRV.R:561–568,675–676` uses one shared
   scalar draw by default, consistent with the main text's scalar scale model.
   The duplicated simulation implementation also defaults to shared draws.
3. **CKF sequencing depth.** `supplementary_body.tex:152` gives 20,000,000;
   `ckf_simulation_09012026.R:73–75,145` uses 400 × 10,000 = 4,000,000.
4. **FDR averaging.** `ckf_simulation_09012026.R:43–49` computes
   `FP/(TP+FP)`, producing NaN for zero discoveries. In
   `ckf_tables_09012026.R:49–60`, Sparse SSRV and its arms drop these values
   when averaging. Lines 62–66 replace CKF NaNs with zero before averaging.
   Thus the table can compare different estimands: for Sparse SSRV, averaging
   excludes zero-discovery replicates; for CKF, those replicates contribute
   zero. The size of the effect has not been quantified because saved fits
   are not included. Distinguish failed fits from valid zero-discovery fits
   before correcting and regenerating the tables.
5. **Missing figure generator.** `supplementary_body.tex:1846` references
   `real/sensitivity_3combined 09012026 mean5.pdf`. No matching filename or
   generator was found in the supplied bundle. Locate its source before
   claiming complete reproduction of all manuscript figures.
6. **Rights and provenance.** The authors' original code now has an
   [MIT License](../LICENSE). Complete [the third-party inventory](../THIRD_PARTY.md),
   including the MiMIX source revision and unresolved dataset terms. MiMIX's
   upstream MIT notice has been restored, and Justin D. Silverman has authorized
   the included fido-derived excerpt under MIT. The supplied bundle had no
   original license file.
7. **R environment.** Supplied versions are historical records, with no
   `renv.lock` or exact GitHub revisions. Do not create a lockfile from this
   review machine and call it the original analysis environment. Capture a
   resolved environment from a validated reproduction.

## Manuscript output handoff

Names differ between generated artifacts and manuscript inputs. These are
candidate correspondences from source inspection, not numerically verified
equivalences. Inspect and compare contents before copying or renaming.

| Generated table stem (`manuscript_tables/`) | Manuscript table name |
|---|---|
| `T1_dataset_summary` | `real_data_summary` |
| `T2_da_results_hybrid` | `real_data_results_hybrid` |
| `T5_theta_perp_recovery_hybrid` | `real_theta_perp_recovery_main` |
| `T5_theta_perp_recovery_stammler` | `stammler_all_pairs_theta_perp_recovery` |
| `S1_ssrv_variants_tpfp` | `real_scale_models_bws` |

See `manuscript_tables.R:743–757` for generated names. The supplementary Study A
plot referenced at `supplementary_body.tex:1769` uses `_supp.png`, whereas
`simulation_rerun/plot_studyA_combined.R` generates a `_pos20` scenario suffix.
The existing manuscript outputs were not overwritten during this review.

## Operational findings

- The bundle is about 15 MB; its largest file is 8,111,375 bytes. No nested
  Git repositories were found.
- The original README's opening description refers to outputs that are not
  bundled. Only raw inputs, source code, and a legacy Stämmler RDS fit are
  present. The new README makes this explicit.
- Stage 1–3 wrappers previously ended successfully even after failed R jobs.
  They now aggregate failures while preserving job order and continuation of
  independent work. Stage 2 still stops if either prerequisite fails.
- Stage 4 and its supplementary runner still mask some failures. Inspect logs
  and method completeness before accepting simulation aggregates. The inherited
  supplementary runner removes a limited number of incomplete replicate
  directories to rerun them; read its preflight before resuming an old run.
  Its comments also describe a completion marker that can precede MiMIX output.
- Prasse's script retains a guarded cluster library path at line 115. It has
  no effect if that directory does not exist. Commented cluster paths remain
  as source history.
- `.DS_Store` files and redundant `.ipynb_checkpoints` were omitted from the
  release copy. Original file hashes remain in `source-manifest.tsv`.
- GitHub CLI authentication reports failure for the configured account.
  Authentication and the destination account/organization must be resolved
  before creating the remote and pushing.

## Checks performed here

- Parsed all 44 supplied R files successfully with R 4.6.1; the new environment
  checker also parses successfully (45 total).
- Validated Bash syntax for stages 1–3 and zsh syntax for stage 4 and the
  supplementary runner.
- Exercised stages 1–3 with mock Rscript commands: all-success, injected
  failures, continued independent jobs, and stage 2 prerequisite abort.
- Checked the local environment against the recorded dependencies. Many
  required packages are absent, installed versions differ, and Julia is
  unavailable. Full analyses were not run.
- Reviewed source for credential patterns, absolute host paths, and nested
  repositories; no credential-like source findings were identified. This is
  not an audit of participant-level dataset privacy or consent.
- Compared source hashes to confirm the supplied R, Julia, and raw input files
  remain unchanged. The original instructions are preserved separately in
  `STUDENT_README.md`; their earlier reproduction claims were not reverified.

## GitHub handoff

After resolving the author decisions above, authenticate `gh`, choose the
destination, and create the remote from this directory. For example, replacing
`OWNER` with the selected account or organization:

```sh
gh auth login
gh repo create OWNER/sparse-da-paper-code --public --source=. --remote=origin --push
```

This command publishes immediately. Review the staged/committed file list and
resolved licensing first. Add the repository URL to `CITATION.cff` and the
paper's code-availability statement, then tag the validated paper release.
