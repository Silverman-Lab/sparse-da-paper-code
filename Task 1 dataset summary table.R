################################################################################
## TASK 1 -- Table 1: the characteristics of every real dataset.
##
##   Rscript "Task 1 dataset summary table.R"
##
## One row per analysed subset: scenario, study, sample material, dimensions,
## citation, sequencing method, how truth was defined, and the signal sparsity
## the study offers.
##
## Nothing is transcribed by hand. Dimensions come from the scale-model record of
## the run the results tables report, so they describe the same fit; signal
## sparsity is the size of that run's truth set over D, read from the truth sets
## `real data results.R` saves. Editorial columns -- scenario, sample, citation,
## sequencing method, truth definition -- are properties of the studies and are
## declared in DATASET_META, in task_common.R.
##
## OUTPUT   manuscript_tables/T1_dataset_summary.tex
##
## The same table manuscript_tables.R emits under this name; both build it from
## tab_dataset_summary() below, so the two cannot drift.
################################################################################

ROOT_CANDIDATES <- local({  # manuscript_code: the folder holding .manuscript_code_root (README.md)
  p <- normalizePath(getwd(), mustWork = FALSE)
  while (!file.exists(file.path(p, ".manuscript_code_root"))) {
    if (dirname(p) == p) stop("run from inside the manuscript_code folder (see README.md)")
    p <- dirname(p)
  }
  p
})
ROOT <- ROOT_CANDIDATES[dir.exists(ROOT_CANDIDATES)][1]
if (is.na(ROOT)) stop("no project root found")
source(file.path(ROOT, "task_common.R"))

check_results_freshness()
check_jin_gold()

################################################################################
## T1 -- summary of all datasets
################################################################################
tab_dataset_summary <- function() {
  ## Expand first, render second: a dataset with `subsets` contributes one row
  ## per subset, and the scenario multirow has to span ROWS, not datasets.
  specs <- list()
  for (k in names(DATASET_META)) {
    m <- DATASET_META[[k]]
    if (is.null(m$subsets)) {
      specs[[length(specs) + 1]] <- list(key = k, meta = m, pair = m$pair,
                                         sample = m$sample)
    } else {
      for (sub in m$subsets) {
        specs[[length(specs) + 1]] <- list(key = k, meta = m, pair = sub$pair,
                                           sample = sub$sample)
      }
    }
  }
  scen <- vapply(specs, function(z) z$meta$scenario, character(1))

  rows <- character(0); seen <- character(0)
  for (i in seq_along(specs)) {
    z <- specs[[i]]; m <- z$meta
    d <- dims_for(z$key, z$pair)
    first <- !(m$scenario %in% seen)
    if (first && length(seen)) rows <- c(rows, "        \\midrule")
    sc <- if (first) scenario_cell(m$scenario, sum(scen == m$scenario)) else ""
    seen <- c(seen, m$scenario)

    ## Sparsity is truths / D. Withheld where the study's reference cannot
    ## support one -- see check_jin_gold().
    nt <- n_truth_of(z$key)
    sparsity <- if (z$key %in% SPARSITY_WITHHELD || is.na(nt)) DASH
                else formatC(nt / d$D, format = "f", digits = 3)

    ## Dataset and Sample break before their trailing parenthetical, so each sets
    ## as two lines -- "McMurrough et al" / "(2014)", "Pig" / "(fecal)".
    rows <- c(rows, row_tex(sc, stack_paren(m$label),
                            stack_paren(z$sample),
                            sprintf("%d $\\times$ %d", d$D, d$N),
                            paste0("\\cite{", m$cite, "}"),
                            m$seq, m$truth_def, sparsity))
  }
  c("\\begin{table}[!tbp]",
    "    \\centering",
    "    \\resizebox{\\textwidth}{!}{%",
    "    \\begin{tabular}{c c c c c c c c}",
    "        \\toprule",
    "        Scenario & Dataset & Sample & \\shortstack{Dimension \\\\ (D$\\times$N)} & Reference & \\shortstack{Sequencing \\\\ method} & Truth Definition & \\shortstack{Signal Sparsity \\\\ (Truths / D)} \\\\",
    "        \\cmidrule(r){1-1}\\cmidrule(lr){2-2}\\cmidrule(lr){3-3}\\cmidrule(lr){4-4}\\cmidrule(lr){5-5}\\cmidrule(lr){6-6}\\cmidrule(lr){7-7}\\cmidrule(l){8-8}",
    rows,
    "        \\bottomrule",
    "    \\end{tabular}%",
    "    }",
    paste0("    \\caption{The real datasets used to benchmark the methods. Dimensions are read from the ",
           "scale-model record of the run each row reports, so they describe the fit the results tables ",
           "score. Signal sparsity is the size of the study's truth set divided by the number of features ",
           "retained after filtering. \\citet{jin2022quantsinglecell} is reported on its pooled cOTU ",
           "run, the two cecal locations analysed together: cOTUs are the representation BarBIQ ",
           "quantifies in absolute terms, and pooling is what gives the measured-scale reference enough ",
           "power to call any feature at an adjusted significance level -- it calls one pooled and none ",
           "at either location alone. That study's ASV and OTU tables are conventional 16S gene-copy ",
           "measurements with no measured scale, so they carry no truth set of the kind this column ",
           "counts; they and the two locations separately are reported in ",
           "Table~\\ref{tab:jin_dimensions}.}"),
    "    \\label{tab:real_studies_summary}",
    "\\end{table}")
}

################################################################################
## EMIT
################################################################################
T1 <- tab_dataset_summary()
write_table(T1, "T1_dataset_summary")

## Echo when run non-interactively, so `Rscript "Task 1 ..."` prints paste-ready
## LaTeX. Silent when sourced by manuscript_tables.R, which prints its own bundle.
if (!interactive() && !task_quiet()) cat(T1, sep = "\n")
