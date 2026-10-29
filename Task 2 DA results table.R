################################################################################
## TASK 2 -- Table 2: differential abundance results on the real datasets.
##
##   Rscript "Task 2 DA results table.R"
##
## The MAIN results table. Each cell is the number of false positives over the
## number of true positives in parentheses, scored against the study's truth set.
## Sparse SSRV appears once, under the HYBRID rule at the DEFAULT bandwidth
## (nrd0) -- the method as the manuscript proposes it -- beside every competing
## method. All nine Sparse SSRV variants are in the supplementary tables that
## Task 3 builds.
##
## \citet{jin2022quantsinglecell} is included, on its POOLED cOTU run: pooling
## the two cecal locations is what gives its measured-scale reference enough
## power to call a feature, so that run can be scored where neither location can
## on its own. Its ASV and OTU representations carry no measured scale at all and
## stay in the Jin subset tables.
##
## The companion table that reported the two scale models separately is GONE. It
## showed two of the nine Sparse SSRV settings; the variant table Task 3 builds
## shows all nine, including those two, so keeping both meant maintaining one
## table that was a strict subset of another.
##
## SOURCE   real_data_results.RData -- the datasets x methods matrix of "FP/(TP)"
##          cells that `real data results.R` builds and saves. This file does not
##          refit anything, so if that script has not been re-run since the last
##          analysis change these numbers are stale; check_results_freshness()
##          compares timestamps and warns.
##
## OUTPUT   manuscript_tables/T2_da_results_hybrid.tex
##
## The same table manuscript_tables.R emits under that name; both build it from
## tab_da_results() below, so the two cannot drift.
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
check_truth_sizes()
HYBRID_MISMATCH    <- check_hybrid_consistency()
TP_BOUND_VIOLATION <- check_tp_bound()

################################################################################
## T2 -- differential abundance results
##
## Sparse SSRV appears once, under the hybrid rule at the default bandwidth --
## the method as the manuscript proposes it. All nine settings are in the variant
## table Task 3 builds.
################################################################################
tab_da_results <- function() {
  keys  <- setdiff(names(DATASET_META), DROP_FROM_RESULTS)
  cols  <- c("sparse_SSRV", METHODS_OTHER)
  scen  <- vapply(DATASET_META[keys], `[[`, character(1), "scenario")

  rows <- character(0); seen <- character(0)
  for (k in keys) {
    m <- DATASET_META[[k]]; d <- dims_of(k)
    first <- !(m$scenario %in% seen)
    if (first && length(seen)) rows <- c(rows, "        \\midrule")
    sc <- if (first) scenario_cell(m$scenario, sum(scen == m$scenario), block = TRUE) else ""
    seen <- c(seen, m$scenario)
    if (!(m$short %in% rownames(fptp))) {
      stop("no results row for ", k, " (looked for '", m$short, "' in real_data_results.RData)")
    }
    ## The raw "FP/(TP)" strings first, so the best ratio in the row can be
    ## chosen before any of them is turned into LaTeX.
    raw <- vapply(cols, function(cc) {
      if (!(cc %in% colnames(fptp))) NA_character_ else fptp[m$short, cc]
    }, character(1))
    hit <- best_fptp(raw)
    cells <- vapply(seq_along(cols), function(j) {
      if (is.na(raw[j])) DASH else fptp_cell(raw[j], bold = hit[j])
    }, character(1))
    if (!task_quiet())
      message(sprintf("  %-12s best FP:TP -> %s", m$short,
                      paste(cols[hit], collapse = ", ")))
    ## Two lines, year underneath -- as in T1. A ONE-line \shortstack in a row
    ## whose other cells are two-line \begin{tabular}[t] blocks sits on the row
    ## BASELINE, i.e. the lower of the two lines, so the dataset name dropped
    ## away from the numbers it labels and left a gap above it. Two lines fill
    ## the row and the name lands level with the false-positive figure.
    rows <- c(rows, row_tex(sc, stack_paren(m$label, block = TRUE),
                            tex_block(c(sprintf("%d $\\times$ %d", d$D, d$N),
                                        sprintf("(%d)", n_truth_of(k))),
                                      sep = " \\\\ ", align = "c"),
                            cells))
  }

  spec <- paste0("c c c ", paste(rep("c", length(cols)), collapse = " "))
  hdr  <- paste0("        & & \\shortstack{D $\\times$ N \\\\ (Truth)} & ",
                 ## tex_escape: METHOD_LABEL is shared with a ggplot figure and so holds
                 ## plain text, but "ALDEx2_def" carries an underscore that LaTeX reads as
                 ## a subscript. Escaped here rather than in the vector.
                 paste(rot(tex_escape(unname(METHOD_LABEL[cols]))), collapse = " & "), " \\\\")
  arm_note <- paste0(" Sparse SSRV is reported under the hybrid rule at the default bandwidth, that ",
                     "is, the scale model the selection rule chose for that dataset. ",
                     "Table~\\ref{tab:ssrv_variants_tpfp} reports all nine combinations of scale-model ",
                     "rule and bandwidth, the two scale models separately among them.")
  c("\\begin{table}[!tbp]",
    "    \\centering",
    "    \\renewcommand{\\arraystretch}{1.2}",
    "    \\resizebox{\\textwidth}{!}{%",
    paste0("    \\begin{tabular}{", spec, "}"),
    "        \\toprule",
    hdr,
    "        \\midrule",
    rows,
    "        \\bottomrule",
    "    \\end{tabular}%",
    "    }",
    paste0("    \\caption{Differential abundance results on the real datasets. Each cell gives the number ",
           "of false positives over the number of true positives in parentheses, scored against the ",
           "study's truth set, whose size is given beside the dimensions. \\textbf{Bold} marks the best ",
           "false-to-true positive ratio in each row. Every row attains a ratio of zero, so ties on the ",
           "ratio are broken by the largest number of true positives recovered at that ratio; methods ",
           "recovering none of the truth set are not eligible, their ratio being undefined. Where several ",
           "methods remain tied they are all marked. \\citet{jin2022quantsinglecell} ",
           "is scored on its pooled cOTU run, the two cecal locations analysed together: its ",
           "measured-scale reference calls one cOTU pooled and none at either location alone, so only ",
           "the pooled run supports a truth set. Its ASV and OTU representations have no measured scale ",
           "and are reported separately in the Jin tables.", arm_note,
           " A dash marks a figure that was withheld rather than computed.}"),
    "    \\label{tab:real_studies_result}",
    "\\end{table}")
}

################################################################################
## EMIT
################################################################################
T2 <- tab_da_results()
write_table(T2, "T2_da_results_hybrid")

## The old T3_da_results_by_arm.tex is deliberately not written any more. Remove
## the stale file so a \\input of it cannot go on typesetting a table nothing
## rebuilds -- the failure mode that keeps a deleted table alive in a manuscript.
old_T3 <- file.path(DIR_TAB, "T3_da_results_by_arm.tex")
if (file.exists(old_T3)) {
  file.remove(old_T3)
  message("Removed superseded ", old_T3,
          " (its two arms are columns of S1_ssrv_variants_tpfp.tex now)")
}

if (!interactive() && !task_quiet()) cat(T2, sep = "\n")
