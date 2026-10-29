################################################################################
## TASK 3 -- Supplementary tables: ALL NINE Sparse SSRV variants.
##
##   Rscript "Task 3 sparse SSRV variant tables.R"
##
## Table 2 reports Sparse SSRV once, under the hybrid rule at the default
## bandwidth. This is the full grid behind that one column: three scale-model
## rules (hybrid, bootstrap, Laplace) crossed with three KDE bandwidths (nrd0,
## nrd, SJ-ste), on every analysis in the study -- the five scored datasets, the
## fifteen Stammler pool pairs, and Jin's six subsets.
##
## Why the grid and not two separate sensitivities: the Laplace arm's scale
## variance is O(h^2) in the bandwidth, so the two knobs are not independent.
## Reading a bandwidth sweep at the hybrid rule alone hides that the rule may be
## switching arms as h changes, and reading the arms at one bandwidth alone hides
## how much of an arm's behaviour is the bandwidth. The hybrid columns therefore
## carry a superscript naming the arm the rule picked AT THAT BANDWIDTH; a hybrid
## cell always duplicates the arm cell its superscript names.
##
## OUTPUTS
##   manuscript_tables/S1_ssrv_variants_tpfp.tex   the scored datasets, Jin's
##                                                 pooled run among them
##   manuscript_tables/S3_ssrv_variants_jin.tex    Jin's SUBSETS -- the two cecal
##                                                 locations separately and the
##                                                 ASV / OTU representations --
##                                                 as discovery counts
##
## The fifteen-pool-pair table is NOT emitted. It carried the same fifteen
## comparisons at the same nine settings as the Task 5 figure, which shows them
## as true positive rate and false discovery rate and is far easier to read; two
## presentations of one result is one too many.
##   real_data_ssrv_variants.RData / .csv                   the tidy frame behind
##                                                          all three, and behind
##                                                          the Task 5 figure
##
## JIN'S POOLED RUN IS SCORED; ITS SUBSETS ARE NOT. Pooling the two cecal
## locations is what gives the measured-scale reference enough power to call a
## feature -- one cOTU pooled, none at either location alone -- so the pooled run
## joins the scored table while the stratified locations and the scale-free ASV
## and OTU representations stay in S3 as discovery counts, which is the most they
## support.
##
## Nothing is refitted. The default bandwidth comes from each study's
## sparse_arm_sig.rds and the other two from its saved sensitivity fit; see the
## variant block of task_common.R for what is read from where.
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
## BUILD
################################################################################
message("\nReading the nine variants of every analysis ...")
V_SCORED   <- ssrv_variants_scored()
V_STAMMLER <- ssrv_variants_stammler_pairs()
V_JIN      <- ssrv_variants_jin()

## A variant that is missing everywhere means a study was never re-run under the
## sensitivity sweep, and the column would be an empty stripe in every table.
missing_everywhere <- function(d, what) {
  m <- tapply(d$n_called, d$variant, function(z) all(is.na(z)))
  if (any(m)) {
    warning(what, ": no record at all for ", paste(names(m)[m], collapse = ", "),
            " -- those columns are dashed. Re-run the DA scripts so ",
            "check_sparse_bw_sens() records every bandwidth.",
            call. = FALSE, immediate. = TRUE)
  }
  invisible(names(m)[m])
}
missing_everywhere(V_SCORED,   "scored datasets")
missing_everywhere(V_STAMMLER, "Stammler pool pairs")
missing_everywhere(V_JIN,      "Jin subsets")

################################################################################
## CROSS-CHECK AGAINST TABLE 2
##
## Five of the nine variants already appear in real_data_results.RData, put there
## by `real data results.R` from the same records: sparse_SSRV, sparse_SSRV_B and
## sparse_SSRV_L at the default bandwidth, and the hybrid rule at the two other
## bandwidths. They must agree cell for cell. When they do not, one of the two
## files is reading a run the other is not -- which is the single failure mode
## these tables cannot survive, since the supplementary grid exists precisely to
## explain the main table's column.
################################################################################
check_against_main <- function(d) {
  map <- c(`hybrid|nrd0`      = "sparse_SSRV",
           `bootstrap|nrd0`   = "sparse_SSRV_B",
           `laplace|nrd0`     = "sparse_SSRV_L",
           `hybrid|nrd`       = "sparse_SSRV_nrd",
           `hybrid|SJ-ste`    = "sparse_SSRV_SJ-ste")
  bad <- character(0); skipped <- character(0)
  for (k in setdiff(names(DATASET_META), DROP_FROM_RESULTS)) {
    short <- DATASET_META[[k]]$short
    if (!(short %in% rownames(fptp))) next
    for (v in names(map)) {
      cc <- map[[v]]
      if (!(cc %in% colnames(fptp))) { skipped <- union(skipped, cc); next }
      r <- d[d$dataset == k & d$variant == v, ]
      if (!nrow(r) || is.na(r$TP[1])) next
      want <- split_fptp(fptp[short, cc])
      if (is.na(want[1])) next                      # withheld in the main table
      got <- c(as.character(r$FP[1]), as.character(r$TP[1]))
      if (!identical(want, got)) {
        bad <- c(bad, sprintf("%s / %s: this table %s/(%s), real_data_results %s/(%s)",
                              k, v, got[1], got[2], want[1], want[2]))
      }
    }
  }
  if (length(skipped)) {
    message("Not cross-checked -- absent from real_data_results.RData: ",
            paste(sort(skipped), collapse = ", "))
  }
  if (length(bad)) {
    warning("the variant grid disagrees with the main results table:\n  ",
            paste(bad, collapse = "\n  "),
            "\nRe-run `real data results.R`, which reads the same records, so both ",
            "tables describe one fit.", call. = FALSE, immediate. = TRUE)
  } else {
    message("Every variant that Table 2 also reports agrees with it cell for cell")
  }
  invisible(bad)
}
MAIN_DISAGREEMENT <- check_against_main(V_SCORED)

################################################################################
## RENDERING
##
## One row per analysis, nine cells across, grouped three-and-three by rule.
################################################################################

## The nine columns in reading order: rule-major, so each rule's three
## bandwidths sit together under one spanning header.
VARIANT_COLS <- as.vector(t(outer(SSRV_RULES, SSRV_BWS, paste, sep = "|")))

## Column group headers, and the cmidrules that separate them. `lead` is the
## number of non-variant columns to the left.
variant_header <- function(lead, lead_labels) {
  starts <- lead + seq(1, by = length(SSRV_BWS), length.out = length(SSRV_RULES))
  span <- paste0(paste(rep("", lead), collapse = " & "),
                 if (lead) " & " else "",
                 paste(sprintf("\\multicolumn{%d}{c}{%s}", length(SSRV_BWS),
                               unname(SSRV_RULE_LABEL[SSRV_RULES])),
                       collapse = " & "), " \\\\")
  rules <- paste(mapply(function(a, i)
      sprintf("\\cmidrule(%s){%d-%d}", if (i == length(starts)) "l" else "lr",
              a, a + length(SSRV_BWS) - 1),
      starts, seq_along(starts)), collapse = "")
  bws <- paste0("        ",
                paste(c(lead_labels, rep(unname(SSRV_BW_LABEL[SSRV_BWS]),
                                         length(SSRV_RULES))), collapse = " & "),
                " \\\\")
  c(paste0("        ", span), paste0("        ", rules), bws)
}

## One cell. `metric` picks what the number is: "tpfp" for a scored analysis,
## "n" for a discovery count. Hybrid cells carry the arm the rule chose at that
## bandwidth, because a hybrid figure with no arm attached cannot be read.
variant_cell <- function(r, metric) {
  if (!nrow(r)) return(DASH)
  mark <- ""
  if (identical(r$rule[1], "hybrid") && !is.na(r$selected[1])) {
    mark <- paste0("\\textsuperscript{", toupper(substr(r$selected[1], 1, 1)), "}")
  }
  if (identical(metric, "n")) {
    if (is.na(r$n_called[1])) return(DASH)
    return(paste0(r$n_called[1], mark))
  }
  if (is.na(r$TP[1]) || is.na(r$FP[1])) return(DASH)
  tex_block(paste0(r$FP[1], mark, " \\\\ (", r$TP[1], ")"), align = "c")
}

## The nine cells of one analysis.
variant_row_cells <- function(d, dataset, pair, metric) {
  vapply(VARIANT_COLS, function(v) {
    r <- d[d$dataset == dataset & d$variant == v &
           ((is.na(pair) & is.na(d$pair)) |
            (!is.na(pair) & !is.na(d$pair) & d$pair == pair)), ]
    variant_cell(r, metric)
  }, character(1))
}

## Common table scaffolding, so the three tables differ only in their rows.
variant_table <- function(rows, lead_spec, lead_labels, caption, label,
                          arraystretch = "1.2") {
  spec <- paste0(lead_spec, " ", paste(rep("c", length(VARIANT_COLS)), collapse = " "))
  c("\\begin{table}[!tbp]",
    "    \\centering",
    paste0("    \\renewcommand{\\arraystretch}{", arraystretch, "}"),
    "    \\resizebox{\\textwidth}{!}{%",
    paste0("    \\begin{tabular}{", spec, "}"),
    "        \\toprule",
    variant_header(length(lead_labels), lead_labels),
    "        \\midrule",
    rows,
    "        \\bottomrule",
    "    \\end{tabular}%",
    "    }",
    paste0("    \\caption{", caption, "}"),
    paste0("    \\label{", label, "}"),
    "\\end{table}")
}

## The sentence that explains the superscripts, shared by all three captions.
HYBRID_NOTE <- paste0(
  "The hybrid rule selects one of the two scale models for each fit, so every hybrid entry ",
  "duplicates the entry of the arm its superscript names -- \\textsuperscript{B} bootstrap, ",
  "\\textsuperscript{L} Laplace -- in the same bandwidth column. Where the superscript changes ",
  "down a column, the selection rule is switching arms as the bandwidth changes.")

################################################################################
## S1 -- the scored datasets
################################################################################
tab_variants_scored <- function(d) {
  keys <- setdiff(names(DATASET_META), DROP_FROM_RESULTS)
  scen <- vapply(DATASET_META[keys], `[[`, character(1), "scenario")
  rows <- character(0); seen <- character(0)
  for (k in keys) {
    m <- DATASET_META[[k]]; dm <- dims_of(k)
    first <- !(m$scenario %in% seen)
    if (first && length(seen)) rows <- c(rows, "        \\midrule")
    sc <- if (first) scenario_cell(m$scenario, sum(scen == m$scenario)) else ""
    seen <- c(seen, m$scenario)
    rows <- c(rows, row_tex(sc, paste0("\\shortstack{", m$label, "}"),
                            tex_block(c(sprintf("%d $\\times$ %d", dm$D, dm$N),
                                        sprintf("(%d)", n_truth_of(k))),
                                      sep = " \\\\ ", align = "c"),
                            variant_row_cells(d, k, m$pair, "tpfp")))
  }
  variant_table(
    rows,
    lead_spec   = "c c c",
    lead_labels = c("Scenario", "Dataset", "\\shortstack{D $\\times$ N \\\\ (Truth)}"),
    caption = paste0(
      "Sparse SSRV under all nine combinations of scale-model rule and kernel bandwidth, on the ",
      "datasets of Table~\\ref{tab:real_studies_result}. Each cell gives the number of false ",
      "positives over the number of true positives in parentheses, scored against the same truth ",
      "set as Table~\\ref{tab:real_studies_result}, whose size is given beside the dimensions. The ",
      "column reported in the main text is Hybrid / nrd0. ", HYBRID_NOTE, " The Laplace scale ",
      "variance is $O(h^{2})$ in the bandwidth $h$, so the bandwidth is not a nuisance setting for ",
      "that arm: it is the scale uncertainty. \\citet{jin2022quantsinglecell} appears on its pooled ",
      "cOTU run, the only one of its subsets whose measured-scale reference calls anything; its two ",
      "cecal locations separately and its ASV and OTU representations are in ",
      "Table~\\ref{tab:ssrv_variants_jin}. A dash marks a fit that was not recorded."),
    label = "tab:ssrv_variants_tpfp")
}

################################################################################
## S2 -- the fifteen Stammler pool pairs
##
## The same fifteen comparisons Figure~\ref{fig:stammler_all_pairs} shows, as
## numbers. Every pair is scored against the same two spike-in taxa, so the truth
## set is stated once in the caption rather than repeated down a column.
################################################################################
tab_variants_stammler <- function(d) {
  rows <- character(0)
  for (a in stammler_pair_analyses()) {
    dm <- dims_for("Stammler2016", a$pair)
    rows <- c(rows, row_tex(paste0(a$ref, " vs ", a$comp),
                            sprintf("%d $\\times$ %d", dm$D, dm$N),
                            variant_row_cells(d, "Stammler2016", a$pair, "tpfp")))
  }
  variant_table(
    rows,
    lead_spec   = "l c",
    lead_labels = c("Pool pair", "D $\\times$ N"),
    caption = paste0(
      "Sparse SSRV under all nine combinations of scale-model rule and kernel bandwidth, on all ",
      "fifteen pairs of \\citet{stammler2016adjusting} spike-in pools. Each cell gives the number ",
      "of false positives over the number of true positives in parentheses. Every pair is scored ",
      "against the same two spike-in taxa, so the truth set has size 2 throughout, and a pair whose ",
      "pools were prepared at similar spike-in levels offers little signal to find. ", HYBRID_NOTE,
      " Figure~\\ref{fig:stammler_all_pairs_variants} shows the same fifteen comparisons as true ",
      "positive rate and false discovery rate. A dash marks a fit that was not recorded."),
    label = "tab:ssrv_variants_stammler_pairs")
}

################################################################################
## S3 -- Jin's six subsets, as discovery counts
################################################################################
tab_variants_jin <- function(d) {
  rows <- character(0); seen <- character(0)
  subs <- jin_subsets()
  lev  <- vapply(subs, `[[`, character(1), "level")
  for (a in subs) {
    dm <- dims_for("Jin2022", a$pair)
    first <- !(a$level %in% seen)
    if (first && length(seen)) rows <- c(rows, "        \\midrule")
    lv <- if (first) sprintf("\\multirow{%d}{*}{%s}", sum(lev == a$level), a$level) else ""
    seen <- c(seen, a$level)
    rows <- c(rows, row_tex(lv,
                            if (a$location == "prox") "proximal" else "distal",
                            sprintf("%d $\\times$ %d", dm$D, dm$N),
                            variant_row_cells(d, "Jin2022", a$pair, "n")))
  }
  variant_table(
    rows,
    lead_spec   = "c c c",
    lead_labels = c("Representation", "Cecal location", "D $\\times$ N"),
    caption = paste0(
      "Sparse SSRV under all nine combinations of scale-model rule and kernel bandwidth, on the six ",
      "subsets of \\citet{jin2022quantsinglecell}. Entries are the NUMBER OF FEATURES CALLED, not ",
      "true and false positives: the measured-scale reference for that study calls no cOTU at an ",
      "adjusted significance level at either cecal location, so scoring against it would report the ",
      "power of the reference rather than of the methods. ", HYBRID_NOTE, " Only the cOTU rows ",
      "carry a measured scale; the ASV and OTU tables are conventional 16S gene-copy measurements ",
      "and are shown because the three representations of the same samples do not agree. A dash ",
      "marks a fit that was not recorded."),
    label     = "tab:ssrv_variants_jin",
    arraystretch = "1.1")
}

################################################################################
## EMIT
################################################################################
S1 <- tab_variants_scored(V_SCORED)
S3 <- tab_variants_jin(V_JIN)

write_table(S1, "S1_ssrv_variants_tpfp")
write_table(S3, "S3_ssrv_variants_jin")

## The pool-pair TABLE is retired in favour of the Task 5 figure, but V_STAMMLER
## is still built above: it goes into the tidy frame below, which is what Task 5
## reads, and it is what the console grid prints. Only the LaTeX is dropped.
old_S2 <- file.path(DIR_TAB, "S2_ssrv_variants_stammler_pairs.tex")
if (file.exists(old_S2)) {
  file.remove(old_S2)
  message("Removed superseded ", old_S2,
          " (the same fifteen pairs are Figure stammler_all_pairs_ssrv_variants.pdf)")
}

## The tidy frame behind all three, so the Task 5 figure reads the numbers this
## table printed rather than deriving its own, and so the grid can be inspected
## outside LaTeX.
VARIANTS_ALL <- rbind(
  cbind(V_SCORED,   group = "scored",         level = NA_character_, location = NA_character_,
        ref = NA_character_, comp = NA_character_),
  cbind(V_STAMMLER, group = "stammler_pairs", level = NA_character_, location = NA_character_),
  cbind(V_JIN,      group = "jin",            ref = NA_character_, comp = NA_character_))
VARIANTS_ALL <- VARIANTS_ALL[, c("group", "dataset", "pair", "label", "level", "location",
                                 "ref", "comp", "variant", "rule", "bandwidth", "selected",
                                 "n_called", "TP", "FP", "scored")]
saveRDS(VARIANTS_ALL, file.path(ROOT, "real_data_ssrv_variants.RData"))
write.csv(VARIANTS_ALL, file.path(ROOT, "real_data_ssrv_variants.csv"), row.names = FALSE)
message("Wrote ", file.path(ROOT, "real_data_ssrv_variants.RData"), " (",
        nrow(VARIANTS_ALL), " rows)")

if (!task_quiet()) report_stamps()

## A compact console view of the same grid, for reading without compiling LaTeX.
show_grid <- function(d, value = c("tpfp", "n")) {
  value <- match.arg(value)
  cell <- if (value == "tpfp") sprintf("%s/(%s)", d$FP, d$TP) else as.character(d$n_called)
  cell[is.na(d$TP) & value == "tpfp"] <- "--"
  cell[is.na(d$n_called)] <- "--"
  w <- reshape(data.frame(row = ifelse(is.na(d$pair), d$dataset,
                                       paste0(d$dataset, "/", d$pair)),
                          variant = factor(d$variant, levels = VARIANT_COLS),
                          cell = cell, stringsAsFactors = FALSE),
               direction = "wide", idvar = "row", timevar = "variant")
  names(w) <- sub("^cell\\.", "", names(w))
  w
}
if (!task_quiet()) {
  cat("\n== scored datasets, FP/(TP) ==\n");     print(show_grid(V_SCORED),   row.names = FALSE)
  cat("\n== Stammler pool pairs, FP/(TP) ==\n"); print(show_grid(V_STAMMLER), row.names = FALSE)
  cat("\n== Jin subsets, features called ==\n"); print(show_grid(V_JIN, "n"), row.names = FALSE)
}

if (!interactive() && !task_quiet())
  cat("\n", paste(c(S1, "", S3), collapse = "\n"), "\n", sep = "")
