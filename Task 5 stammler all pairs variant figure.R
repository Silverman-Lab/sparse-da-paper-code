################################################################################
## TASK 5 -- Figure: all nine Sparse SSRV variants on all fifteen Stammler pairs.
##
##   Rscript "Task 5 stammler all pairs variant figure.R"
##
## The Task 4 figure compares methods and shows Sparse SSRV once, under the
## hybrid rule at the default bandwidth. This is the same fifteen comparisons
## with the method's own two settings opened up: three scale-model rules (hybrid,
## bootstrap, Laplace) down the rows, three KDE bandwidths (nrd0, nrd, SJ-ste)
## across the columns.
##
## Read it as a grid, not as nine unrelated panels:
##   * down a column   what the scale model costs at a fixed bandwidth
##   * across a row    what the bandwidth costs at a fixed scale model. The
##                     Laplace scale variance is O(h^2), so this row is where a
##                     bandwidth effect should appear if there is one.
##   * the hybrid row  always equals one of the two rows beneath it, cell by
##                     cell; where it changes which one, the selection rule is
##                     switching arms.
##
## Inside each panel, as in Task 4: upper triangle is the true positive rate of
## that pool pair, lower triangle its false discovery rate, both out of the same
## two spike-in taxa.
##
## SOURCE   ssrv_variants_stammler_pairs(), in task_common.R -- the same reader
##          that builds Table~\ref{tab:ssrv_variants_stammler_pairs}, so this
##          figure and that table cannot report different numbers.
##
## OUTPUT   manuscript_figures/stammler_all_pairs_ssrv_variants.pdf
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

suppressMessages({ library(ggplot2); library(dplyr); library(ggnewscale) })

P_TRUTH <- length(STAMMLER_TRUTH)

################################################################################
## READ
################################################################################
V <- ssrv_variants_stammler_pairs()

missing <- tapply(V$n_called, V$variant, function(z) all(is.na(z)))
if (any(missing)) {
  warning("no record at all for ", paste(names(missing)[missing], collapse = ", "),
          " -- those panels will be empty. Re-run the pairs loop of the Stammler ",
          "DA script so check_sparse_bw_sens() records every bandwidth.",
          call. = FALSE, immediate. = TRUE)
}

V <- V %>%
  filter(!is.na(TP)) %>%
  mutate(TPR = TP / P_TRUTH,
         FDR = ifelse(TP + FP == 0, 0, FP / (TP + FP)))

################################################################################
## HOW OFTEN THE SELECTION RULE ACTUALLY MOVES
##
## The whole point of showing the hybrid row beside the two arms is to say
## whether choosing between them matters here. Reported to the console rather
## than onto the figure, which has no room for it, but it belongs in the caption.
################################################################################
sel <- V[V$rule == "hybrid", c("pair", "bandwidth", "selected")]
tab <- table(sel$bandwidth, sel$selected)
cat("\n== which arm the rule selected, by bandwidth (15 pool pairs) ==\n")
print(tab)
flip <- tapply(sel$selected, sel$pair, function(z) length(unique(z[!is.na(z)])) > 1)
if (any(flip, na.rm = TRUE)) {
  message("The selection rule changes arms with the bandwidth on: ",
          paste(names(flip)[which(flip)], collapse = ", "))
} else {
  message("The selection rule picks the same arm at all three bandwidths on every pair")
}

################################################################################
## PLOT
##
## Same mirror layout as Task 4: each pair appears once above the diagonal
## carrying its TPR and once below carrying its FDR.
################################################################################
mirror <- bind_rows(
  transform(V, X = ref,  Y = comp, Value = TPR, MetricType = "TPR"),
  transform(V, X = comp, Y = ref,  Value = FDR, MetricType = "FDR"))

mirror <- mirror %>%
  mutate(X = factor(X, levels = STAMMLER_POOLS),
         Y = factor(Y, levels = STAMMLER_POOLS),
         rule = factor(unname(SSRV_RULE_LABEL[rule]),
                       levels = unname(SSRV_RULE_LABEL[SSRV_RULES])),
         bandwidth = factor(bandwidth, levels = SSRV_BWS)) %>%
  filter(!is.na(Value))

p <- ggplot(mirror) +
  geom_tile(data = filter(mirror, MetricType == "TPR"), aes(x = X, y = Y, fill = Value)) +
  scale_fill_gradient(low = "white", high = "#084594", name = "TPR",
                      limits = c(0, 1), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_tile(data = filter(mirror, MetricType == "FDR"), aes(x = X, y = Y, fill = Value)) +
  scale_fill_gradient(low = "white", high = "#cb181d", name = "FDR",
                      limits = c(0, 1), guide = guide_colorbar(order = 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  ## Rows are the scale model, columns the bandwidth. facet_grid rather than
  ## facet_wrap so the two settings stay on their own axes and a row or a column
  ## can be read as one comparison.
  facet_grid(rule ~ bandwidth, drop = FALSE) +
  theme_minimal() +
  labs(x = "Reference Pool", y = "Comparison Pool",
       title = NULL,
       caption = paste0("Rows: scale model. Columns: KDE bandwidth. ",
                        "The hybrid row duplicates, cell by cell, whichever arm ",
                        "the selection rule chose at that bandwidth.")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, colour = "grey30"))

f_out <- file.path(DIR_FIG, "stammler_all_pairs_ssrv_variants.pdf")
ggsave(f_out, p, width = 10, height = 9)
message("Wrote ", f_out)
report_stamps()
