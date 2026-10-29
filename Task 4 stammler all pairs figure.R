################################################################################
## TASK 4 -- Figure: every method on all fifteen Stammler pool pairs.
##
##   Rscript "Task 4 stammler all pairs figure.R"
##
## \citet{stammler2016adjusting} spiked six pools at known, different levels, so
## every one of the 15 pool-vs-pool contrasts is a labelled problem with the same
## two true positives. One panel per method; inside a panel the upper triangle is
## the true positive rate of that comparison and the lower triangle its false
## discovery rate, so a method that finds the spike-ins without inventing
## anything is dark blue above the diagonal and white below it.
##
## Sparse SSRV appears ONCE here, under the hybrid rule at the default bandwidth
## -- the method as the manuscript proposes it -- so this figure compares methods
## rather than settings. All nine Sparse SSRV variants are in the Task 5 figure
## and in Table~\ref{tab:ssrv_variants_stammler_pairs}.
##
## SOURCES, per pair, resolved to the NEWEST datestamp that exists on disk:
##   res_tpfp_<ref>_<comp>_<stamp>.RData   every competing method
##   sparse_arm_sig.rds                    the per-arm record, used only to check
##                                         that the sparse_SSRV row above really
##                                         is the selected arm of the same fit
##   stammler_mrckf_out_all_pairs.RDS      CKF, for runs whose res_tpfp predates
##                                         CKF being folded in
##
## OUTPUT   manuscript_figures/stammler_all_pairs_methods.pdf
##
## Supersedes 2016_stammler_microbiome_micehuman/all_pairs_results.R, which named
## one datestamp in the file and so went on plotting the March fits after the
## August re-run, and which carried five Sparse SSRV panels among its sixteen.
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

## Total positives per comparison: the two spike-in taxa, every pair.
P_TRUTH <- length(STAMMLER_TRUTH)

## ALDEx2_original is the un-tuned call that the results tables do not report.
DROP_METHODS <- "ALDEx2_original"

################################################################################
## READ
################################################################################
PAIRS <- stammler_pair_analyses()

rates <- function(TP, FP) {
  data.frame(TPR = TP / P_TRUTH,
             FDR = ifelse(TP + FP == 0, 0, FP / (TP + FP)))
}

## Every competing method, from each pair's own results object.
method_results <- do.call(rbind, lapply(PAIRS, function(a) {
  if (is.na(a$tpfp)) {
    warning("no res_tpfp for ", a$label, " under any known datestamp -- that ",
            "column of every panel will be blank", call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  m <- as.data.frame(readRDS(a$tpfp)$aldex_gold)
  m$Method <- rownames(m)
  cbind(m[, c("TP", "FP", "Method")], Reference = a$ref, Comparison = a$comp,
        rates(m$TP, m$FP), row.names = NULL)
}))

## CKF, only for pairs whose res_tpfp does not already carry it. The August runs
## fold CKF in; the March ones do not, and adding it unconditionally from the
## saved MRCKF object would give those pairs two CKF rows -- which geom_tile
## paints on top of each other, so the duplicate would never be visible.
have_ckf <- with(method_results[method_results$Method == "CKF", ],
                 paste(Reference, Comparison))
need_ckf <- Filter(function(a) !(paste(a$ref, a$comp) %in% have_ckf), PAIRS)
if (length(need_ckf)) {
  f <- file.path(DIR_STAMMLER, "stammler_mrckf_out_all_pairs.RDS")
  if (file.exists(f)) {
    ck <- readRDS(f)
    ## Positional: the saved list is in combn() order, the same order PAIRS is
    ## built in. Checked rather than assumed, because a silently misaligned CKF
    ## panel would attribute one comparison's result to another.
    if (length(ck) != length(PAIRS)) {
      warning("stammler_mrckf_out_all_pairs.RDS holds ", length(ck),
              " entries but there are ", length(PAIRS), " pairs -- CKF is ",
              "withheld rather than matched up by position",
              call. = FALSE, immediate. = TRUE)
    } else {
      idx <- match(vapply(need_ckf, `[[`, character(1), "pair"),
                   vapply(PAIRS, `[[`, character(1), "pair"))
      add <- do.call(rbind, Map(function(a, i) {
        mn <- colMeans(ck[[i]])
        cbind(data.frame(TP = mn[1], FP = mn[2], Method = "CKF",
                         Reference = a$ref, Comparison = a$comp),
              rates(mn[1], mn[2]), row.names = NULL)
      }, need_ckf, idx))
      method_results <- rbind(method_results, add)
      message("CKF taken from stammler_mrckf_out_all_pairs.RDS for ",
              length(need_ckf), " pair(s) whose results file predates it")
    }
  }
}
method_results <- method_results[!(method_results$Method %in% DROP_METHODS), ]

################################################################################
## THE HYBRID SPARSE SSRV ROW MUST BE THE SELECTED ARM OF THE SAME FIT
##
## The sparse_SSRV row above comes from the pair's res_tpfp; the arm record comes
## from record_scale_model() partway through the same script. A DA run that was
## interrupted between the two leaves a fresh arm record beside a stale method
## table, both under one datestamp, and the panel then shows a fit nothing else
## does. Compared here because the filename cannot show it.
################################################################################
check_sparse_row <- function(d) {
  v <- ssrv_variants_stammler_pairs()
  v <- v[v$variant == "hybrid|nrd0", ]
  bad <- character(0)
  for (a in PAIRS) {
    r <- d[d$Method == "sparse_SSRV" & d$Reference == a$ref & d$Comparison == a$comp, ]
    q <- v[v$pair == a$pair, ]
    if (!nrow(r) || !nrow(q) || is.na(q$TP[1])) next
    if (r$TP[1] != q$TP[1] || r$FP[1] != q$FP[1]) {
      bad <- c(bad, sprintf("%s: results file %d/(%d) but selected arm %d/(%d)",
                            a$label, r$FP[1], r$TP[1], q$FP[1], q$TP[1]))
    }
  }
  if (length(bad)) {
    warning("the Sparse SSRV panel does not match the recorded selected arm on:\n  ",
            paste(bad, collapse = "\n  "),
            "\nThose pairs' method table and per-arm record came from different runs; ",
            "re-run the pairs loop of the Stammler DA script to completion.",
            call. = FALSE, immediate. = TRUE)
  } else {
    message("Sparse SSRV matches the recorded selected arm on every pool pair")
  }
  invisible(bad)
}
SPARSE_ROW_MISMATCH <- check_sparse_row(method_results)

################################################################################
## PLOT
##
## Panel order follows Table~\ref{tab:real_studies_result}: Sparse SSRV first,
## then the competitors in the order that table lists them. Twelve panels fill a
## four-column grid exactly.
################################################################################
PANEL_ORDER <- c("sparse_SSRV", METHODS_OTHER)
PANEL_LABEL <- c(sparse_SSRV = unname(METHOD_LABEL[["sparse_SSRV"]]),
                 METHOD_LABEL[METHODS_OTHER])

unknown <- setdiff(unique(method_results$Method), PANEL_ORDER)
if (length(unknown)) {
  warning("dropping method(s) with no panel: ", paste(unknown, collapse = ", "),
          " -- add them to METHOD_LABEL in task_common.R to plot them",
          call. = FALSE, immediate. = TRUE)
}

## The mirror layout: the same pair appears once above the diagonal carrying its
## TPR and once below carrying its FDR, so one square grid shows both.
mirror <- bind_rows(
  transform(method_results, X = Reference, Y = Comparison, Value = TPR, MetricType = "TPR"),
  transform(method_results, X = Comparison, Y = Reference, Value = FDR, MetricType = "FDR"))

mirror <- mirror %>%
  filter(Method %in% PANEL_ORDER) %>%
  mutate(X = factor(X, levels = STAMMLER_POOLS),
         Y = factor(Y, levels = STAMMLER_POOLS),
         Method = factor(unname(PANEL_LABEL[Method]),
                         levels = unname(PANEL_LABEL[PANEL_ORDER])))

## Blank rather than paint a tile with no value behind it -- geom_tile draws NA
## in its own grey, which is easy to read as a low-but-real number.
mirror <- mirror[!is.na(mirror$Value), ]

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
  facet_wrap(~ Method, ncol = 4, drop = FALSE) +
  theme_minimal() +
  labs(x = "Reference Pool", y = "Comparison Pool") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

f_out <- file.path(DIR_FIG, "stammler_all_pairs_methods.pdf")
ggsave(f_out, p, width = 12, height = 9)
message("Wrote ", f_out)
report_stamps()
