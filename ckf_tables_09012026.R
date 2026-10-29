################################################################################
## ckf_tables_09012026.R -- the two FDR / Power tables for the sign-balance
## study, built from the saved output of ckf_simulation_09012026.R.
##
##   manuscript_tables/T3_ckf_sign_balance.tex        CKF (both targets) + sparse SSRV
##   manuscript_tables/S4_ckf_sign_balance_arms.tex   the same, plus each scale model
##
## Reads only the saved RDS, so the tables can be rebuilt without refitting.
################################################################################

suppressMessages(library(dplyr))

#ROOT <- "/storage/work/wpg5129/sparse_project"   # HPC
ROOT <- local({  # manuscript_code: the folder holding .manuscript_code_root (README.md)
  p <- normalizePath(getwd(), mustWork = FALSE)
  while (!file.exists(file.path(p, ".manuscript_code_root"))) {
    if (dirname(p) == p) stop("run from inside the manuscript_code folder (see README.md)")
    p <- dirname(p)
  }
  p
})
ROOT <- path.expand(ROOT)
RES  <- file.path(ROOT, "variable selection on CLR", "results with sparsity 0.20")
TAB  <- file.path(ROOT, "manuscript_tables")

## Same pair the simulation script takes. SFX keeps a non-default pi from
## overwriting the pi = 0.2 .tex files and their labels.
PI   <- as.numeric(Sys.getenv("CKF_PI",  "0.20"))
TAG  <- Sys.getenv("CKF_TAG", "09012026")
SFX  <- if (TAG == "09012026") "" else paste0("_", sub("_09012026$", "", TAG))
DTAX <- 400

df_sparse <- readRDS(file.path(RES, paste0("df_sparse_", TAG, ".RData")))
df_sparse$fdr   <- as.numeric(df_sparse$fdr)
df_sparse$power <- as.numeric(df_sparse$power)
df_arm    <- readRDS(file.path(RES, paste0("sparse_arm_record_", TAG, ".RData")))

## from the data, so a run at a different replicate count cannot ship a
## caption that misstates it
NREP <- nrow(df_sparse) / 3L

## tpr = TP / num.spiked, so this recovers the truth-set size the fits actually
## used. Checked against PI because these tables have no other way to notice
## that CKF_PI and CKF_TAG disagree about which run the RDS came from.
NSPIKE <- unique(round(df_arm$TP_boot / df_arm$tpr_boot, 6))
NSPIKE <- NSPIKE[is.finite(NSPIKE)]
stopifnot(length(NSPIKE) == 1, isTRUE(all.equal(NSPIKE, PI * DTAX)))

df_sparse_summary <- df_sparse %>%
  group_by(effect_dist) %>%
  summarize(fdr_mean = mean(fdr, na.rm = TRUE),
            power_mean = mean(power, na.rm = TRUE))

df_arm_summary <- df_arm %>%
  group_by(effect_dist) %>%
  summarize(prop_bootstrap = mean(scale_model == "bootstrap", na.rm = TRUE),
            prop_laplace   = mean(scale_model == "laplace",   na.rm = TRUE),
            fdr_bootstrap  = mean(fdr_boot, na.rm = TRUE),
            tpr_bootstrap  = mean(tpr_boot, na.rm = TRUE),
            fdr_laplace    = mean(fdr_lap,  na.rm = TRUE),
            tpr_laplace    = mean(tpr_lap,  na.rm = TRUE))

avg_ckf <- function(f){
  out.list <- readRDS(file.path(RES, paste0(f, "_", TAG, ".RData")))
  out.list <- lapply(out.list, function(mat){ mat[is.nan(mat)] <- 0; mat })
  Reduce(`+`, out.list) / length(out.list)
}
df_balanced    <- avg_ckf("out_balanced")
df_1dec4inc    <- avg_ckf("out_1dec4inc")
df_allincrease <- avg_ckf("out_allincrease")

scen <- c("Symmetric", "80% Positive", "All Positive")
meth <- c("CKF target FDR 0.1", "CKF target FDR 0.05", "Sparse SSRV",
          "Sparse SSRV (bootstrap arm)", "Sparse SSRV (laplace arm)")
lev  <- c("Symmetric (Sum = 0)", "80% Positive", "All Positive")
sp.i <- match(lev, df_sparse_summary$effect_dist)
ar.i <- match(lev, df_arm_summary$effect_dist)

fdr_df <- matrix(NA, nrow = 5, ncol = 3, dimnames = list(meth, scen))
fdr_df[1:2,1] <- df_balanced[,"FDP_S"]
fdr_df[1:2,2] <- df_1dec4inc[,"FDP_S"]
fdr_df[1:2,3] <- df_allincrease[,"FDP_S"]
fdr_df[3,] <- df_sparse_summary$fdr_mean[sp.i]
fdr_df[4,] <- df_arm_summary$fdr_bootstrap[ar.i]
fdr_df[5,] <- df_arm_summary$fdr_laplace[ar.i]

power_df <- matrix(NA, nrow = 5, ncol = 3, dimnames = list(meth, scen))
power_df[1:2,1] <- df_balanced[,"power"]
power_df[1:2,2] <- df_1dec4inc[,"power"]
power_df[1:2,3] <- df_allincrease[,"power"]
power_df[3,] <- df_sparse_summary$power_mean[sp.i]
power_df[4,] <- df_arm_summary$tpr_bootstrap[ar.i]
power_df[5,] <- df_arm_summary$tpr_laplace[ar.i]

cell <- matrix(paste0(formatC(fdr_df, digits = 2, format = "f"), " / ",
                      formatC(power_df, digits = 2, format = "f")),
               nrow = 5, dimnames = list(meth, scen))

## LaTeX -----------------------------------------------------------------------
TEX_ROW <- c("CKF target FDR 0.1", "CKF target FDR 0.05", "Sparse SSRV",
             "\\quad bootstrap scale model", "\\quad Laplace scale model")

ckf_table <- function(rows, caption, label, arraystretch = "1.2") {
  body <- vapply(rows, function(i)
    paste0("        ", TEX_ROW[i], " & ", paste(cell[i, ], collapse = " & "), " \\\\"),
    character(1))
  c("\\begin{table}[!tbp]",
    "    \\centering",
    paste0("    \\renewcommand{\\arraystretch}{", arraystretch, "}"),
    "    \\begin{tabular}{l c c c}",
    "        \\toprule",
    "        & Symmetric & 80\\% Positive & All Positive \\\\",
    "        \\midrule",
    body,
    "        \\bottomrule",
    "    \\end{tabular}",
    paste0("    \\caption{", caption, "}"),
    paste0("    \\label{", label, "}"),
    "\\end{table}")
}

SETTING <- paste0(
  "Each cell reports the mean false discovery rate and the mean power over ",
  NREP, " simulated datasets, separated by a slash (FDR / Power). Data are ",
  "generated by \\texttt{SparseDOSSA2} with $N = 100$, $D = 400$, ",
  "$\\pi = ", PI, "$ spiked taxa (", NSPIKE, " of ", DTAX,
  ") and read depth $\\lambda = 10^{4} D$. The three ",
  "columns differ only in the signs of the spiked log fold changes: symmetric ",
  "(equal numbers up and down, so the true effects sum to zero), 80\\% ",
  "positive, and all positive.")

t1 <- ckf_table(
  1:3,
  paste0("Empirical FDR and power as the signs of the effect sizes become ",
         "unbalanced. ", SETTING, " CKF is the compositional knockoff filter at ",
         "each of its two target FDR levels; sparse SSRV uses the hybrid ",
         "scale-model rule at a nominal level of 0.05. CKF's FDR rises with the ",
         "imbalance and exceeds its own target in the all-positive column, ",
         "while sparse SSRV is unaffected."),
  paste0("tab:ckf_sign_balance", SFX))

sel <- formatC(as.numeric(df_arm_summary$prop_bootstrap[ar.i]) * 100,
               digits = 1, format = "f")
t2 <- ckf_table(
  1:5,
  paste0("Empirical FDR and power as the signs of the effect sizes become ",
         "unbalanced, with sparse SSRV decomposed into its two scale models. ",
         SETTING, " The two indented rows are the bootstrap and Laplace scale ",
         "models scored on their own, on the same fits; the sparse SSRV row is ",
         "whichever one the width rule selected in each replicate. That rule ",
         "chose the bootstrap arm in ", sel[1], "\\%, ", sel[2], "\\% and ",
         sel[3], "\\% of replicates across the three columns, which is why the ",
         "sparse SSRV row nearly duplicates it. Taken alone, the Laplace scale ",
         "model exceeds the nominal 0.05 level in every column."),
  paste0("tab:ckf_sign_balance_arms", SFX))

dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
writeLines(t1, file.path(TAB, paste0("T3_ckf_sign_balance", SFX, ".tex")))
writeLines(t2, file.path(TAB, paste0("S4_ckf_sign_balance_arms", SFX, ".tex")))

print(as.data.frame(cell))
cat("\nwrote", file.path(TAB, paste0("T3_ckf_sign_balance", SFX, ".tex")), "\n")
cat("wrote", file.path(TAB, paste0("S4_ckf_sign_balance_arms", SFX, ".tex")), "\n")
