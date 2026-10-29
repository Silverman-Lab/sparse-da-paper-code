################################################################################
## shared.R -- load everything simulation_rerun needs, sourcing the UNCHANGED
## files read-only from ../simulation so there is exactly one copy to maintain.
##
## Nothing here writes to ../simulation.
##
## Only two files are copied rather than sourced, because they had to change:
##   src/sparse_SSRV_rerun.R    edits (a) and (b)
##   src/run_methods_rerun.R    edit (c)   [Study A / B only]
##
## shared_setup(need) controls how much is loaded. The pilot does not touch a
## single DA method, and attaching ALDEx2 / DESeq2 / ANCOMBC / baySeq costs
## several seconds PER WORKER on a PSOCK cluster -- with 7 workers over ~1850
## replicates that is not free, and it also pulls in packages that mask base
## generics for no reason.
################################################################################

shared_setup <- function(need = c("pilot", "methods")) {
  need <- match.arg(need)

  if (!exists("HERE", inherits = TRUE))
    stop("source config.R before shared.R")

  suppressWarnings(suppressMessages({
    library(SparseDOSSA2)
    # rdirichlet: dirmult's, NOT MCMCpack's. Study A's worker
    # (sparse_analysis_trycatch.R:11) attaches dirmult and never MCMCpack, so
    # that is the one the published results used. The two draw differently, so
    # mixing them would silently make the pilot non-comparable to Study A.
    # (simulation/scale_choice/run_scale_choice.R uses MCMCpack -- which is one
    # reason its numbers are not directly comparable to the main sweep.)
    library(dirmult)
  }))

  # sim_sparse_auto() and its KDE helpers. Self-contained: no library(), no
  # source(), no hardcoded path of its own.
  source(file.path(HERE, "src", "sparse_SSRV_rerun.R"))
  source(file.path(HERE, "src", "simulate_one.R"))

  if (need == "methods") {
    suppressWarnings(suppressMessages({
      library(phyloseq); library(ALDEx2); library(DESeq2); library(modeest)
      library(ANCOMBC);  library(S4Vectors); library(LinDA); library(limma)
      library(edgeR);    library(baySeq); library(numDeriv)
      library(MRCKF);    library(knockoff); library(compositions)
      # stringr is used by the ANCOMBC2 retry handler. The original worker never
      # attached it and only worked because something else pulled it in.
      library(stringr)
    }))

    source(file.path(OLD_SIM, "src", "01a_main_functions.R"))
    source(file.path(OLD_SIM, "src", "01b_helper_functions.R"))
    source(file.path(OLD_SIM, "src", "02_simulations.R"))
    source(file.path(OLD_SIM, "src", "common", "mrckf_fit.R"))
    source(file.path(OLD_SIM, "src", "common", "tidy_metrics.R"))
    source(file.path(OLD_SIM, "src", "common", "io_compat.R"))
    source(file.path(HERE,    "src", "run_methods_rerun.R"))

    # The three SSRV arms are three METHODS downstream, so they each need a
    # log2 multiplier. theta is already log2, hence 1. Appended here rather than
    # edited into tidy_metrics.R so ../simulation stays untouched.
    METHOD_SCALE <<- c(METHOD_SCALE,
                       "Sparse_SSRV_L" = 1,
                       "Sparse_SSRV_B" = 1)
  }

  invisible(TRUE)
}
