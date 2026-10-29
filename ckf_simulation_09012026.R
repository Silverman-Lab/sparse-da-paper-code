################################################################################
## CKF vs sparse SSRV as the signs of the effect sizes become unbalanced.
## 09012026 rerun of ckf_simulation.R: sparse_auto() is now taken from
## sparse_SSRV.R -- the version the 09012026 real-data scripts use -- instead of
## the copy that used to be defined at the top of this file. Everything not
## needed for the FDR / Power table has been dropped.
##
## Recorded per replicate, alongside the selected (hybrid) result: which scale
## model the width rule picked, and the TPR / FDR of EACH arm on its own.
################################################################################

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

## pi and the tag it writes under. Defaults reproduce the 09012026 run; set both
## together to sweep pi. The guard is deliberate -- a sweep that reused the tag
## would silently overwrite another pi's results, leaving a directory whose
## files no longer come from one design.
PI  <- as.numeric(Sys.getenv("CKF_PI",  "0.20"))
TAG <- Sys.getenv("CKF_TAG", "09012026")
if (!isTRUE(all.equal(PI, 0.20)) && TAG == "09012026")
  stop("set CKF_TAG as well: pi = ", PI, " must not be written under the pi = 0.2 tag")

library(SparseDOSSA2)
library(dirmult)
library(dplyr)
library(MRCKF)
library(doParallel)
library(foreach)

source(file.path(ROOT, "sparse_SSRV.R"))

## fdr is NaN when an arm selects nothing, as in the original script
score_arm <- function(tab, truth, num.spiked){
  if(is.null(tab) || all(is.na(tab$sig.q)))
    return(c(TP = NA_real_, FP = NA_real_, tpr = NA_real_, fdr = NA_real_))
  sig <- rownames(tab)[which(tab$sig.q)]
  TP <- sum(sig %in% truth)
  FP <- sum(!sig %in% truth)
  c(TP = TP, FP = FP, tpr = TP/num.spiked, fdr = FP/(TP+FP))
}

ARM_NA <- data.frame(scale_model = NA_character_, n_laplace_fail = NA_integer_,
                     TP_boot = NA_real_, FP_boot = NA_real_, tpr_boot = NA_real_, fdr_boot = NA_real_,
                     TP_lap  = NA_real_, FP_lap  = NA_real_, tpr_lap  = NA_real_, fdr_lap  = NA_real_,
                     stringsAsFactors = FALSE)

sparse_out_alltypes <- list()
sparse_out_alltypes[[1]] <- matrix(NA, nrow = 1000, ncol = 4)
sparse_out_alltypes[[2]] <- matrix(NA, nrow = 1000, ncol = 4)
sparse_out_alltypes[[3]] <- matrix(NA, nrow = 1000, ncol = 4)
arm_rec_alltypes <- list()
out.list_j <- list()

names(sparse_out_alltypes) <- c("Symmetric (Sum = 0)", "All Positive", "80% Positive")

cl <- makeCluster(5)
registerDoParallel(cl)
## sourced once per worker rather than once per replicate
clusterExport(cl, c("ROOT", "score_arm", "ARM_NA"))
clusterEvalQ(cl, source(file.path(ROOT, "sparse_SSRV.R")))

for(j in 1:3){
  n <- 100
  d <- 400
  lambda <- d*10000
  pi <- PI

  S <- 2000
  pval <- 0.05

  out.list <- foreach(m = 1:1000, .packages = c("SparseDOSSA2", "MASS", "MRCKF", "dirmult")) %dopar% {

    seed.num <- m
    set.seed(seed.num)
    Stool <-
      Stool_mod <-
      SparseDOSSA2:::Stool
    feature_param_template <- cbind("pi0" = Stool$EM_fit$fit$pi0,
                                    "mu" = Stool$EM_fit$fit$mu,
                                    "sigma" = Stool$EM_fit$fit$sigma)
    feature_param <-
      SparseDOSSA2:::generate_featureParam(
        F_fit = Stool$F_fit,
        param_original = feature_param_template,
        n_feature = d)

    # Change pi = 0 to remove generation of inflated zero
    Stool_mod$EM_fit$fit$pi0 <- 0 #feature_param[, "pi0"]
    Stool_mod$EM_fit$fit$mu <- feature_param[, "mu"]
    Stool_mod$EM_fit$fit$sigma <- feature_param[, "sigma"]


    set.seed(seed.num)
    # Just to get names of taxa
    sim <- SparseDOSSA2(template= Stool_mod,
                        n_sample=n, n_feature = d,
                        new_features=T,
                        median_read_depth=50000)

    num.spiked <- pi * d

    if(j == 1){
      # balanced
      lfcs <- 5*rbeta((num.spiked/2),1,3)+1
      lfcs <- c(lfcs, -lfcs)
      subset_features <- sample(rownames(sim$params$feature_param), length(lfcs))

    }else if(j == 2){
      # imbalanced : all increasing
      lfcs <- 5*rbeta(num.spiked,1,3)+1
      lfcs <- abs(lfcs)
      subset_features <- sample(rownames(sim$params$feature_param), length(lfcs))

    }else if(j == 3){
      # imbalanced ; 4/5 increasing 1/5 decreasing
      lfcs <- 5*rbeta(num.spiked,1,3)+1
      lfcs <- lfcs * c(-1,1,1,1,1)
      subset_features <- sample(rownames(sim$params$feature_param), length(lfcs))
    }

    mat_metadata <- as.matrix(data.frame(disease=c(rep(0, n/2), rep(1, n/2))))
    metadata_df <- data.frame(
      metadata_datum=rep(1, length(subset_features)),
      feature_spiked=subset_features,
      associated_property=rep("abundance", length(subset_features)),
      effect_size=lfcs)

    set.seed(seed.num)
    sim <- SparseDOSSA2(template= Stool_mod,
                        n_sample=n, n_feature = d,
                        new_features=T,
                        metadata_matrix=mat_metadata,
                        spike_metadata=metadata_df,
                        median_read_depth= lambda)

    W <- t(sim$simulated_data+1)
    Y <- mat_metadata

    out <- matrix(NA, nrow = 2, ncol = 4)
    colnames(out) <- c("target FDR", "num_selected", "FDP_S", "power")
    out[,1] <- c(0.1, 0.05)
    l <- 1
    for(q in c(0.1, 0.05)){
      resultSingle = MRCKF(W, Y, q, penalty = "lasso")
      selectionSetS = resultSingle$S

      num_selected <- length(selectionSetS)
      FDP.S <- sum(!selectionSetS %in% metadata_df$feature_spiked)/length(selectionSetS)
      power <- (num_selected-sum(!selectionSetS %in% metadata_df$feature_spiked))/num.spiked

      out[l,2:4] <- cbind(num_selected = num_selected, FDP_S = FDP.S, power=power)
      l <- l+1
    }

    Y <- sim$simulated_data
    S <- 2000
    D <- nrow(Y)
    N <- ncol(Y)
    X <- mat_metadata[,1]

    resam.array <- array(NA, dim = c(S,D,N))
    for(n in 1:N){
      dirichlet_samples <- rdirichlet(S, Y[,n] + 0.5)
      resam.array[, , n] <- (log2(dirichlet_samples))
    }

    ## tryCatch: the current sparse_auto() stops when the laplace arm keeps fewer
    ## than 2 draws, which would otherwise kill the whole 1000-replicate run.
    sp <- tryCatch({
      sparse_out <- sparse_auto(dat.array = resam.array, Y, X, D, N, S,
                                filtered.ind = rep(FALSE, D))
      truth <- metadata_df$feature_spiked
      arms  <- sparse_arms(sparse_out)

      sb <- score_arm(arms$bootstrap, truth, num.spiked)
      sl <- score_arm(arms$laplace,   truth, num.spiked)

      list(sparse = score_arm(arms$hybrid, truth, num.spiked),
           rec = data.frame(scale_model    = attr(sparse_out, "scale_model"),
                            n_laplace_fail = attr(sparse_out, "n_laplace_fail"),
                            TP_boot = sb["TP"], FP_boot = sb["FP"],
                            tpr_boot = sb["tpr"], fdr_boot = sb["fdr"],
                            TP_lap  = sl["TP"], FP_lap  = sl["FP"],
                            tpr_lap  = sl["tpr"], fdr_lap  = sl["fdr"],
                            row.names = NULL, stringsAsFactors = FALSE))
    }, error = function(e) list(sparse = rep(NA_real_, 4), rec = ARM_NA))

    return(list(out = out, sparse = sp$sparse, rec = sp$rec))
  }

  for (m in 1:1000) {
    out.list_j[[m]] <- out.list[[m]]$out
    sparse_out_alltypes[[j]][m, ] <- out.list[[m]]$sparse
  }

  arm_rec <- do.call(rbind, lapply(out.list, `[[`, "rec"))
  arm_rec <- cbind(effect_dist = names(sparse_out_alltypes)[j],
                   rep = seq_len(nrow(arm_rec)), arm_rec)
  arm_rec_alltypes[[j]] <- arm_rec

  if(j == 1){
    out_balanced <- out.list_j
    saveRDS(out_balanced, file.path(RES, paste0("out_balanced_", TAG, ".RData")))

  }else if(j == 2){
    out_allincrease <- out.list_j
    saveRDS(out_allincrease, file.path(RES, paste0("out_allincrease_", TAG, ".RData")))

  }else if(j == 3){
    out_1dec4inc <- out.list_j
    saveRDS(out_1dec4inc, file.path(RES, paste0("out_1dec4inc_", TAG, ".RData")))
  }
}

stopCluster(cl)


df_sparse <- as.data.frame(rbind(cbind(sparse_out_alltypes[[1]], "Symmetric (Sum = 0)"),
                                 cbind(sparse_out_alltypes[[2]], "All Positive"),
                                 cbind(sparse_out_alltypes[[3]], "80% Positive")))
colnames(df_sparse) <- c("TP", "FP", "power", "fdr", "effect_dist")
saveRDS(df_sparse, file.path(RES, paste0("df_sparse_", TAG, ".RData")))

df_sparse$fdr <- as.numeric(df_sparse$fdr)
df_sparse$power <- as.numeric(df_sparse$power)

df_sparse_summary <- df_sparse %>%
  group_by(effect_dist) %>%
  summarize(fdr_mean = mean(fdr, na.rm = TRUE),
            power_mean = mean(power, na.rm = TRUE),
            n_fdr = sum(!is.na(fdr)))
df_sparse_summary


## per-replicate scale-model record -------------------------------------------
df_arm <- do.call(rbind, arm_rec_alltypes)
saveRDS(df_arm, file.path(RES, paste0("sparse_arm_record_", TAG, ".RData")))
write.csv(df_arm, file.path(RES, paste0("sparse_arm_record_", TAG, ".csv")),
          row.names = FALSE)

df_arm_summary <- df_arm %>%
  group_by(effect_dist) %>%
  summarize(prop_bootstrap = mean(scale_model == "bootstrap", na.rm = TRUE),
            prop_laplace   = mean(scale_model == "laplace",   na.rm = TRUE),
            fdr_bootstrap  = mean(fdr_boot, na.rm = TRUE),
            tpr_bootstrap  = mean(tpr_boot, na.rm = TRUE),
            fdr_laplace    = mean(fdr_lap,  na.rm = TRUE),
            tpr_laplace    = mean(tpr_lap,  na.rm = TRUE),
            n_fdr_boot     = sum(!is.na(fdr_boot)),
            n_fdr_lap      = sum(!is.na(fdr_lap)),
            n_laplace_fail = mean(n_laplace_fail, na.rm = TRUE))
df_arm_summary

## Both FDR / Power tables. Reads the RDS written above, so it can also be run
## on its own to rebuild the tables without refitting.
source(file.path(ROOT, "ckf_tables_09012026.R"))
