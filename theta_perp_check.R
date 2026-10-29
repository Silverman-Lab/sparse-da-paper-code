################################################################################
## theta_perp_check.R -- recompute ONLY theta^perp for the three studies with
## external quantification, without re-running the full DA pipelines.
##
## Why this exists. The DA scripts now record theta^perp themselves (see
## record_theta_perp() in scale_model_record.R), but each of them also fits
## twelve competitor methods, and two of them shell out to Julia for MiMIX. When
## the question is only "what scale does the method estimate, and how does it
## compare with the measurement?", that is hours of work for a number that takes
## minutes. This file is the short path: same data prep, same Dirichlet
## resampling, same sparse_auto(), nothing else.
##
##   Rscript theta_perp_check.R              # all three
##   PERP_WHICH=zemb Rscript theta_perp_check.R
##   PERP_OUT=/somewhere Rscript theta_perp_check.R
##
## Output goes to ROOT/theta_perp_check/ by DEFAULT, deliberately NOT into the
## study directories: those hold the records the DA scripts write, and a file
## from this script sitting among them would be indistinguishable from one
## produced by a full run. Point PERP_OUT at a study directory only if you mean
## to overwrite its record.
##
## These numbers are the same posterior sampled independently, not draw-for-draw
## identical to a full DA run: the full scripts consume random numbers in ALDEx2
## and elsewhere before reaching sparse_auto(), so the streams diverge. The
## Monte Carlo error of the mean is reported (mc_se) and is small next to the
## posterior sd -- but the DA scripts' own records remain the canonical ones.
##
## The data prep below is copied from each DA script rather than sourced, because
## those scripts fit models as a side effect of being sourced. Copying means it
## can drift when a DA script's filter changes -- and it did: Zemb moved from
## `mean > 5` to 10% prevalence, and this file went on estimating theta^perp for
## 198 features while the analysis fitted 411. A comment recording the expected
## dimensions was not enough to catch that, so fit_perp() now CHECKS D and N
## against the study's scale_model_record.rds and refuses to record a run built
## on a different feature set.
################################################################################

ROOT <- local({  # manuscript_code: the folder holding .manuscript_code_root (README.md)
  p <- normalizePath(getwd(), mustWork = FALSE)
  while (!file.exists(file.path(p, ".manuscript_code_root"))) {
    if (dirname(p) == p) stop("run from inside the manuscript_code folder (see README.md)")
    p <- dirname(p)
  }
  p
})
suppressMessages({
  library(phyloseq); library(dplyr); library(tidyverse); library(data.table)
  library(stringr);  library(dirmult)
})
source(file.path(ROOT, "sparse_SSRV.R"))
source(file.path(ROOT, "scale_model_record.R"))

OUT   <- Sys.getenv("PERP_OUT",   file.path(ROOT, "theta_perp_check"))
WHICH <- Sys.getenv("PERP_WHICH", "all")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(1234)

## One study (or one level/location within a study): build the Dirichlet array
## exactly as the DA script does, fit, record. `obs` must be a DENSITY -- total
## abundance per unit sample mass -- see the note above record_theta_perp().
## Studies whose prep no longer matches their DA script are SKIPPED, not fatal.
##
## One study halting the run used to cost the other four: the Zemb block sits
## first, so when its filter drifted, Vandeputte, Jin and Prasse never ran and
## theta_perp_check.csv kept whatever it held from the previous run -- silently,
## because the file was still there. Each drifted study is now reported at the
## point of failure, listed again at the end, and the script exits non-zero, so
## nothing is quietly skipped and the studies that ARE in step still get fitted.
##
## Set PERP_STRICT=1 to restore the halt-on-first-drift behaviour.
.perp_skipped <- character(0)
PERP_STRICT   <- nzchar(Sys.getenv("PERP_STRICT"))

fit_perp <- function(Y0, X0, S, dataset, pair = NA_character_, obs = NA_real_,
                     matched = NULL, notes = NA) {
  D <- nrow(Y0); N <- ncol(Y0)
  tag <- paste0(dataset, if (is.na(pair)) "" else paste0("/", pair))
  ## The feature set must be the one the DA script fitted, or theta^perp here
  ## describes a different model from the one every table reports.
  rec_f <- file.path(ROOT, "real_data_scale_model.RData")
  if (file.exists(rec_f)) {
    rec <- readRDS(rec_f)
    hit <- rec[rec$dataset == dataset &
               ((is.na(pair) & is.na(rec$pair)) |
                (!is.na(pair) & !is.na(rec$pair) & rec$pair == pair)), ]
    if (nrow(hit) && (hit$D[1] != D || hit$N[1] != N)) {
      msg <- sprintf("%s: this file builds %d x %d but the analysed run is %d x %d -- the data prep here has drifted from that study's DA script, or that script did not finish writing its record.",
                     tag, D, N, hit$D[1], hit$N[1])
      if (PERP_STRICT) stop(msg)
      .perp_skipped <<- c(.perp_skipped, msg)
      warning(msg, "\n  SKIPPED -- no theta^perp recorded for it this run.",
              call. = FALSE, immediate. = TRUE)
      return(invisible(NULL))
    }
  }
  cat(sprintf("\n### %s%s : D=%d N=%d S=%d  (n1=%d n0=%d)\n", dataset,
              if (is.na(pair)) "" else paste0(" / ", pair), D, N, S, sum(X0 == 1), sum(X0 == 0)))
  arr <- array(NA_real_, dim = c(S, D, N))
  for (n in 1:N) arr[,,n] <- log2(rdirichlet(S, as.vector(Y0[,n] + 0.5)))
  o <- sparse_auto(dat.array = arr, Y0, X0, D, N, S, filtered.ind = rep(FALSE, D))
  rm(arr); gc(verbose = FALSE)
  record_theta_perp(o, dataset = dataset, dir = OUT, pair = pair,
                    S = S, N = N, obs = obs, matched = matched, notes = notes)
}

## ------------------------------------------------------------------- ZEMB ---
if (WHICH %in% c("all", "zemb")) {
  ## Same parser the DA scripts use, so the sample set cannot drift between the
  ## two. It drops the no-feces control tube (oz1802-19) at the source.
  source(file.path(ROOT, "zemb_parse.R"))
  zemb     <- parse_2020_zemb_microbiologyOpen_spike()

  ## SAMPLE DROP, matching ZEMB_DROP_SAMPLES in sparse_DA_09012026.R: the 100 ul
  ## extreme-dose tube. Counts, metadata and qPCR are subset TOGETHER, because
  ## `scale` below is positional -- dropping a sample from one table only is
  ## exactly how the scale and the counts come apart without saying so.
  ## character(0) since 2026-09-01, matching sparse_DA_09012026.R: the 100 ul
  ## extreme-dose tube is back in and this is the full 9-sample design. This
  ## file must follow that script exactly -- it estimates theta^perp FOR that
  ## fit, so a sample set that disagrees is silently a different quantity.
  ZEMB_DROP_SAMPLES <- character(0)
  if (length(ZEMB_DROP_SAMPLES)) {
    stopifnot(ZEMB_DROP_SAMPLES %in% colnames(zemb$counts))
    keep_s <- !(colnames(zemb$counts) %in% ZEMB_DROP_SAMPLES)
    message("  [Zemb] dropping extreme-dose sample(s): ",
            paste(ZEMB_DROP_SAMPLES, collapse = ", "),
            "  (", ncol(zemb$counts), " -> ", sum(keep_s), " samples)")
    zemb$counts   <- zemb$counts[ , keep_s, drop = FALSE]
    zemb$metadata <- zemb$metadata[keep_s, , drop = FALSE]
    zemb$scale    <- zemb$scale[   keep_s, , drop = FALSE]
  }

  metadata <- cbind(zemb$metadata, spiked = zemb$metadata$added_Coli > 0)
  phylo <- phyloseq(otu_table(t(zemb$counts), taxa_are_rows = FALSE), sample_data(metadata))

  ## FEATURE FILTER, matching sparse_DA_09012026.R: rowMeans(raw counts) >= 5.
  ##
  ## This block has had to follow FIVE rules now -- `mean > 5`, 10% prevalence,
  ## 0.05% relative abundance, >=5 counts in >=2 samples, and back to `mean >= 5`
  ## -- and every time it lagged the DA script by one change it estimated
  ## theta^perp on a different feature set than the analysis fitted (209 against
  ## 256, the last time). If the DA script's filter moves again, THIS IS THE
  ## FIRST PLACE TO FOLLOW IT.
  ##
  ## Counts, not proportions -- do NOT reintroduce transform_sample_counts() here.
  ## The DA script thresholds the raw otu_table, and a relative threshold keeps a
  ## different feature set in a shallow sample than an absolute one does.
  ##
  ## Seq135 is the spiked-in internal standard, not part of the community, and
  ## the DA script removes it with setdiff() after filtering. Removed here too,
  ## or the feature counts cannot match.
  MIN_MEAN <- 5
  Mc <- as(otu_table(phylo), "matrix")
  if (!taxa_are_rows(phylo)) Mc <- t(Mc)             # taxa x samples
  keep_ids <- taxa_names(phylo)[rowMeans(Mc) >= MIN_MEAN]
  keep_ids <- setdiff(keep_ids, "Seq135")
  red <- prune_taxa(keep_ids, phylo)
  message(sprintf("  [Zemb] mean count >= %g: kept %d of %d taxa",
                  MIN_MEAN, ntaxa(red), ntaxa(phylo)))

  Y0 <- t(otu_table(red)); X0 <- as.numeric(sample_data(red)$spiked == "TRUE")

  ## DENSITY, not copies per tube -- this contrast loads 11-14 mg into the spiked
  ## tubes and 38-55 mg into the unspiked ones (they come from a different
  ## experiment in the source paper), so the per-tube total measures the pipetting.
  dens <- zemb$scale$total16S_per_mg        # divided by mgfeces in the parser
  ## Named, so perp_matched_obs() can check the measurement against the counts by
  ## sample name instead of trusting position. The parser already asserts
  ## identical(colnames(counts), rownames(scale)), which is what makes the names
  ## meaningful here rather than decorative.
  names(dens) <- rownames(zemb$scale)
  stopifnot(length(dens) == ncol(Y0), all(is.finite(dens)))

  ## theta^perp as fitted is the log2 fold change of the total load of the KEPT
  ## features, so the qPCR -- which counts the whole community -- has to be
  ## projected onto them before bias and coverage mean anything. Seq135 is out of
  ## the denominator as well as the numerator: it is spiked-in standard DNA, not
  ## community, and the qPCR total is derived from it rather than including it.
  ##
  ## Zemb is the study where this correction does almost nothing -- r_n is
  ## 0.98-0.99, worth -0.005 log2 -- and that is the point of running it here:
  ## a filter that drops rare features drops almost no reads with them. What
  ## moves the number on the other studies is tax_glom(), which Zemb has none of.
  M_all_z <- as(otu_table(phylo), "matrix")            # samples x taxa
  M_all_z <- M_all_z[, setdiff(colnames(M_all_z), "Seq135"), drop = FALSE]
  m_zemb  <- perp_matched_obs(total = dens, X = X0, Y_all = M_all_z,
                              taxa_are_rows = FALSE, keep = keep_ids)
  fit_perp(Y0, X0, 50000, "Zemb2020",
           obs = mean(log2(dens[X0 == 1])) - mean(log2(dens[X0 == 0])),
           matched = m_zemb,
           notes = "qPCR: average16S copies per mg feces (density)")
}

## ------------------------------------------------------------- VANDEPUTTE ---
if (WHICH %in% c("all", "vandeputte")) {
  md <- read.csv(file.path(ROOT, "Vandeputte", "data", "41586_2017_BFnature24460_MOESM10_ESM.csv")) %>%
    dplyr::rename("CellCount" = "Average.cell.count..per.gram.of.frozen.feces.")
  OTU  <- fread(file.path(ROOT, "Vandeputte", "data", "OTU_nochim.csv"))[-c(1:40), ]
  taxa <- fread(file.path(ROOT, "Vandeputte", "data", "taxa_assignments.csv"), header = TRUE) %>%
    column_to_rownames("V1")
  taxa_id <- paste0("sp", 1:nrow(taxa)); rownames(taxa) <- taxa_id; taxa <- as.matrix(taxa)
  OTU <- as.data.frame(OTU) %>% column_to_rownames("V1")
  colnames(OTU) <- taxa_id; rownames(OTU) <- paste0("sa", 1:nrow(OTU))
  phylo <- phyloseq(otu_table(OTU, taxa_are_rows = FALSE), sample_data(md), tax_table(taxa))
  red <- filter_taxa(tax_glom(phylo, taxrank = "Genus"),
                     function(v) sum(v > 1) > .2 * length(v), prune = TRUE)
  Y0 <- t(otu_table(red)); X0 <- as.numeric(sample_data(red)$Health.status == "CD")
  ## already a density -- cells per GRAM of frozen faeces
  cc <- sample_data(red)$CellCount
  names(cc) <- sample_names(red)
  stopifnot(length(cc) == ncol(Y0))

  ## Y_kept, NOT keep = taxa_names(red). tax_glom() names each genus after its
  ## most abundant constituent OTU, so those 29 ids all resolve against `phylo`
  ## and select 29 OTUs -- where the model was fitted on 29 GENERA. The retained
  ## fraction would come out far too small, and nothing would say so.
  ##
  ## Y_all is the PRE-GLOM table: 2041 OTUs down to 195 genera, and the OTUs with
  ## no genus assignment that tax_glom() discards are real cells the cytometer
  ## counted. That choice is most of the correction here (glom -0.408 of it on
  ## Prasse, versus +0.080 for the prevalence filter alone).
  ##
  ## Caveat that belongs with this study specifically: r_n is a READ fraction,
  ## and the assay counts CELLS, so 16S copy-number variation and extraction bias
  ## enter the projection. Vandeputte is the weakest of the four on that footing.
  m_van <- perp_matched_obs(total = cc, X = X0,
                            Y_all = as(otu_table(phylo), "matrix"),
                            taxa_are_rows = FALSE, Y_kept = Y0,
                            kept_taxa_are_rows = TRUE)
  fit_perp(Y0, X0, 50000, "Vandeputte2017",
           obs = mean(log2(cc[X0 == 1])) - mean(log2(cc[X0 == 0])),
           matched = m_van,
           notes = "flow cytometry: cells per gram of frozen faeces (density)")
}

## -------------------------------------------------------------------- JIN ---
if (WHICH %in% c("all", "jin")) {
  local <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates/")
  dat <- readxl::read_xlsx(paste0(local, "Supplementary Data 8.xlsx"),
                           sheet = "Sequencing-determined counts", skip = 1)
  cnt <- as.matrix(dat[, 2:56]); rownames(cnt) <- dat$`cOTU-ID`
  scl <- readxl::read_xlsx(paste0(local, "Supplementary Data 8.xlsx"),
                           sheet = "Absolute total abundance", skip = 1)
  meta <- as.data.frame(matrix(NA, nrow = nrow(scl), ncol = 5))
  rownames(meta) <- scl$Sample
  colnames(meta) <- c("diet", "subject", "location", "sample_type", "technical_replicate")
  meta$diet <- substr(scl$Sample, 1, 2); meta$subject <- substr(scl$Sample, 3, 3)
  meta$location <- substr(scl$Sample, 4, 7)
  meta$sample_type <- unlist(lapply(strsplit(rownames(meta), "-"), function(z) z[2]))
  meta$sample_type[is.na(meta$sample_type)] <- "cell"
  cnt  <- cnt[, colnames(cnt) %in% rownames(meta)]
  meta <- meta[rownames(meta) %in% colnames(cnt), ]
  as_count_mat <- function(d, idcol) {
    m <- as.matrix(d[, setdiff(colnames(d), idcol), drop = FALSE]); mode(m) <- "numeric"
    m[is.na(m)] <- 0; rownames(m) <- d[[idcol]]; m
  }
  asv_counts <- as_count_mat(readxl::read_xlsx(paste0(local, "Supplementary Data 5.xlsx"),
                                               sheet = "VA group", skip = 1), "ASV-ID")
  otu_counts <- as_count_mat(readxl::read_xlsx(paste0(local, "Supplementary Data 7.xlsx"),
                                               sheet = "VA group", skip = 1), "OTU-ID")

  OTU_all <- t(cnt)
  metadata_all <- meta[match(rownames(OTU_all), rownames(meta)), ]
  ## by sample NAME, not position -- out$scale interleaves the -ecDNA rows
  si <- match(rownames(metadata_all), scl$Sample); stopifnot(!anyNA(si))
  ## already a density -- copies per mg
  scale_all <- scl$`Absolute total abudance (copies/mg)`[si]

  MIN_COUNT <- 5; MIN_SAMPLES <- 3
  for (level in c("cOTU", "ASV", "OTU")) {
    has_scale  <- identical(level, "cOTU")
    COUNTS_ALL <- if (has_scale) OTU_all else t(if (level == "ASV") asv_counts else otu_counts)
    ## "both" is the POOLED cOTU run -- the same eight mice at both cecal
    ## locations, 16 samples -- and it is the run `real data results.R` reports
    ## as the primary Jin result. Only cOTU has it, because only cOTU carries a
    ## measured scale; the key is the bare "both", matching the DA script's PAIR.
    for (loc in c(if (has_scale) "both", "prox", "dist")) {
      key <- if (has_scale) loc else paste0(level, "_", loc)
      keep_names <- rownames(metadata_all)[
        metadata_all$diet %in% c("VS", "VD") & metadata_all$sample_type == "cell" &
        (if (identical(loc, "both")) TRUE else metadata_all$location == loc)]
      keep_names <- intersect(keep_names, rownames(COUNTS_ALL))
      ph  <- phyloseq(otu_table(COUNTS_ALL[keep_names, , drop = FALSE], taxa_are_rows = FALSE),
                      sample_data(metadata_all[keep_names, , drop = FALSE]))
      red <- filter_taxa(ph, function(v) sum(v >= MIN_COUNT) >= MIN_SAMPLES, prune = TRUE)
      Y0 <- t(otu_table(red)); X0 <- as.numeric(sample_data(red)$diet == "VD")
      sk <- if (has_scale) scale_all[match(keep_names, rownames(metadata_all))] else NULL
      if (has_scale) names(sk) <- keep_names
      obs <- if (has_scale) {
        stopifnot(length(sk) == ncol(Y0))
        mean(log2(sk[X0 == 1])) - mean(log2(sk[X0 == 0]))
      } else NA_real_
      ## The BarBIQ total is copies/mg over the WHOLE community; theta^perp is
      ## fitted on the kept cOTUs only. `keep =` rather than `Y_kept =` because
      ## Jin does no agglomeration -- the fitted features are literally rows of
      ## `ph`, so the id lookup is both valid and a check that they are.
      ##
      ## Y_all is the sample-subset table, not COUNTS_ALL: the community whose
      ## fraction is being taken has to be the community in these samples.
      matched <- if (has_scale)
        perp_matched_obs(total = sk, X = X0,
                         Y_all = as(otu_table(ph), "matrix"), taxa_are_rows = FALSE,
                         keep = taxa_names(red)) else NULL
      fit_perp(Y0, X0, 10000, "Jin2022", pair = key, obs = obs, matched = matched,
               notes = if (has_scale) "BarBIQ absolute total abundance, copies/mg (density)"
                       else NA_character_)
    }
  }
}

## ----------------------------------------------------------------- PRASSE ---
##
## Prasse measured total 16S by qPCR (QPCR_16S), so theta^perp can be scored here
## as it is for Zemb, Vandeputte and Jin. Two things make this study different
## from those three, and both belong in the caption of anything that reports it.
##
## 1. IT IS NOT A DENSITY. The other three report total load per unit mass -- per
##    mg for Zemb and Jin, per gram for Vandeputte. Prasse reports copies PER
##    SAMPLE and gives no mass column to divide by, so its theta^perp is a
##    per-tube quantity and confounds load with how much soil went into the tube.
##    Recorded, but not comparable with the other three on the same axis.
##
## 2. THE QSeq PASS IS CIRCULAR. QSeq(TFW, "QPCR_16S") rescales every sample by
##    its own qPCR total, so the counts that pass feeds the model ALREADY carry
##    the measurement being scored against. Its bias against `obs` is therefore
##    not evidence the scale model recovered anything. The RAW pass is the honest
##    test: raw reads in, qPCR held out as an independent measurement. Both are
##    fitted so the difference is visible rather than argued.
##
## The prep below is sparse_DA_09012026_rawVSqseq_mean5.R's, down to the
## re-keying step: tax_glom() names each genus after its most abundant
## constituent OTU and therefore picks a different archetype before and after
## rescaling, so the QSeq genus table is re-keyed onto the RAW table's ids or the
## two passes carry different row names.
if (WHICH %in% c("all", "prasse")) {
  suppressMessages(library(QSeq))
  dp   <- file.path(ROOT, "2015_prasse_appliedEnvMicrobiology_wetlands")
  TFW  <- readRDS(file.path(dp, "Data", "TFW_comdat.RDS"))
  ## fileEncoding is required: the CSV begins with a UTF-8 BOM, so a plain
  ## read.csv() names the first column "X...ID" and meta$ID comes back NULL.
  meta <- read.csv(file.path(dp, "Data", "TFW_metadata.csv"), fileEncoding = "UTF-8-BOM")
  stopifnot(!is.null(meta$ID))
  rownames(meta) <- meta$ID
  sample_data(TFW) <- sample_data(meta)
  sample_names(TFW) <- paste0("Sample_", sample_names(TFW))

  Q.TFW <- QSeq(TFW, "QPCR_16S")
  stopifnot(identical(taxa_names(TFW), taxa_names(Q.TFW)),
            identical(sample_names(TFW), sample_names(Q.TFW)))

  taxa_p <- Q.TFW@tax_table
  rownames(taxa_p) <- paste0("sp", seq_len(nrow(taxa_p)))
  taxa_p <- as.matrix(taxa_p)
  OTU_qseq <- as.data.frame(t(Q.TFW@otu_table@.Data)); colnames(OTU_qseq) <- rownames(taxa_p)
  OTU_raw  <- as.data.frame(t(TFW@otu_table@.Data));   colnames(OTU_raw)  <- rownames(taxa_p)

  ## Kept as its own object rather than inlined into tax_glom(): it is the
  ## community the qPCR measured, and perp_matched_obs() needs it to work out how
  ## much of that community survived agglomeration. 22381 OTUs -> 561 genera.
  ph_raw_pre <- phyloseq(otu_table(OTU_raw, taxa_are_rows = FALSE),
                         Q.TFW@sam_data, tax_table(taxa_p))
  genus_qseq <- tax_glom(phyloseq(otu_table(OTU_qseq, taxa_are_rows = FALSE),
                                  Q.TFW@sam_data, tax_table(taxa_p)), taxrank = "Genus")
  genus_raw  <- tax_glom(ph_raw_pre, taxrank = "Genus")

  G_qseq <- as.character(tax_table(genus_qseq)[, "Genus"])
  G_raw  <- as.character(tax_table(genus_raw)[,  "Genus"])
  stopifnot(setequal(G_qseq, G_raw), !anyDuplicated(G_qseq), !anyDuplicated(G_raw))
  M_qseq <- as.matrix(otu_table(genus_qseq)); colnames(M_qseq) <- G_qseq
  M_qseq <- M_qseq[sample_names(genus_raw), G_raw, drop = FALSE]
  colnames(M_qseq) <- taxa_names(genus_raw)
  genus_qseq <- phyloseq(otu_table(M_qseq, taxa_are_rows = FALSE),
                         sample_data(genus_raw), tax_table(genus_raw))

  ## Filter computed ONCE on the RAW counts and applied to both passes, so the
  ## two differ only by the per-sample rescaling.
  ## rowMeans(raw counts) >= 5, matching sparse_DA_09012026_rawVSqseq_mean5.R.
  ## M_raw_chk is samples x taxa, so the per-taxon mean is colMeans().
  M_raw_chk <- as.matrix(otu_table(genus_raw))
  MIN_MEAN  <- 5
  keep_ids  <- taxa_names(genus_raw)[colMeans(M_raw_chk) >= MIN_MEAN]
  message(sprintf("  [Prasse] mean raw count >= %g: kept %d of %d genera",
                  MIN_MEAN, length(keep_ids), ntaxa(genus_raw)))

  for (INPUT in c("raw", "QSeq")) {
    ph <- prune_taxa(keep_ids, if (INPUT == "raw") genus_raw else genus_qseq)
    Y0 <- t(otu_table(ph))
    Xc <- as.character(sample_data(ph)$TYPE)
    stopifnot(identical(sort(unique(Xc)), c("NAT", "REST")))

    ## THE CONTRAST DIRECTION IS COPIED, NOT CHOSEN.
    ##
    ## The DA script codes the groups as `X_0 <- as.numeric(X == tail(X,1))` --
    ## group 1 is whichever TYPE the LAST sample happens to carry. In this
    ## phyloseq ordering that is NAT, so sparse_auto() is fitted on
    ## NAT - REST, and theta^perp here has to use the same coding or it would
    ## describe a model nobody fitted.
    ##
    ## Worth knowing before quoting a sign: the SAME DA script computes its
    ## feature-level lfc_prop / lfc_abs the other way round, from
    ## grpT = REST and grpC = NAT. So the posterior mean of theta and the
    ## lfc_abs column in prasse_tax_lookup_*.csv are opposite in sign, and
    ## neither is wrong on its own terms. Reproduced faithfully rather than
    ## quietly corrected -- flipping it here would only hide the disagreement.
    ##
    ## The assertion pins the coding: if a future re-ordering of the samples
    ## changes which TYPE lands last, the contrast silently reverses, and this
    ## is where that gets caught rather than in a sign in the manuscript.
    grp1 <- tail(Xc, 1)
    X0   <- as.numeric(Xc == grp1)
    stopifnot(identical(grp1, "NAT"))
    message(sprintf("  [Prasse/%s] contrast is %s - %s  (n1 = %d, n0 = %d), copied from the DA script",
                    INPUT, grp1, setdiff(c("NAT", "REST"), grp1), sum(X0 == 1), sum(X0 == 0)))

    qpcr <- as.numeric(sample_data(ph)$QPCR_16S)
    names(qpcr) <- sample_names(ph)
    stopifnot(length(qpcr) == ncol(Y0), all(is.finite(qpcr)), all(qpcr > 0))

    ## THE RETAINED FRACTION IS COMPUTED FROM RAW COUNTS IN BOTH PASSES.
    ##
    ## QSeq() rescales every sample by its own qPCR total -- one multiplier per
    ## SAMPLE, applied to every taxon in it -- so it cancels out of a within-sample
    ## ratio and r_n is identical for the two inputs. Taking the numerator from
    ## the QSeq table and the denominator from the raw OTU table would not cancel:
    ## it would divide rescaled counts by unrescaled ones and read the qPCR total
    ## itself as a retained fraction. `prasse_kept_raw` is therefore the RAW
    ## filtered table in both iterations of this loop.
    ##
    ## This is the study where the projection is largest: 22381 OTUs collapse to
    ## 561 genera and r_n falls to 0.117-0.342, worth -0.328 log2 -- enough to
    ## flip the sign of the bias. Almost all of it is tax_glom() discarding OTUs
    ## with no genus assignment, not the prevalence filter (+0.080 on its own),
    ## so the matched reference here inherits the annotation database's coverage.
    prasse_kept_raw <- t(otu_table(prune_taxa(keep_ids, genus_raw)))
    stopifnot(identical(colnames(prasse_kept_raw), sample_names(ph)))
    m_prasse <- perp_matched_obs(total = qpcr, X = X0,
                                 Y_all = as(otu_table(ph_raw_pre), "matrix"),
                                 taxa_are_rows = FALSE,
                                 Y_kept = prasse_kept_raw, kept_taxa_are_rows = TRUE)

    fit_perp(Y0, X0, 10000, "Prasse2015", pair = paste0("mean5_", INPUT),
             matched = m_prasse,
             ## Same orientation as X0, so bias and coverage compare like with like.
             obs = mean(log2(qpcr[X0 == 1])) - mean(log2(qpcr[X0 == 0])),
             notes = paste0("qPCR: 16S copies PER SAMPLE (not a density -- no mass ",
                            "column to divide by); contrast ", grp1, " - ",
                            setdiff(c("NAT", "REST"), grp1), "; input = ", INPUT,
                            if (INPUT == "QSeq") ", CIRCULAR: counts already rescaled by this same qPCR total" else ""))
  }
}

## --------------------------------------------------------------- STAMMLER ---
##
## Added 2026-09-05. This study was absent from this file with no note saying
## why, and its absence propagated: the scale-recovery table's caption claimed
## only four studies had a direct measurement of total load, when Stammler has
## one too and Task 6 was already drawing a line from it.
##
## THE MEASUREMENT IS PER TUBE, NOT A DENSITY. The study reports 16S rDNA copies
## per SAMPLE and there is no mass column to divide by, so theta^perp here is a
## per-tube quantity while Zemb, Jin and Vandeputte are per mg or per gram. For a
## defined mock community that is defensible -- the pools ARE the experimental
## quantity -- but it is not the same measurement as the others, and `notes`
## records that rather than leaving it to be inferred.
##
## SIXTEEN FITS, and the duplicate is deliberate. The DA script records the
## headline Pool1-vs-Pool6 comparison under the bare tag "Pool6" and the fifteen
## pool pairs under "<ref>_vs_<comp>". Pool1_vs_Pool6 is therefore the same
## comparison as Pool6, fitted twice under two names. Both are produced here
## because both tags are read downstream: the main scale-recovery table looks up
## the run DATASET_META names ("Pool6") and the supplementary one enumerates the
## pairs.
if (WHICH %in% c("all", "stammler")) {
  source(file.path(ROOT, "2016_stammler_microbiome_micehuman", "parse.R"), local = TRUE)
  st <- parse_2016_stammler_microbiome_micehuman()

  ## Matching MIN_MEAN in sparse_DA_09012026.R, and recomputed per comparison
  ## because a mean over samples depends on which samples are in the subset --
  ## the same reason the DA script's pair loop recomputes it.
  STAM_MIN_MEAN <- 5

  stam_fit <- function(ref_pool, target_pool, tag) {
    i_ref <- which(st$metadata$Treatment == ref_pool)
    i_tgt <- which(st$metadata$Treatment == target_pool)
    idx   <- c(i_ref, i_tgt)
    OTU   <- t(st$counts[, idx])                      # samples x taxa
    meta  <- st$metadata[idx, ]
    sc    <- as.numeric(st$scale$`16S rDNA copies per sample`[idx])
    ph    <- phyloseq(otu_table(OTU, taxa_are_rows = FALSE),
                      sample_data(meta), tax_table(as.matrix(st$tax)))

    genus <- tax_glom(ph, taxrank = "Genus")
    Mp <- as(otu_table(genus), "matrix")
    if (!taxa_are_rows(genus)) Mp <- t(Mp)            # taxa x samples
    red <- prune_taxa(taxa_names(genus)[rowMeans(Mp) >= STAM_MIN_MEAN], genus)

    Y0 <- as(otu_table(red), "matrix")
    if (!taxa_are_rows(red)) Y0 <- t(Y0)              # taxa x samples
    ## Contrast copied from the DA script: X_0 <- as.numeric(X == target_pool),
    ## so group 1 is the COMPARISON pool and the estimate is comp - ref. Explicit
    ## there, unlike Prasse's tail(X, 1), so there is nothing to pin with an
    ## assertion beyond the pools being the two expected ones.
    Xc <- as.character(sample_data(red)$Treatment)
    stopifnot(setequal(unique(Xc), c(ref_pool, target_pool)))
    X0 <- as.numeric(Xc == target_pool)
    stopifnot(length(sc) == ncol(Y0), all(is.finite(sc)), all(sc > 0))

    m_stam <- perp_matched_obs(total = sc, X = X0,
                               Y_all = OTU, taxa_are_rows = FALSE,
                               Y_kept = Y0, kept_taxa_are_rows = TRUE)

    fit_perp(Y0, X0, 10000, "Stammler2016", pair = tag, matched = m_stam,
             ## Same orientation as X0, so bias and coverage compare like with like.
             obs = mean(log2(sc[X0 == 1])) - mean(log2(sc[X0 == 0])),
             notes = paste0("16S rDNA copies PER SAMPLE (not a density -- no mass ",
                            "column to divide by); contrast ", target_pool, " - ", ref_pool))
  }

  ## The pools, and every unordered pair of them, derived here rather than taken
  ## from task_common.R: this file deliberately sources only sparse_SSRV.R and
  ## scale_model_record.R, and pulling in the task registry to read one constant
  ## would drag the whole dataset registry and its consistency checks in with it.
  ## combn() over the sorted pools reproduces STAMMLER_PAIRS exactly.
  ## Pool* only: the metadata also carries a SpikeKO and a noSpike tube, one
  ## sample each, which are controls rather than members of the dose series and
  ## are not part of any pair the DA script fits.
  stam_pools <- sort(grep("^Pool[0-9]+$",
                          unique(as.character(st$metadata$Treatment)), value = TRUE))
  stopifnot(identical(stam_pools, paste0("Pool", 1:6)))

  ## The headline comparison, under the tag the DA script and DATASET_META use.
  stam_fit("Pool1", "Pool6", "Pool6")
  ## Then every pool pair, keyed as the DA script keys them.
  for (p in combn(stam_pools, 2, simplify = FALSE))
    stam_fit(p[1], p[2], paste0(p[1], "_vs_", p[2]))
}

## ---------------------------------------------------------------- SUMMARY ---
rec <- collect_perp_records(OUT)
if (!is.null(rec)) {
  cat("\n\n########## estimated theta^perp -- all runs, all arms ##########\n")
  cat("## log2 fold change of TOTAL DENSITY, treatment - control\n")
  cat("## obs   = the assay's own value, over the WHOLE community\n")
  cat("## obs_S = the same value projected onto the features fitted, which is\n")
  cat("##         what theta^perp estimates; obs_S = obs + dlog2_r, exactly\n\n")
  print(rec[, c("dataset", "pair", "arm", "selected", "mean", "sd", "lower", "upper",
                "obs", "bias", "covered")], row.names = FALSE, digits = 4)
  cat("\n## the same runs against the filter-matched reference\n\n")
  print(rec[, c("dataset", "pair", "arm", "D_all", "D_kept", "r_min", "r_max",
                "dlog2_r", "obs_S", "bias_S", "covered_S")], row.names = FALSE, digits = 4)
  s <- score_theta_perp(rec)
  if (!is.null(s)) {
    cat("\n## accuracy against the external measurement, by scale model\n\n")
    print(s, row.names = FALSE, digits = 4)
    cat("\n## the same, per study\n\n")
    print(score_theta_perp(rec, by = c("dataset", "arm")), row.names = FALSE, digits = 4)
  }
  write.csv(rec, file.path(OUT, "theta_perp_check.csv"), row.names = FALSE)
  cat("\nwritten: ", file.path(OUT, "theta_perp_check.csv"), "\n", sep = "")
}

## Everything that could be fitted has been, and the CSV above is written, so a
## drifted study costs only its own row. Failing here rather than earlier keeps
## the exit code honest: a run that skipped a study is not a successful run, and
## the record it just wrote is missing that study.
if (length(.perp_skipped)) {
  stop("theta^perp was NOT recorded for ", length(.perp_skipped), " run(s):\n  ",
       paste(.perp_skipped, collapse = "\n  "),
       "\nEither re-run those DA scripts to completion, or bring the prep in this ",
       "file back in step with them. Every other run above was fitted and saved.")
}
