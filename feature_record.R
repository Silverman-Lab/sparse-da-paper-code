################################################################################
## Feature-level bookkeeping, shared by the per-dataset analysis scripts.
##
## The analysis scripts score the methods (TP/FP against a scale-aware gold
## standard). That is a COUNT. The identity of the significant taxa survives only
## as bare feature ids scattered across a dozen differently shaped tables, which
## is enough to rank methods and not enough to say anything biological.
##
## These helpers pull the ids out of whatever table a method happened to
## produce, attach the lineage, and reshape everything into one long record plus
## a feature x method incidence matrix.
##
## Sourced by:
##   2015_prasse_appliedEnvMicrobiology_wetlands/sparse_DA_03212026.R
##   2022_Jin_natureComm_technicalReplicates/sparse_DA_03212026.R
##
## Each script supplies its own EFFECT_COL, because both the method names and the
## contrast column names are dataset-specific (ANCOM-BC reports "TYPEREST" for
## Prasse and "dietVS" for Jin).
################################################################################

## Lineage for a vector of feature ids. tax_mat is the taxonomy keyed by the same
## ids used everywhere downstream; ids with no match come back as NA rows rather
## than silently dropping out of the record.
lineage_of <- function(ids, tax_mat){
  ids <- as.character(ids)
  m <- as.data.frame(tax_mat[match(ids, rownames(tax_mat)), , drop = FALSE],
                     stringsAsFactors = FALSE)
  rownames(m) <- NULL
  cbind(category = ids, m, stringsAsFactors = FALSE)
}

## Some pipelines (RDP-style classifiers, as in Jin) glue the bootstrap
## confidence onto the name -- "Vampirovibrio(0.71)" -- and write unassigned
## ranks as "*(-)". Left alone those strings do not group, do not sort and do not
## match anything external. This splits each rank into a bare name and its
## confidence, and turns "*" into NA. Taxonomies that are already clean (Prasse)
## pass through unchanged apart from gaining all-NA confidence columns, so it is
## safe to call unconditionally.
clean_taxonomy <- function(tax_mat, ranks = colnames(tax_mat), keep_conf = TRUE){
  tax_mat <- as.matrix(tax_mat)
  nm <- tax_mat[, ranks, drop = FALSE]
  cf <- nm
  has_paren <- !is.na(nm) & grepl("\\([^()]*\\)$", nm)
  cf[!has_paren] <- NA_character_
  cf[has_paren]  <- sub(".*\\(([^()]*)\\)$", "\\1", nm[has_paren])
  nm[has_paren]  <- sub("\\([^()]*\\)$", "", nm[has_paren])
  nm[!is.na(nm) & nm %in% c("*", "")] <- NA_character_   # unassigned rank
  cf <- suppressWarnings(apply(cf, 2, as.numeric))       # "-" -> NA
  out <- as.data.frame(nm, stringsAsFactors = FALSE)
  if (keep_conf) {
    cf <- as.data.frame(cf, stringsAsFactors = FALSE)
    colnames(cf) <- paste0(ranks, "_conf")
    out <- cbind(out, cf)
  }
  ## data.frame, not matrix: as.matrix() would coerce the numeric confidences
  ## back to strings. lineage_of() indexes either kind by rowname, so this stays
  ## a drop-in replacement for the raw taxonomy table.
  rownames(out) <- rownames(tax_mat)
  out
}

## The lowest rank actually assigned for each row, plus which rank that was --
## the label you want on a figure when Genus is NA for half the taxa.
lowest_assigned <- function(tax_df, ranks){
  ranks <- ranks[ranks %in% colnames(tax_df)]
  lab <- rep(NA_character_, nrow(tax_df)); at <- rep(NA_character_, nrow(tax_df))
  for (r in ranks) {                       # ranks ordered broad -> fine
    v <- as.character(tax_df[[r]])
    ok <- !is.na(v) & nzchar(v)
    lab[ok] <- v[ok]; at[ok] <- r
  }
  data.frame(taxon_label = lab, taxon_rank = at, stringsAsFactors = FALSE)
}

## One row per (method, significant feature), lineage attached. Works unchanged
## on the per-method tables and on the gold-standard tables -- anything that is a
## data.frame with a `category` column.
##
## Each method's own statistic is recorded VERBATIM and labelled with its source
## column, never rescaled or sign-flipped, because the sign conventions differ
## between packages: DESeq2/limma/LinDA contrast the treatment group against the
## reference, ANCOM-BC may report the reverse, aldex.lfc() computes
## levels[1] - levels[2], and sparse_SSRV codes its design off the caller's X_0.
## Comparing `stat` across methods is therefore meaningless -- use the lfc_prop /
## lfc_abs columns the scripts compute once from the counts, under a single
## documented convention, for direction.
sig_features_long <- function(list_tax, tax_mat, effect_col = NULL){
  eff_fallback  <- c("lfc", "log2FoldChange", "logFC", "effect", "mean")
  padj_fallback <- c("padj", "adj.P.Val", "FDR", "FWER.DE", "q.val", "q_val", "we.eBH")
  pick <- function(tab, cands){
    hit <- cands[cands %in% colnames(tab)]
    if (length(hit)) hit[1] else NA_character_
  }

  out <- lapply(names(list_tax), function(nm){
    tab <- as.data.frame(list_tax[[nm]], stringsAsFactors = FALSE, check.names = FALSE)
    ## A method that ran and found nothing is not an error: it contributes no
    ## rows here and an all-zero column in the incidence matrix.
    if (nrow(tab) == 0 || is.null(tab[["category"]])) return(NULL)
    ## An explicit NA in effect_col means "this method reports no effect size"
    ## (ANCOM reports detection only; baySeq reports a direction string and a log
    ## posterior, neither of which is an effect size) -- respect it rather than
    ## letting the fallback invent one.
    named <- !is.null(effect_col) && nm %in% names(effect_col)
    ec <- if (named) unname(effect_col[[nm]]) else pick(tab, eff_fallback)
    if (!is.na(ec) && !(ec %in% colnames(tab))) ec <- pick(tab, eff_fallback)
    pc <- pick(tab, padj_fallback)
    data.frame(
      Method   = nm,
      lineage_of(tab[["category"]], tax_mat),
      stat     = if (is.na(ec)) NA_real_ else suppressWarnings(as.numeric(tab[[ec]])),
      stat_col = if (is.na(ec)) NA_character_ else ec,
      padj     = if (is.na(pc)) NA_real_ else suppressWarnings(as.numeric(tab[[pc]])),
      padj_col = if (is.na(pc)) NA_character_ else pc,
      stringsAsFactors = FALSE, check.names = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) return(out)
  rownames(out) <- NULL
  out
}

## feature x method incidence (1 = called significant), with the lineage and the
## reference log fold changes carried along. method_levels fixes the column set
## and its order, so a method that found nothing still appears as a zero column
## instead of vanishing from the table.
sig_features_matrix <- function(sig_long, tax_lookup, method_levels){
  ids <- unique(sig_long$category)
  inc <- matrix(0L, nrow = length(ids), ncol = length(method_levels),
                dimnames = list(ids, method_levels))
  hit <- cbind(match(sig_long$category, ids), match(sig_long$Method, method_levels))
  hit <- hit[!is.na(hit[, 1]) & !is.na(hit[, 2]), , drop = FALSE]
  inc[hit] <- 1L
  out <- cbind(tax_lookup[match(ids, tax_lookup$category), , drop = FALSE],
               as.data.frame(inc, check.names = FALSE),
               stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
