################################################################################
## zemb_parse.R -- THE parser for Zemb et al. (2020), in one place.
##
## Why this file exists. The same 15-line function was copy-pasted into six
## scripts (2020_zemb_microOpen_spike/{parse,sparse_DA,sparse_DA_02252025,
## sparse_DA_03212026,scale_samples}.R, real data LFC diagnosis*.R,
## real_data_scale_model_choice.R, real_data_scale_model_diagnosis.R), differing
## only in how the directory was spelled. That is fine until the parse has to
## CHANGE -- and it does, see below -- at which point "fix the parse function"
## means finding six of them and getting all six right.
##
## Every live script now sources this instead. Sourcing it defines
## parse_2020_zemb_microbiologyOpen_spike() with the same name, the same return
## shape (list of scale / metadata / counts) and the same no-argument call, so
## the call sites did not have to change.
##
##------------------------------------------------------------------------------
## WHAT CHANGED: the no-feces control tubes are dropped here.
##
## The metadata has 24 tubes; 12 were sequenced. Zemb et al. (2020,
## MicrobiologyOpen 9:e977, doi:10.1002/mbo3.977) describe two experiments:
##
##   Experiment 2  the E. coli series -- 15 tubes of 9.3-16 mg feces, doses
##                 0/1/5/10/100 ul of a 10^7 cells/ul suspension in triplicate,
##                 internal standard 2.7x10^7 copies per tube, ONE replicate per
##                 dose sequenced.                        -> tubes 1, 4, 7, 10, 13
##   Experiment 1  varying the STANDARD -- 2.8e6, 2.8e7, 1.4e8, 2.8e8 copies on
##                 38, 55, 38 and 43 mg of feces.         -> tubes 21, 22, 24, 23
##
## Tubes 16-20 belong to neither and appear in no published figure. They are lab
## controls, and the metadata says so without needing the paper: 16-18 carry
## feces with added_std = NA (no-standard controls), 19 is 0 mg feces + no
## standard + 10^7 E. coli (spike-only), 20 is 0 mg feces + standard + no E. coli
## (standard-only).
##
## Of those, only tube 19 was sequenced, and it is not a sample of anything: it is
## the spike suspension put through the pipeline on its own. 97.2% of its reads
## are Seq1 -- the spiked E. coli -- against a ~21% top-taxon share in every real
## sample. Its qPCR is undefined BY CONSTRUCTION and says so: Std_in_tube and
## average16S(copies/mg) are both #DIV/0! and efficiency is NA, because there is
## no internal standard to recover and no mass to normalise by.
##
## A tube with no feces has no microbial density, so it can be neither scored
## against a measured total nor allowed to set the reference for taxa. Dropping it
## HERE rather than in each analysis means counts, metadata and the qPCR table are
## subset TOGETHER at the source -- the three are aligned by construction and no
## downstream script has to keep a positional index in step by hand, which was the
## real hazard of dropping it later.
##
## The rule is `mgfeces <= 0`, not a hard-coded tube name: that is the actual
## criterion (no feces -> no density -> not a sample) and it stays correct if the
## sequenced set ever changes. On the current data it removes exactly oz1802-19.
## Pass drop_controls = FALSE to get all 12 sequenced tubes back.
##
##------------------------------------------------------------------------------
## STILL OPEN, and deliberately NOT decided here. The contrast the DA scripts run
## is `added_Coli > 0`, and its two arms come from the two DIFFERENT experiments
## above: TRUE is tubes 4/7/10/13 (Experiment 2), FALSE is tube 1 (Experiment 2,
## dose 0) plus tubes 21-24 (Experiment 1, 38-55 mg). That is where the 3-5 fold
## feces-mass difference between the arms comes from. The within-Experiment-2
## contrast would be tube 1 versus tubes 4/7/10/13 -- 1 vs 4, all at 11-14 mg with
## the same standard. Changing that is an analysis decision, not a parsing one.
################################################################################

library(stringr)

## Where the data lives. Auto-detected on the same pattern as `real data
## results.R`, first existing candidate wins, so this file works unedited from a
## cluster script and a laptop script -- and from scripts that spell their own
## root ROOT, PROJ, or not at all.
ZEMB_DIR_CANDIDATES <- file.path(local({  # manuscript_code: the folder holding .manuscript_code_root (README.md)
  p <- normalizePath(getwd(), mustWork = FALSE)
  while (!file.exists(file.path(p, ".manuscript_code_root"))) {
    if (dirname(p) == p) stop("run from inside the manuscript_code folder (see README.md)")
    p <- dirname(p)
  }
  p
}), "2020_zemb_microOpen_spike")

## Kept as a named constant for documentation and for anything that wants to ask
## "which tubes are controls?" without re-deriving it. The DROP is by mgfeces,
## not by this list -- see the header.
ZEMB_CONTROL_TUBES <- c("oz1802-16", "oz1802-17", "oz1802-18",   # no internal standard
                        "oz1802-19", "oz1802-20")                # no feces (19 spike-only, 20 standard-only)

parse_2020_zemb_microbiologyOpen_spike <- function(dir = NULL, drop_controls = TRUE) {

  ## `paths`/`local` in the old copies; nothing ever passed them. An explicit
  ## `dir` still wins, so a caller with the data somewhere else is not stuck.
  if (is.null(dir)) {
    dir <- ZEMB_DIR_CANDIDATES[dir.exists(ZEMB_DIR_CANDIDATES)][1]
    if (is.na(dir))
      stop("none of the candidate Zemb directories exist:\n  ",
           paste(ZEMB_DIR_CANDIDATES, collapse = "\n  "),
           "\npass dir = \"...\" explicitly", call. = FALSE)
  }
  ## No setwd(). The old copies changed the working directory as a side effect of
  ## parsing; nothing depended on it -- the two scripts that read a relative path
  ## later set their own directory first -- and a parser silently moving the
  ## process is how a later write lands in the wrong study's folder.

  metadata <- read.csv(file.path(dir, "zemb_metadata.csv"), row.names = 1)
  counts   <- read.csv(file.path(dir, "zemb_counts.csv"),   row.names = 1)
  scale    <- read.csv(file.path(dir, "zemb_qPCR.csv"),     row.names = 1)

  colnames(counts) <- paste0("oz1802-",
    unlist(lapply(str_split(colnames(counts), "Tube"), function(x) x[2])))
  counts   <- counts[, colnames(counts) %in% rownames(metadata)]
  metadata <- metadata[rownames(metadata) %in% colnames(counts), ]
  scale    <- scale[rownames(scale) %in% colnames(counts), ]

  ## The old copies called all.equal() and threw the answer away, so a
  ## misalignment printed nothing and carried on. These stop.
  stopifnot(identical(rownames(metadata), rownames(scale)),
            identical(colnames(counts),   rownames(scale)))

  ## ---------------------------------------------------------------------------
  ## total16S_per_mg -- the DENSITY, and the column every analysis should use.
  ##
  ## average16S(copies_in_tube) is the total 16S in the tube, which is the density
  ## MULTIPLIED BY however much feces was weighed in. Sequencing cannot see that
  ## mass: the composition is w_dn = d_dn*mg_n / (mg_n * sum_d d_dn) = d_dn/sum_d
  ## d_dn, so the mass cancels. Any quantity derived from the composition -- a
  ## per-sample scale model, or theta^perp -- therefore has to be paired with the
  ## density, or the pipetting is read as biology.
  ##
  ## On this dataset that is not a rounding error. The two arms of the
  ## added_Coli > 0 contrast come from two different experiments in the source
  ## paper and were weighed out at 11-14 mg versus 38-55 mg, so
  ##
  ##   log2 ratio of copies_in_tube  =  -0.437
  ##   log2 ratio of mg of feces     =  -1.333   <- pipetting alone
  ##   log2 ratio of copies per mg   =  +0.896
  ##
  ## exactly, since tube = density x mg. The per-tube figure has the wrong SIGN.
  ##
  ## Computed here rather than read from the sheet's own average16S(copies/mg)
  ## column: that column is stored as text because the 0 mg tubes hold #DIV/0!,
  ## so reading it gives a character vector that silently poisons any arithmetic
  ## downstream. The two agree to 3 significant figures where both are defined.
  ##
  ## NA, not Inf, where there is no feces -- a tube with no mass has no density,
  ## and Inf would propagate into a log2() and out the other side as a number.
  ## ---------------------------------------------------------------------------
  mg <- suppressWarnings(as.numeric(as.character(metadata$mgfeces)))
  scale$total16S_per_mg <- ifelse(!is.na(mg) & mg > 0,
                                  scale$average16S.copies_in_tube. / mg,
                                  NA_real_)

  ## Drop the no-feces controls -- all three tables together, so they cannot come
  ## apart. mgfeces is numeric in this file; as.numeric() guards a future one that
  ## reads it as text and would otherwise make every comparison NA and drop
  ## nothing silently.
  dropped <- character(0)
  if (isTRUE(drop_controls)) {
    mg   <- suppressWarnings(as.numeric(as.character(metadata$mgfeces)))
    keep <- !is.na(mg) & mg > 0
    dropped <- rownames(metadata)[!keep]
    if (length(dropped)) {
      message("  [Zemb] dropping no-feces control tube(s): ",
              paste(dropped, collapse = ", "),
              "  (", nrow(metadata), " -> ", sum(keep), " samples)")
      metadata <- metadata[keep, , drop = FALSE]
      scale    <- scale[keep, , drop = FALSE]
      counts   <- counts[, keep, drop = FALSE]
    }
  }

  out <- list(scale = scale, metadata = metadata, counts = counts)
  ## Travels with the result so a downstream table can say what was excluded
  ## rather than the exclusion being invisible once the object is saved.
  attr(out, "dropped_controls") <- dropped
  attr(out, "zemb_dir")         <- dir
  out
}
