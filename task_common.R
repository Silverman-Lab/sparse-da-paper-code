################################################################################
## task_common.R -- the plumbing every reporting task shares.
##
## The manuscript's tables and figures are produced by six task scripts:
##
##   Task 1 dataset summary table.R            T1  characteristics of each dataset
##   Task 2 DA results table.R                 T2  TP/FP, Sparse SSRV under the
##                                                 hybrid rule at the default
##                                                 bandwidth, against every other
##                                                 method
##   Task 3 sparse SSRV variant tables.R       S1-S3  all NINE Sparse SSRV
##                                                 variants (3 scale-model rules
##                                                 x 3 KDE bandwidths)
##   Task 4 stammler all pairs figure.R        F4  Stammler pool pairs, one panel
##                                                 per method, Sparse SSRV shown
##                                                 once (hybrid / nrd0)
##   Task 5 stammler all pairs variant figure.R  F5  the same fifteen pairs, one
##                                                 panel per Sparse SSRV variant
##   Task 6 LFC distribution diagnostic.R      F6  where the LFC mode sits
##
## Each of those runs on its own (`Rscript "Task 3 sparse SSRV variant tables.R"`)
## and starts by sourcing this file. Everything that more than one of them needs
## -- paths, the dataset registry, the truth sets, the LaTeX helpers, and the
## reader that pulls all nine Sparse SSRV variants off a saved fit -- lives here
## and nowhere else, so the main table and its supplementary variant cannot come
## to disagree about which run they describe.
##
## NOTHING IS REFITTED anywhere in the task scripts. They read the artifacts the
## DA scripts and `real data results.R` write. If a study has been re-run since
## `real data results.R` last ran, the checks below say so rather than quietly
## reporting the older fit.
################################################################################

################################################################################
## PATHS -- nothing below this block hardcodes an absolute path.
##
## ROOT is auto-detected rather than toggled by hand, so the same files run
## unedited on the cluster and on the laptop. The first candidate that actually
## exists wins; to point somewhere else, assign ROOT yourself before sourcing.
################################################################################
if (!exists("ROOT") || is.null(ROOT) || is.na(ROOT) || !dir.exists(ROOT)) {
  ROOT_CANDIDATES <- local({  # manuscript_code: the folder holding .manuscript_code_root (README.md)
    p <- normalizePath(getwd(), mustWork = FALSE)
    while (!file.exists(file.path(p, ".manuscript_code_root"))) {
      if (dirname(p) == p) stop("run from inside the manuscript_code folder (see README.md)")
      p <- dirname(p)
    }
    p
  })
  ROOT <- ROOT_CANDIDATES[dir.exists(ROOT_CANDIDATES)][1]
  if (is.na(ROOT)) {
    stop("None of the candidate project roots exist:\n  ",
         paste(ROOT_CANDIDATES, collapse = "\n  "))
  }
}
message("ROOT = ", ROOT)

DIR_SELEX      <- file.path(ROOT, "selex")
DIR_ZEMB       <- file.path(ROOT, "2020_zemb_microOpen_spike")
DIR_STAMMLER   <- file.path(ROOT, "2016_stammler_microbiome_micehuman")
DIR_JIN_       <- file.path(ROOT, "2022_Jin_natureComm_technicalReplicates")
DIR_PRASSE     <- file.path(ROOT, "2015_prasse_appliedEnvMicrobiology_wetlands")
DIR_VANDEPUTTE <- file.path(ROOT, "Vandeputte")

DIR_TAB <- file.path(ROOT, "manuscript_tables")
DIR_FIG <- file.path(ROOT, "manuscript_figures")
dir.create(DIR_TAB, showWarnings = FALSE)
dir.create(DIR_FIG, showWarnings = FALSE)

## The stamp the current round of fits carries, used to name outputs. Same value
## the Jin file and `real data results.R` use; RESULT_STAMPS below is what
## actually resolves an INPUT file, so this is a label, not a lookup key.
if (!exists("DATESTAMP")) DATESTAMP <- "09012026"

## read_arm_sig() / arm_tpfp() / collect_scale_records(): the per-arm records
## every DA script writes through record_scale_model(). The variant reader below
## uses read_arm_sig() for the default bandwidth and the saved sensitivity fits
## for the other two.
source(file.path(ROOT, "scale_model_record.R"))

################################################################################
## DATESTAMP RESOLUTION
##
## Every study writes its outputs under the datestamp of the run that produced
## them. Naming one stamp in the reader is how a table goes on reporting a March
## fit after the study was re-run in August -- silently, because the March files
## still exist. stamped() takes the newest stamp that actually exists and records
## which it used; report_stamps() prints the result.
##
## Lifted from `real data results.R`, which owns the same problem. Kept as a
## separate copy rather than sourced, because that file also loads every study's
## results and refreshes the cross-dataset artifacts -- far more than a table
## builder should do on import.
################################################################################
## 09012026 is the rowMeans(raw counts) >= 5 re-run (see the filter header on any
## of the DA scripts); 08242026 is the superseded August prevalence-filter round,
## kept resolvable so a partially re-run project still produces a table AND says
## so through report_stamps() below. 08262026 is here for Prasse alone: its
## August rawVSqseq script was renamed from _08262026_ to _08242026_ AFTER its
## last run, so those outputs still sit under that stamp.
##
## NOTE the interaction with the Prasse pattern further down: it carries a
## _mean5 suffix now, and no _mean5 file exists at any older stamp, so a stale
## Prasse cannot resolve silently -- stamped() runs out of candidates and stops.
## For every other study the pattern is stamp-only and the fallback is live.
RESULT_STAMPS <- c("09012026", "08242026", "08262026", "03212026")
.task_stamp_log <- new.env(parent = emptyenv())

## `pattern` carries one %s where the datestamp goes. `required = FALSE` returns
## NA instead of stopping, for a file a study legitimately may not have.
stamped <- function(dir, pattern, what = pattern, required = TRUE) {
  for (st in RESULT_STAMPS) {
    f <- file.path(dir, sprintf(pattern, st))
    if (file.exists(f)) {
      assign(paste(basename(dir), what), st, envir = .task_stamp_log)
      return(f)
    }
  }
  if (!required) return(NA_character_)
  stop("no file for ", what, " in ", dir, " under any known datestamp (",
       paste(RESULT_STAMPS, collapse = ", "), ")")
}

report_stamps <- function() {
  k <- ls(.task_stamp_log)
  if (!length(k)) return(invisible(NULL))
  d <- data.frame(file = k, stamp = unlist(mget(k, envir = .task_stamp_log)),
                  row.names = NULL, stringsAsFactors = FALSE)
  d <- d[order(d$stamp, d$file), ]
  cat("\n== datestamp actually read, per file ==\n")
  print(d, row.names = FALSE)
  if (length(unique(d$stamp)) > 1) {
    old <- setdiff(unique(d$stamp), RESULT_STAMPS[1])
    warning("this run mixes stamps: ", sum(d$stamp != RESULT_STAMPS[1]), " of ",
            nrow(d), " files are still at ", paste(old, collapse = "/"),
            " because no newer version exists. Re-run those studies to make the ",
            "comparison uniform.", call. = FALSE, immediate. = TRUE)
  }
  invisible(d)
}

################################################################################
## CROSS-DATASET ARTIFACTS, written by `real data results.R`
################################################################################
F_SCALE   <- file.path(ROOT, "real_data_scale_model.RData")
F_RESULTS <- file.path(ROOT, "real_data_results.RData")
F_ARM     <- file.path(ROOT, "real_data_arm_features.RData")
F_TRUTH   <- file.path(ROOT, "real_data_truth_sets.RData")

for (f in c(F_SCALE, F_RESULTS, F_ARM, F_TRUTH)) {
  if (!file.exists(f)) {
    stop("missing artifact:\n  ", f,
         "\nRun `real data results.R` to produce it, then re-run this task.")
  }
}

scale_records <- readRDS(F_SCALE)
fptp          <- readRDS(F_RESULTS)   # datasets x methods, cells "FP/(TP)"
arm_score     <- readRDS(F_ARM)$score
truth_sets    <- readRDS(F_TRUTH)     # keyed "<dataset>|<pair>", pair "NA" if none

## The saved FP/TP matrix is only as current as the last run of
## `real data results.R`. Nothing here can detect a stale FIT, but a results file
## older than the newest fit output is a clear signal. Compared against the DA
## scripts' own outputs, not against the bookkeeping files the reporting scripts
## rewrite, which would make this warn forever.
check_results_freshness <- function() {
  fits <- c(Sys.glob(file.path(ROOT, "*", "sparse_sens_*.RData")),
            Sys.glob(file.path(ROOT, "*", "res_tpfp_*.RData")),
            Sys.glob(file.path(ROOT, "*", "jin_results_all_levels_*.RData")))
  if (!length(fits)) return(invisible(NULL))
  t_res <- file.mtime(F_RESULTS); t_rec <- max(file.mtime(fits), na.rm = TRUE)
  if (is.na(t_res) || !is.finite(t_rec)) return(invisible(NULL))
  if (t_res < t_rec) {
    warning("real_data_results.RData (", format(t_res, "%Y-%m-%d %H:%M"),
            ") predates the newest fit output (", format(t_rec, "%Y-%m-%d %H:%M"),
            "). The results tables may describe an earlier fit; re-run ",
            "`real data results.R`.", call. = FALSE, immediate. = TRUE)
  } else {
    message("real_data_results.RData is newer than every fit output it summarises")
  }
  invisible(NULL)
}

################################################################################
## DATASET METADATA
##
## The editorial columns of T1 -- scenario, sample material, citation key,
## sequencing method, how truth was defined -- are properties of the studies, not
## of any saved object, so they are declared here.
##
## Truth-set SIZES are not declared here. They are read from
## real_data_truth_sets.RData, which `real data results.R` writes from the same
## objects it scores against, keyed "<dataset>|<pair>". Transcribing them is how
## the signal-sparsity column goes stale when a study is re-filtered.
##
## `pair` names the run each dataset's row reports, and is the key into
## real_data_scale_model.RData. It must match the run `real data results.R`
## scored, or T1's dimensions will describe a different fit from T2's counts.
##
## `sens` is the sensitivity-fit file for that run, with %s where the datestamp
## goes. Spelled out per analysis rather than derived: Prasse's two filters put
## the run name AFTER the stamp while Jin puts the subset BEFORE it, and a
## derived pattern that gets that wrong reads another run's file without failing.
################################################################################
DATASET_META <- list(
  list(key = "McMurrough2014", short = "selex", pair = NA_character_,
       dir = DIR_SELEX, sens = "sparse_sens_%s.RData",
       scenario = "Sparse Low Variance", label = "McMurrough et al (2014)",
       sample = "", cite = "mcmurrough2014selex", seq = "16s rRNA",
       truth_def = "Known Truth"),
  list(key = "Zemb2020", short = "zemb", pair = NA_character_,
       dir = DIR_ZEMB, sens = "sparse_sens_%s.RData",
       scenario = "Sparse Low Variance", label = "Zemb et al (2020)",
       sample = "Pig (fecal)", cite = "zemb2020absolute", seq = "16s rRNA",
       truth_def = "Known Truth"),
  list(key = "Stammler2016", short = "stammler", pair = "Pool6",
       dir = DIR_STAMMLER, sens = "sparse_sens_%s.RData",
       scenario = "Sparse Low Variance", label = "St\\\"ammler et al (2016)",
       sample = "Mice (cecum content)", cite = "stammler2016adjusting",
       seq = "16s rRNA", truth_def = "Known Truth"),
  ## ONE row, the POOLED cOTU run. Pooling the two cecal locations is what gives
  ## the measured-scale reference enough power to call anything -- 1 cOTU pooled
  ## against 0 at either location on its own -- so the pooled run is the one that
  ## can be scored, and it is what `real data results.R` reports (JIN_PAIR there).
  ## The sample material is therefore the whole cecum content, not a location.
  ##
  ## The stratified locations and the ASV / OTU representations are not dropped;
  ## they are the subject of the Jin subset tables, which read jin_subsets()
  ## rather than this entry.
  list(key = "Jin2022", short = "jin", pair = "both",
       dir = DIR_JIN_, sens = NA_character_,
       scenario = "Sparse Low Variance", label = "Jin et al (2022)",
       sample = "Mice (cecum content)", cite = "jin2022quantsinglecell",
       ## Plain "16s rRNA", matching every other row. The assay detail
       ## (single-cell, BarBIQ) belongs in the text, not in a column whose other
       ## six cells carry only the sequencing target -- it was the one entry
       ## forcing that column wide.
       seq = "16s rRNA", truth_def = "Scale-informed Truth"),
  ## `pair` is the key the DA script writes into scale_model_record.rds /
  ## sparse_arm_sig.rds: "<filter>_<input>", so it moved with the filter. It must
  ## agree with PRASSE_RECORD_FILTER in `real data results.R`, which is what
  ## keys the truth set these tables score against.
  list(key = "Prasse2015", short = "prasse", pair = "mean5_QSeq",
       dir = DIR_PRASSE, sens = "sparse_sens_%s_mean5_QSeq.RData",
       scenario = "Sparse High Variance", label = "Prasse et al (2015)",
       sample = "Soil", cite = "prasse2015soil", seq = "16s rRNA",
       truth_def = "Scale-informed Truth"),
  list(key = "Vandeputte2017", short = "vandeputte", pair = NA_character_,
       dir = DIR_VANDEPUTTE, sens = "sparse_sens_%s.RData",
       scenario = "Not Sparse High Variance", label = "Vandeputte et al (2017)",
       sample = "Human (fecal)", cite = "vandeputte2017quantitative",
       seq = "16s rRNA", truth_def = "Scale-informed Truth")
)
names(DATASET_META) <- vapply(DATASET_META, `[[`, character(1), "key")

## Key into truth_sets / arm records for one (dataset, pair).
meta_key <- function(dataset, pair = NA_character_) {
  paste0(dataset, "|", if (is.na(pair) || !nzchar(pair)) "NA" else as.character(pair))
}

## Truth set of one analysis; character(0) when the study has none this table can
## score against (Jin at every subset, Stammler's fifteen pool pairs).
truth_for <- function(dataset, pair = NA_character_) {
  k <- meta_key(dataset, pair)
  if (!(k %in% names(truth_sets))) return(character(0))
  as.character(truth_sets[[k]])
}

## The truth-set size of the run a dataset's row reports.
n_truth_of <- function(key) {
  m <- DATASET_META[[key]]
  tr <- truth_for(key, m$pair)
  if (!length(tr)) NA_integer_ else length(tr)
}

## Row lookup: dimensions of a given (dataset, pair).
dims_for <- function(key, pair) {
  r <- scale_records[scale_records$dataset == key &
                     ((is.na(pair) & is.na(scale_records$pair)) |
                      (!is.na(pair) & !is.na(scale_records$pair) &
                       scale_records$pair == pair)), ]
  if (!nrow(r)) {
    stop("no scale-model record for ", key,
         if (is.na(pair)) " (pair NA)" else paste0(" (pair '", pair, "')"))
  }
  r[1, ]
}

## Row lookup: dimensions of the run a dataset's row reports.
dims_of <- function(key) {
  m <- DATASET_META[[key]]
  r <- scale_records[scale_records$dataset == key &
                     ((is.na(m$pair) & is.na(scale_records$pair)) |
                      (!is.na(m$pair) & !is.na(scale_records$pair) &
                       scale_records$pair == m$pair)), ]
  if (!nrow(r)) {
    stop("no scale-model record for ", key,
         if (is.na(m$pair)) " (pair NA)" else paste0(" (pair '", m$pair, "')"),
         " -- the DATASET_META `pair` and the recorded run disagree")
  }
  r[1, ]
}

################################################################################
## JIN'S EMPTY GOLD REFERENCE
##
## The ALDEx2 fit on Jin's measured absolute abundance -- the reference every
## TP/FP figure for that study would rest on -- is EMPTY at both cecal locations,
## and not marginally: refitting gives a minimum BH-adjusted p of 0.179 proximally
## and 0.355 distally. At 4 versus 4 samples the correction over ~150 cOTUs
## leaves it without power. Consequences, all deliberate:
##
##   * T1 dashes the signal sparsity for both Jin rows rather than printing 0/150
##     and 0/153, which would read as "no features change".
##   * T2 and S1 omit Jin: scored against an empty reference every method takes
##     TP = 0 and every call becomes a false positive, which measures the
##     reference's power and not the methods.
##   * S3 reports Jin's nine variants as DISCOVERY COUNTS instead, which is the
##     only thing that subset can say without a truth set.
##
## check_jin_gold() re-reads the reference each run: if it ever becomes non-empty
## the choices above need revisiting, and a silent zero is what would hide that.
################################################################################
## prox/dist are the two STRATIFIED subsets. Since 2026-08-28 the DA script also
## fits the POOLED cOTU run, which `real data results.R` reports as the primary
## Jin result -- and pooling is exactly what gives the measured-scale reference
## enough power to call anything, so a check that looked only at prox and dist
## would keep reporting "empty" after the premise had stopped holding.
LOC_KEYS  <- c("prox", "dist")
JIN_RUNS  <- c("both", LOC_KEYS)

## The pooled run drops the location from its filename; the stratified ones carry
## it. Same rule the DA script writes them under.
jin_gold_set <- function(loc) {
  f <- file.path(DIR_JIN_, if (identical(loc, "both"))
                             paste0("tax_scale_", DATESTAMP, ".RData")
                           else paste0("tax_scale_", loc, "_", DATESTAMP, ".RData"))
  if (!file.exists(f)) return(NULL)
  g <- readRDS(f)
  if (is.null(g) || !nrow(g)) character(0) else as.character(g$category)
}

check_jin_gold <- function() {
  n <- vapply(JIN_RUNS, function(l) {
    g <- jin_gold_set(l); if (is.null(g)) NA_integer_ else length(g)
  }, integer(1))
  if (all(is.na(n))) {
    message("Jin gold reference: not found on disk -- treated as empty")
    return(invisible(NULL))
  }
  ## The arrangement this now guards is the reverse of the old one: Jin IS
  ## reported, on the strength of the pooled reference. So the failure to catch
  ## is that reference going empty -- at which point Jin's sparsity, its row in
  ## the results tables and its truth set all stop meaning anything.
  n_both <- n[["both"]]
  if (is.na(n_both) || n_both == 0) {
    warning("the Jin POOLED gold reference is empty (",
            paste(sprintf("%s: %s", JIN_RUNS, ifelse(is.na(n), "absent", n)), collapse = ", "),
            "). Jin's sparsity in Table 1 and its row in the scored results tables both ",
            "rest on it being non-empty; scored against an empty reference every method ",
            "takes TP = 0 and every call becomes a false positive. Withhold Jin again ",
            "(SPARSITY_WITHHELD / DROP_FROM_RESULTS) or find out why it emptied.",
            call. = FALSE, immediate. = TRUE)
  } else {
    message("Jin pooled gold reference calls ", n_both, " cOTU(s)",
            sprintf(" (stratified: %s)", paste(sprintf("%s %s", LOC_KEYS,
                    ifelse(is.na(n[LOC_KEYS]), "absent", n[LOC_KEYS])), collapse = ", ")),
            " -- Jin is reported on the pooled run")
  }
  invisible(n)
}

## Signal sparsity is withheld where a study's reference cannot support one.
## Jin came off this list once the POOLED run gave it a non-empty reference; the
## stratified locations still could not support one, which is why the pooled run
## is the row Table 1 shows.
SPARSITY_WITHHELD <- character(0)

## Datasets dropped from the SCORED results tables. Empty now: Jin's pooled run
## is scored like every other study. Its ASV and OTU representations still are
## not -- they carry no measured scale at all -- and appear only in the Jin
## subset tables.
DROP_FROM_RESULTS <- character(0)

## The methods T2 compares, in column order, and how they are labelled.
METHOD_LABEL <- c(
  sparse_SSRV   = "Sparse SSRV",
  sparse_SSRV_B = "Sparse SSRV (bootstrap)",
  sparse_SSRV_L = "Sparse SSRV (Laplace)",
  ## "ALDEx2_def", not "ALDEx2". This column is aldex(gamma = 0.5) -- ALDEx2 with
  ## its DEFAULT SCALE MODEL, not vanilla ALDEx2. Printed as bare "ALDEx2" it
  ## read as the scale-free run, which is a different fit that lives under the
  ## key ALDEx2_original and is not reported here.
  ##
  ## Kept as a PLAIN string: this vector also labels the Task 4 figure panels,
  ## where a LaTeX-escaped underscore would render literally. Task 2 escapes it
  ## on the way into the table instead.
  ALDEx2_def    = "ALDEx2_def",
  DESeq2        = "DESeq2",
  edgeR         = "edgeR",
  limma         = "limma",
  baySeq        = "baySeq",
  MiMIX         = "MiMIX",
  `ANCOM-BC2`   = "ANCOM-BC2",
  `ANCOM-BC`    = "ANCOM-BC",
  ANCOM         = "ANCOM",
  LinDA         = "LinDA",
  CKF           = "CKF"
)
METHODS_OTHER <- setdiff(names(METHOD_LABEL),
                         c("sparse_SSRV", "sparse_SSRV_B", "sparse_SSRV_L"))

################################################################################
## LATEX HELPERS
##
## Character-for-character the definitions in
## 2022_Jin.../jin_DA_biological_interpretation.R, which owns them for the Jin
## tables. Repeated here so a task script that has nothing to do with Jin does
## not have to source that file -- sourcing it re-runs the whole Jin analysis and
## rewrites its outputs, which is a surprising side effect of building Table 1.
################################################################################
DASH <- "--"

tex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  for (ch in c("&", "%", "$", "#", "_", "{", "}")) {
    x <- gsub(ch, paste0("\\", ch), x, fixed = TRUE)
  }
  x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("^", "\\textasciicircum{}", x, fixed = TRUE)
  x
}

row_tex <- function(...) paste0("        ", paste(c(...), collapse = " & "), " \\\\")

## A stack of entries in one cell, TOP-aligned: a nested tabular with [t] puts
## its first line on the outer row's baseline.
tex_block <- function(entries, sep = " \\\\[3pt] ", align = "l") {
  if (!length(entries)) return(DASH)
  paste0("\\begin{tabular}[t]{@{}", align, "@{}}",
         paste(entries, collapse = sep),
         "\\end{tabular}")
}

rot <- function(x) paste0("\\rotatebox{75}{", x, "}")

## A scenario spanning several dataset rows.
## TWO lines, not one per word. Splitting on every space set "Sparse Low
## Variance" as three lines and "Not Sparse High Variance" as four, which made
## the scenario column taller than the rows it spans. The qualifier is always the
## last two words ("Low Variance", "High Variance"), so those go on line two and
## whatever precedes them on line one -- "Sparse" / "Low Variance",
## "Not Sparse" / "High Variance".
## `block` as in stack_paren(): TRUE when the table's other cells are
## tabular[t] blocks. It matters only in the SINGLE-row case -- \multirow does
## its own vertical centring across the rows it spans, so the n > 1 branch is
## unaffected -- but that single-row case is exactly where a \shortstack sat a
## line above the figures beside it, leaving a gap over the one-dataset
## categories.
## `lines_per_row` is how many TEXT lines one table row occupies -- 2 in both
## T1 and T2, whose every cell is a two-line stack. \multirow centres its content
## over nrows x \baselineskip, so passing the number of TABLE rows centres it
## over half the block it actually spans and the label rides up: with four
## two-line datasets, \multirow{4} put "Sparse / Low Variance" beside the first
## and second rows instead of the middle. Passing n * lines_per_row centres it
## over the true height.
scenario_cell <- function(text, n, block = FALSE, lines_per_row = 2L) {
  w <- strsplit(text, " ")[[1]]
  parts <- if (length(w) <= 2) w
           else c(paste(utils::head(w, -2), collapse = " "),
                  paste(utils::tail(w, 2), collapse = " "))
  txt <- paste(parts, collapse = " \\\\ ")
  if (n == 1L) {
    return(if (block) tex_block(parts, sep = " \\\\ ", align = "c")
           else paste0("\\shortstack{", txt, "}"))
  }
  paste0("\\multirow{", n * lines_per_row, "}{*}{\\shortstack{", txt, "}}")
}

## A label with a trailing parenthetical, set on two lines: "McMurrough et al"
## over "(2014)", "Pig" over "(fecal)". Keeps the Dataset and Sample columns
## narrow enough that the table does not need \resizebox to fit the text block.
## Text with no trailing parenthetical -- "Soil" -- is returned on one line, and
## an empty string stays empty rather than becoming an empty \shortstack.
## `block = TRUE` emits a \begin{tabular}[t] block instead of a \shortstack.
## The two are NOT interchangeable in a row that mixes them: \shortstack centres
## its lines on the row baseline while tabular[t] hangs them from the top, so a
## two-line \shortstack sits a full line above a two-line tabular[t] beside it.
## T1's cells are all \shortstack and agree with each other; T2's numeric cells
## are tabular[t], so its text cells must be too or the dataset name floats above
## the figures it labels.
stack_paren <- function(text, block = FALSE) {
  if (!nzchar(text)) return("")
  m <- regmatches(text, regexec("^(.*[^ ]) +(\\([^()]*\\))$", text))[[1]]
  parts <- if (length(m) == 3L) c(m[2], m[3]) else text
  if (block) tex_block(parts, sep = " \\\\ ", align = "c")
  else paste0("\\shortstack{", paste(parts, collapse = " \\\\ "), "}")
}

## "FP/(TP)" as saved -- summarize_tpfp() stores TP in column 1 and FP in column
## 2, and `real data results.R` pastes them in the order FP then (TP) -- split
## back into its two numbers. Withheld cells are "--".
split_fptp <- function(cell) {
  if (is.na(cell) || cell %in% c("--", "NA/(NA)")) return(c(NA_character_, NA_character_))
  m <- regmatches(cell, regexec("^(.*)/\\((.*)\\)$", cell))[[1]]
  if (length(m) != 3) return(c(NA_character_, NA_character_))
  m[2:3]
}

## Small counts spelled out, for captions that name how many rows they describe.
## Written out so a caption can COUNT its own table rather than assert a number
## that goes stale the next time a study is added -- T5 said "the three studies"
## while listing four, for one release.
english_count <- function(n) {
  w <- c("zero", "one", "two", "three", "four", "five", "six", "seven",
         "eight", "nine", "ten")
  if (is.na(n) || n < 0 || n > 10) as.character(n) else w[n + 1L]
}

## One cell of a results table: false positives over true positives in
## parentheses.
## `bold` sets both numbers AND the parentheses in \textbf, used to mark the best
## false-to-true positive ratio in a row. The parentheses are inside the \textbf
## deliberately: a bold numeral between upright parentheses reads as a typo at
## table size.
fptp_cell <- function(cell, bold = FALSE) {
  v <- split_fptp(cell)
  if (is.na(v[1])) return(DASH)
  fp <- tex_escape(v[1]); tp <- tex_escape(v[2])
  if (isTRUE(bold))
    tex_block(paste0("\\textbf{", fp, "} \\\\ \\textbf{(", tp, ")}"), align = "c")
  else
    tex_block(paste0(fp, " \\\\ (", tp, ")"), align = "c")
}

## Which cells in one results row carry the best false-to-true positive ratio.
##
## Returns a logical over `cells`. The rule, and why it is not just min(FP/TP):
##
##   * TP = 0 is EXCLUDED, not scored as a tie at zero. A method that found none
##     of the truth set has no meaningful ratio, and 0/(0) would otherwise win
##     outright -- "found nothing" is not the best result in the row. Four of the
##     six rows contain such a cell.
##   * Every row's minimum is FP = 0, so the ratio alone cannot separate the
##     winners: on Stammler nine of thirteen methods reach it. Ties on the ratio
##     are therefore broken by the MOST true positives -- the most of the truth
##     set recovered at the same (zero) cost. That resolves five of the six rows
##     to a single cell.
##   * Whatever is still tied after that is all marked. Stammler stays a nine-way
##     tie because nine methods each recover both true positives with no false
##     ones; nothing in the data distinguishes them, so the table should not
##     pretend otherwise.
best_fptp <- function(cells) {
  v  <- vapply(cells, split_fptp, character(2))
  fp <- suppressWarnings(as.numeric(v[1, ]))
  tp <- suppressWarnings(as.numeric(v[2, ]))
  ratio <- ifelse(is.na(fp) | is.na(tp) | tp <= 0, NA_real_, fp / tp)
  if (all(is.na(ratio))) return(rep(FALSE, length(cells)))
  ok  <- !is.na(ratio) & ratio == min(ratio, na.rm = TRUE)
  ok & !is.na(tp) & tp == max(tp[ok], na.rm = TRUE)
}

## The same, from the two numbers rather than from a saved string.
tpfp_cell <- function(tp, fp) {
  if (is.na(tp) || is.na(fp)) return(DASH)
  tex_block(paste0(fp, " \\\\ (", tp, ")"), align = "c")
}

## TRUE when a task script is being sourced by another script (manuscript_tables.R)
## rather than run on its own. Task scripts still WRITE their outputs; they just
## stop echoing LaTeX and diagnostic grids over their caller's own output. Same
## idiom as the `jin.quiet` option the Jin file uses.
task_quiet <- function() isTRUE(getOption("task.quiet"))

## Write a table and say where it went.
write_table <- function(lines, name, dir = DIR_TAB) {
  f <- file.path(dir, paste0(name, ".tex"))
  writeLines(lines, f)
  message("Wrote ", f)
  invisible(f)
}

################################################################################
## CONSISTENCY CHECKS -- run by the task scripts that report the numbers checked
################################################################################

## No arm may report more true positives than the declared truth set holds.
check_truth_sizes <- function() {
  bad <- character(0)
  for (k in names(DATASET_META)) {
    n <- n_truth_of(k)
    if (is.na(n)) next
    s <- arm_score[arm_score$dataset == k, ]
    if (!nrow(s)) next
    mx <- max(c(s$TP_laplace, s$TP_bootstrap), na.rm = TRUE)
    if (mx > n) bad <- c(bad, sprintf("%s: %d TP recorded against a declared truth set of %d",
                                      k, mx, n))
  }
  if (length(bad)) {
    warning("declared truth-set sizes look wrong:\n  ", paste(bad, collapse = "\n  "),
            "\nRe-run `real data results.R` so the saved truth sets match the scores.",
            call. = FALSE, immediate. = TRUE)
  } else {
    message("Declared truth-set sizes are consistent with the recorded TP counts")
  }
  invisible(NULL)
}

## sparse_SSRV is by construction whichever arm the selection rule chose, so its
## cell must equal that arm's cell. When it does not, the method-level TP/FP
## object and the per-arm record were produced by different runs of that study's
## DA script, and the hybrid column describes a fit nothing else here describes.
check_hybrid_consistency <- function() {
  bad <- character(0); noted <- character(0)
  for (k in names(DATASET_META)) {
    m <- DATASET_META[[k]]
    if (!(m$short %in% rownames(fptp))) next
    d <- dims_of(k)
    arm_col <- if (identical(d$scale_model, "bootstrap")) "sparse_SSRV_B" else "sparse_SSRV_L"
    hyb <- fptp[m$short, "sparse_SSRV"]; arm <- fptp[m$short, arm_col]
    if (arm %in% c("--", "NA/(NA)")) next            # deliberately withheld
    if (identical(hyb, arm)) next
    ## Two different causes, and only one is fixable by re-running. The method row
    ## comes from that study's own res_tpfp, scored against the gold reference its
    ## DA script built; the arm rows come from arm_tpfp(), scored against the truth
    ## set `real data results.R` declares. When those two references differ, the
    ## rows are not comparable and no re-run will reconcile them.
    empty_ref <- identical(k, "Jin2022") && all(lengths(lapply(LOC_KEYS, jin_gold_set)) == 0)
    line <- sprintf("%s: hybrid %s but selected arm (%s) %s", k, hyb, d$scale_model, arm)
    if (k %in% DROP_FROM_RESULTS || empty_ref) {
      noted <- c(noted, paste0(line,
        if (empty_ref) "  [its gold reference is empty, so the two columns are scored against different references]" else "",
        if (k %in% DROP_FROM_RESULTS) "  [not reported in the scored tables]" else ""))
    } else {
      bad <- c(bad, line)
    }
  }
  if (length(noted)) {
    message("Hybrid/arm mismatch in datasets this table does not report:\n  ",
            paste(noted, collapse = "\n  "))
  }
  if (length(bad)) {
    warning("the hybrid Sparse SSRV column does not match the selected arm for:\n  ",
            paste(bad, collapse = "\n  "),
            "\nRe-run those studies' DA scripts, then `real data results.R`, so the ",
            "method-level and per-arm records come from one fit.",
            call. = FALSE, immediate. = TRUE)
  } else {
    message("Hybrid Sparse SSRV matches the selected arm in every reported dataset")
  }
  invisible(bad)
}

## No method may report more true positives than the truth set holds. This
## catches a truth set redefined without DATASET_META following it.
check_tp_bound <- function() {
  bad <- character(0)
  for (k in names(DATASET_META)) {
    m <- DATASET_META[[k]]
    n <- n_truth_of(k)
    if (is.na(n) || !(m$short %in% rownames(fptp))) next
    for (cc in colnames(fptp)) {
      tp <- suppressWarnings(as.numeric(split_fptp(fptp[m$short, cc])[2]))
      if (!is.na(tp) && tp > n) {
        bad <- c(bad, sprintf("%s / %s: %s true positives against a truth set of %d",
                              k, cc, format(tp), n))
      }
    }
  }
  if (length(bad)) {
    warning("more true positives than the truth set can contain:\n  ",
            paste(bad, collapse = "\n  "),
            "\nThat study's method-level TP/FP object was scored against a different, ",
            "larger gold reference than the truth set the Sparse SSRV arms are scored ",
            "against, so those columns are not comparable.",
            call. = FALSE, immediate. = TRUE)
  } else {
    message("No method reports more true positives than its truth set holds")
  }
  invisible(bad)
}

################################################################################
## THE NINE SPARSE SSRV VARIANTS
##
## Sparse SSRV has two knobs that the manuscript reports separately, and their
## nine combinations are what the supplementary tables and figures examine:
##
##   SCALE MODEL   how theta^perp's uncertainty is obtained
##                   bootstrap   resampling spread of the LFC mode
##                   laplace     closed-form curvature of the KDE at its mode
##                   hybrid      whichever of the two the selection rule picked
##                               -- so a hybrid cell ALWAYS duplicates one of the
##                               other two in the same bandwidth column
##
##   BANDWIDTH     the KDE bandwidth the mode and its curvature are read at
##                   nrd0 (default), nrd, SJ-ste
##
## sigma^2 for the laplace arm is O(h^2) -- exactly h^2 when every LFC coincides
## -- so the bandwidth IS the Laplace scale uncertainty, and the two knobs are not
## independent axes of nuisance: this grid is the sensitivity analysis for the
## estimated scale itself.
##
## WHERE THE NUMBERS COME FROM. Nothing here refits anything.
##
##   nrd0           <dir>/sparse_arm_sig.rds, written by record_scale_model()
##                  inside the DA script at the moment of the main fit. It stores
##                  the two arms' significant feature IDs and which arm was
##                  selected, so all three rules come off it directly.
##
##   nrd, SJ-ste    the saved sensitivity fit, <dir>/sparse_sens_*.RData, one
##                  element per bandwidth. check_sparse_bw_sens() re-attaches
##                  attr "arm_out" (both arms' full tables) and attr
##                  "scale_model" (the arm the rule picked AT THAT BANDWIDTH) to
##                  each element, so the same three rules come off those too.
##
## The hybrid set is derived the same way at every bandwidth -- take the selected
## arm's significant features -- rather than trusting the already-subset table the
## sensitivity file stores, and the two are compared: `out[out$sig.q == "TRUE", ]`
## keeps an all-NA row for any feature whose sig.q is NA, so that table can carry
## rows the arm table does not. A disagreement is reported, never silently taken.
################################################################################
SSRV_RULES <- c("hybrid", "bootstrap", "laplace")
SSRV_BWS   <- c("nrd0", "nrd", "SJ-ste")
SSRV_VARIANTS <- as.vector(t(outer(SSRV_RULES, SSRV_BWS, paste, sep = "|")))

## How the two settings are spelled in tables and on facet strips. The bandwidth
## map is an identity today; it is the one place to rename, e.g. to
## "nrd0 (default)", without hunting through three files.
SSRV_RULE_LABEL <- c(hybrid = "Hybrid", bootstrap = "Bootstrap", laplace = "Laplace")
SSRV_BW_LABEL   <- c(nrd0 = "nrd0", nrd = "nrd", `SJ-ste` = "SJ-ste")

## Feature IDs an arm table calls significant. Two absences that must not be
## confused, and the return value is what keeps them apart:
##
##   NULL           no table at all -- that arm was never recorded
##   character(0)   the arm was recorded and called nothing
##
## NA sig.q is "not called", never "called": an unusable arm must not manufacture
## discoveries.
.sig_ids <- function(tab) {
  if (is.null(tab)) return(NULL)
  if (!nrow(tab))   return(character(0))
  s <- tab$sig.q
  rownames(tab)[!is.na(s) & s]
}

## The same distinction for the bare id vectors in sparse_arm_sig.rds: a record
## written before an arm was stored has the field missing, and as.character(NULL)
## would turn that into "called nothing".
.arm_ids <- function(e, nm) if (is.null(e[[nm]])) NULL else as.character(e[[nm]])

## All nine variants for ONE analysis.
##
##   dir      the study directory holding sparse_arm_sig.rds
##   dataset  the name record_scale_model() recorded under
##   pair     the run key inside that record (NA when the study recorded one run)
##   sens     path to the sensitivity fit, or NA to skip the non-default
##            bandwidths
##
## Returns a named list "<rule>|<bandwidth>" -> character vector of called
## features. A variant with no record on disk is ABSENT from the list rather than
## character(0), so "not fitted" and "fitted, found nothing" stay distinguishable.
ssrv_variant_calls <- function(dir, dataset, pair = NA_character_, sens = NA_character_,
                               quiet = FALSE) {
  out <- list()
  ## Which arm the rule picked, per bandwidth. Returned as an attribute because a
  ## hybrid cell is only interpretable next to it: hybrid ALWAYS duplicates one of
  ## the other two arms, and which one it duplicates is the result.
  picked <- setNames(rep(NA_character_, length(SSRV_BWS)), SSRV_BWS)

  e <- suppressWarnings(read_arm_sig(dir, dataset, pair))
  if (is.null(e)) {
    if (!quiet)
      warning("no per-arm record for ", meta_key(dataset, pair), " in ", dir,
              " -- the default-bandwidth variants are withheld", call. = FALSE)
  } else {
    out[["bootstrap|nrd0"]] <- .arm_ids(e, "bootstrap")
    out[["laplace|nrd0"]]   <- .arm_ids(e, "laplace")
    if (!is.null(e$selected) && e$selected %in% names(e)) {
      picked[["nrd0"]] <- e$selected
      out[["hybrid|nrd0"]] <- .arm_ids(e, e$selected)
    }
  }

  if (!is.na(sens) && file.exists(sens)) {
    s <- readRDS(sens)
    for (bw in names(s)) {
      tab <- s[[bw]]
      ao  <- attr(tab, "arm_out")
      sm  <- attr(tab, "scale_model")
      if (is.null(ao)) {
        if (!quiet)
          warning(meta_key(dataset, pair), " at bandwidth '", bw,
                  "' carries no per-arm tables -- re-run that DA script so ",
                  "check_sparse_bw_sens() records them", call. = FALSE)
        next
      }
      out[[paste0("bootstrap|", bw)]] <- .sig_ids(ao$bootstrap)
      out[[paste0("laplace|", bw)]]   <- .sig_ids(ao$laplace)
      if (!is.null(sm) && sm %in% names(ao)) {
        picked[[bw]] <- sm
        hyb <- .sig_ids(ao[[sm]])
        out[[paste0("hybrid|", bw)]] <- hyb
        rn <- rownames(tab); rn <- rn[!is.na(rn)]
        if (!setequal(rn, hyb) && !quiet) {
          warning(meta_key(dataset, pair), " at bandwidth '", bw,
                  "': the stored hybrid table lists ", length(rn),
                  " features but the selected arm (", sm, ") calls ", length(hyb),
                  ". Reporting the arm.", call. = FALSE)
        }
      }
    }
  } else if (!is.na(sens) && !quiet) {
    warning("no sensitivity fit at ", sens, " -- only the default bandwidth is ",
            "available for ", meta_key(dataset, pair), call. = FALSE)
  }

  out <- out[intersect(SSRV_VARIANTS, names(out))]
  attr(out, "selected") <- picked
  out
}

## The same, scored, as a tidy frame: one row per (analysis, variant).
##
## `truth` may be character(0), in which case TP/FP are NA and only `n_called`
## carries information -- the Jin case, where the measured-scale reference is
## empty and scoring against it would measure that reference's power.
ssrv_variant_scores <- function(dir, dataset, pair = NA_character_,
                                sens = NA_character_, truth = character(0),
                                label = dataset, quiet = FALSE) {
  calls  <- ssrv_variant_calls(dir, dataset, pair, sens, quiet = quiet)
  picked <- attr(calls, "selected")
  scored <- length(truth) > 0
  do.call(rbind, lapply(SSRV_VARIANTS, function(v) {
    ids  <- calls[[v]]
    have <- !is.null(ids)
    z    <- strsplit(v, "|", fixed = TRUE)[[1]]
    data.frame(
      dataset  = dataset, pair = if (is.na(pair)) NA_character_ else pair,
      label    = label, variant = v, rule = z[1], bandwidth = z[2],
      selected = unname(picked[[z[2]]]),
      n_called = if (have) length(ids) else NA_integer_,
      TP = if (have && scored) sum(ids %in% truth)  else NA_integer_,
      FP = if (have && scored) sum(!ids %in% truth) else NA_integer_,
      scored   = scored,
      stringsAsFactors = FALSE)
  }))
}

################################################################################
## THE THREE VARIANT TABLES, AS TIDY FRAMES
##
## Built here rather than in the task scripts, because Task 3 (tables) and Task 5
## (the Stammler figure) have to report the same numbers, and the surest way for
## a supplementary table and a supplementary figure to disagree is for each to
## read the fits its own way.
##
## One row per (analysis, variant). `scored` says whether TP/FP mean anything:
## FALSE for Jin, whose measured-scale reference is empty, where only `n_called`
## is informative.
################################################################################

## The five datasets Table 2 scores, at the run each of them reports.
ssrv_variants_scored <- function() {
  keys <- setdiff(names(DATASET_META), DROP_FROM_RESULTS)
  do.call(rbind, lapply(keys, function(k) {
    m <- DATASET_META[[k]]
    ssrv_variant_scores(m$dir, k, m$pair,
                        sens  = jin_or_meta_sens(k, m$pair, m$sens),
                        truth = truth_for(k, m$pair),
                        label = m$label)
  }))
}

## Jin's six subsets -- three feature representations at two cecal locations --
## as DISCOVERY COUNTS. See check_jin_gold() for why they cannot be scored.
ssrv_variants_jin <- function() {
  do.call(rbind, lapply(jin_subsets(), function(a)
    cbind(ssrv_variant_scores(a$dir, a$key, a$pair, sens = a$sens,
                              truth = truth_for(a$key, a$pair), label = a$label),
          level = a$level, location = a$location, row.names = NULL)))
}

## Stammler's fifteen pool pairs, scored against the two spike-ins.
ssrv_variants_stammler_pairs <- function() {
  do.call(rbind, lapply(stammler_pair_analyses(), function(a)
    cbind(ssrv_variant_scores(a$dir, a$key, a$pair, sens = a$sens,
                              truth = STAMMLER_TRUTH, label = a$label),
          ref = a$ref, comp = a$comp, row.names = NULL)))
}

################################################################################
## THE ANALYSIS REGISTRY
##
## Every fit the variant tables and figures report, in one place: the five scored
## datasets, Jin's six unscored subsets, and Stammler's fifteen pool pairs.
##
## `sens` is resolved through stamped(), so a study re-run under a newer stamp is
## picked up without editing anything here.
################################################################################

## Jin names its sensitivity fit by subset ("sparse_sens_ASV_prox_%s.RData"),
## which the single `sens` pattern in DATASET_META cannot express; every other
## study has one pattern for its one reported run.
jin_or_meta_sens <- function(key, pair, pattern) {
  if (identical(key, "Jin2022")) {
    ## The pooled run drops the subset from its filename -- sparse_sens_<stamp> --
    ## the same asymmetry its gold reference has. Every stratified subset carries
    ## it: sparse_sens_ASV_prox_<stamp>, and so on.
    pat <- if (identical(pair, "both")) "sparse_sens_%s.RData"
           else paste0("sparse_sens_", pair, "_%s.RData")
    return(stamped(DIR_JIN_, pat, paste0("sparse_sens ", pair), required = FALSE))
  }
  if (is.na(pattern)) return(NA_character_)
  stamped(DATASET_META[[key]]$dir, pattern, "sparse_sens", required = FALSE)
}

## The five datasets T2 scores, plus Jin's two cOTU rows -- one entry per row of
## T1, in T1's order.
main_analyses <- function() {
  out <- list()
  for (k in names(DATASET_META)) {
    m <- DATASET_META[[k]]
    subs <- if (is.null(m$subsets)) list(list(pair = m$pair, sample = m$sample)) else m$subsets
    for (s in subs) {
      out[[length(out) + 1]] <- list(
        key = k, short = m$short, pair = s$pair, dir = m$dir,
        label = m$label, scenario = m$scenario,
        sens = jin_or_meta_sens(k, s$pair, m$sens))
    }
  }
  out
}

## Jin's six subsets: three feature representations at two cecal locations. The
## record key is the bare location for cOTU and "<level>_<location>" otherwise,
## which is how the DA script wrote it.
JIN_LEVELS <- c("cOTU", "ASV", "OTU")
jin_subsets <- function() {
  out <- list()
  for (lv in JIN_LEVELS) for (lo in LOC_KEYS) {
    pr <- if (identical(lv, "cOTU")) lo else paste0(lv, "_", lo)
    out[[length(out) + 1]] <- list(
      key = "Jin2022", pair = pr, level = lv, location = lo, dir = DIR_JIN_,
      label = paste0(lv, ", ", if (lo == "prox") "proximal" else "distal"),
      sens = stamped(DIR_JIN_, paste0("sparse_sens_", pr, "_%s.RData"),
                     paste0("sparse_sens ", pr), required = FALSE))
  }
  out
}

## Stammler's fifteen pool pairs. The truth set is the two spike-ins, taken from
## the Pool6 entry of the saved truth sets rather than retyped, so it cannot
## drift from the set `real data results.R` scores the reported pair against.
STAMMLER_POOLS <- paste0("Pool", 1:6)
STAMMLER_PAIRS <- combn(STAMMLER_POOLS, 2, simplify = FALSE)
STAMMLER_TRUTH <- truth_for("Stammler2016", "Pool6")
if (!length(STAMMLER_TRUTH)) {
  stop("no truth set recorded for Stammler2016|Pool6 -- re-run `real data results.R`")
}

## The pair loop of sparse_DA_*.R keys its records "<ref>_vs_<comp>" and names
## its files "<prefix>_<ref>_<comp>_<stamp>.RData".
stammler_pair_analyses <- function() {
  lapply(STAMMLER_PAIRS, function(p) list(
    key = "Stammler2016", pair = paste0(p[1], "_vs_", p[2]),
    ref = p[1], comp = p[2], dir = DIR_STAMMLER,
    label = paste0(p[1], " vs ", p[2]),
    sens = stamped(DIR_STAMMLER, paste0("sparse_sens_", p[1], "_", p[2], "_%s.RData"),
                   paste0("sparse_sens ", p[1], "_", p[2]), required = FALSE),
    tpfp = stamped(DIR_STAMMLER, paste0("res_tpfp_", p[1], "_", p[2], "_%s.RData"),
                   paste0("res_tpfp ", p[1], "_", p[2]), required = FALSE)))
}

message("task_common.R loaded")
