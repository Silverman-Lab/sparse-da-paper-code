library(readr)
parse_2016_stammler_microbiome_micehuman <- function(paths = NULL) {
  ## Was a hardcoded cluster path, so this parser could not run on a laptop at
  ## all. Same auto-detection every other script in the project uses; `paths`
  ## still wins when a caller supplies one.
  local <- if (!is.null(paths)) paths else {
    cand <- file.path(local({  # manuscript_code: the folder holding .manuscript_code_root (README.md)
      p <- normalizePath(getwd(), mustWork = FALSE)
      while (!file.exists(file.path(p, ".manuscript_code_root"))) {
        if (dirname(p) == p) stop("run from inside the manuscript_code folder (see README.md)")
        p <- dirname(p)
      }
      p
    }), "2016_stammler_microbiome_micehuman")
    hit <- cand[dir.exists(cand)][1]
    if (is.na(hit)) stop("no Stammler data directory found among:\n  ",
                         paste(cand, collapse = "\n  "))
    paste0(hit, "/")
  }
  
  dat <- read_table(paste0(local, "40168_2016_175_MOESM9_ESM.txt"), skip = 1)
  counts <- as.matrix(dat[,2:37])
  rownames(counts) <- dat$OTU_ID
  
  tax <- as.matrix(dat[,38:44])
  rownames(tax) <- dat$OTU_ID
  
  scale <- readxl::read_xlsx(paste0(local, "40168_2016_175_MOESM1_ESM.xlsx"), skip = 5, col_names = FALSE, sheet = "Tabelle1")[,c(1,5)]
  colnames(scale) <- c("Sample", "16S rDNA copies per sample")
  
  idmap <- read.csv(paste0(local, "40168_2016_175_MOESM11_ESM.csv"), sep = "\t")
  
  metadata <- as.data.frame(read_table(paste0(local, "40168_2016_175_MOESM7_ESM.txt"))[-c(23,24),])
  rownames(metadata) <- metadata$SampleID
  
  samID <- idmap$SampleID[match(scale$Sample, idmap$SampleNumber)]
  scale <- scale[samID %in% rownames(metadata),]
  scale$SampleID <- samID[samID %in% rownames(metadata)]
  
  counts <- counts[,match(rownames(metadata), colnames(counts))]
  
  #remove samples that doesn't match
  counts <- counts[,colnames(counts) %in% rownames(metadata)]
  metadata <- metadata[metadata$SampleID %in% colnames(counts),]
  rownames(metadata) <- metadata$SampleID
  scale <- scale[scale$SampleID %in% colnames(counts),]
  
  all.equal(colnames(counts), rownames(metadata))
  all.equal(colnames(counts), scale$SampleID)
  all.equal(rownames(metadata), scale$SampleID)
  
  return(list(scale=scale, metadata=metadata, counts=counts, tax=tax))
}

