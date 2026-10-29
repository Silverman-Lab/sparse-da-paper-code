################################################################################
## plot_conventions.R -- the method colours, line types and legend order shared
## by every Study A figure.
##
## ONE definition, sourced by plot_studyA_combined.R and plot_runtime.R. They were
## inline in plot_studyA_combined.R until 2026-08-29; a second figure needed them
## and copying was not an option, because copying is precisely how they went
## wrong before. The original
##   simulation/src/06072025/Results combined updated.R
## redefines method_colors FOUR times and line_types twice, and only the last of
## each produced the published figures -- so a reader who copied from the wrong
## block got a figure that looked plausible and was not the published one.
##
## Values are that file's FINAL definitions (method_colors line 72, line_types
## line 89, neworder line 122), unchanged.
##
## The one rename from the original: "CKF" is called "MRCKF" in this pipeline.
##
## ANY figure keyed by method sources this. If a method is added, it is added
## HERE and nowhere else.
################################################################################

grey_color <- "grey80"

## Colour carries the method FAMILY, not the method: ANCOM and ANCOM-BC share
## plum1, and everything without a scale model shares grey80. That is deliberate
## -- it is what makes the scale-aware methods legible as a group -- and it is
## why the legend has to key on colour AND line type together. Neither aesthetic
## identifies a method on its own.
method_colors <- c(
  "Sparse_SSRV" = "red",
  "ALDEx2_clr"  = "skyblue1",
  "ALDEx2_def"  = "skyblue3",
  "ANCOM"       = "plum1",
  "ANCOM-BC"    = "plum1",
  "ANCOM-BC2"   = "plum3",
  "LinDA"       = "olivedrab3",
  "MiMIX"       = "yellow2",
  "DESeq2"      = grey_color,
  "baySeq"      = grey_color,
  "edgeR"       = grey_color,
  "limma"       = grey_color,
  "MRCKF"       = grey_color)          # "CKF" in the original

line_types <- c(
  "ALDEx2_clr"  = "dashed",
  "ALDEx2_def"  = "dotdash",
  "ANCOM"       = "dotdash",
  "ANCOM-BC"    = "dotted",
  "ANCOM-BC2"   = "longdash",
  "DESeq2"      = "dashed",
  "LinDA"       = "dotdash",
  "MiMIX"       = "dashed",
  "MRCKF"       = "solid",             # "CKF" in the original
  "baySeq"      = "dotted",
  "edgeR"       = "dotdash",
  "limma"       = "longdash",
  "Sparse_SSRV" = "solid")

## Legend order. ANCOM and ANCOM-BC share a colour and ANCOM-BC/limma share a
## line type, so the pairing of the two aesthetics is what separates them --
## which is why the legend must key on both.
neworder <- c("Sparse_SSRV", "ALDEx2_clr", "ALDEx2_def", "ANCOM", "ANCOM-BC",
              "ANCOM-BC2", "LinDA", "MiMIX", "DESeq2", "baySeq", "edgeR",
              "limma", "MRCKF")

stopifnot(setequal(names(method_colors), neworder),
          setequal(names(line_types),    neworder))
