#!/usr/bin/env Rscript
# Checks that this machine can run the analyses. It does not run any analysis
# or install anything.
#
#   Rscript run/check_environment.R                  # parse scripts, check packages
#   Rscript run/check_environment.R --syntax-only    # only parse scripts
#
# Exits with status 1 if a script fails to parse or something is missing.

args <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(args, "--syntax-only")) > 0) {
  stop("Usage: Rscript run/check_environment.R [--syntax-only]")
}
if (!file.exists(".manuscript_code_root")) {
  stop("Run this from the repository root.")
}

# ---- 1. Check that every R script parses -------------------------------------
r_files <- list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
parse_failures <- character()
for (f in r_files) {
  tryCatch(parse(f), error = function(e) {
    parse_failures <<- c(parse_failures, f)
    message(f, ": ", conditionMessage(e))
  })
}
cat(length(r_files), "R scripts parsed;", length(parse_failures), "failed\n")

if ("--syntax-only" %in% args) {
  quit(status = if (length(parse_failures) > 0) 1 else 0)
}

# ---- 2. Compare installed R packages with the versions we used ---------------
cat("\nThis R:", R.version.string, "\n")
cat("Original runs used R 4.4.3 and Julia 1.10\n\n")

deps <- read.delim("dependencies.tsv", colClasses = "character")
installed <- installed.packages()
deps$installed_version <- installed[match(deps$package, rownames(installed)), "Version"]
deps$status <- ifelse(is.na(deps$installed_version), "MISSING",
               ifelse(deps$installed_version == deps$recorded_version,
                      "same", "different version"))
print(deps, row.names = FALSE)

# ---- 3. Check for the programs the run scripts call --------------------------
programs <- Sys.which(c("bash", "zsh", "julia"))
cat("\nPrograms (blank = not found):\n")
print(programs)

missing_something <- any(deps$status == "MISSING") || any(programs == "")
cat("\nNote: this does not check Julia packages, and different package versions",
    "can change results.\n")
quit(status = if (length(parse_failures) > 0 || missing_something) 1 else 0)
