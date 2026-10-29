#!/usr/bin/env Rscript
# Read-only checks; does not source analysis files, install packages, or run fits.
args <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(args, "--syntax-only"))) stop("Usage: Rscript run/check_environment.R [--syntax-only]")
if (!file.exists(".manuscript_code_root")) stop("Run from the repository root")
files <- list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
failed <- character()
for (f in files) {
  tryCatch(parse(f), error = function(e) {
    failed <<- c(failed, f)
    message(f, ": ", conditionMessage(e))
  })
}
cat(length(files), "R scripts parsed; failures:", length(failed), "\n")
if ("--syntax-only" %in% args) quit(status = as.integer(length(failed) > 0))

cat(R.version.string, "\nRecorded analysis environment: R 4.4.3, Julia 1.10\n")
deps <- read.delim("dependencies.tsv", colClasses = "character")
installed <- installed.packages()
deps$installed_version <- installed[match(deps$package, rownames(installed)), "Version"]
deps$status <- ifelse(is.na(deps$installed_version), "MISSING",
                      ifelse(deps$installed_version == deps$recorded_version, "matches", "version differs"))
print(deps, row.names = FALSE)
commands <- Sys.which(c("bash", "zsh", "julia"))
cat("\nExecutables:\n")
print(commands)
missing <- anyNA(deps$installed_version) || any(!nzchar(commands))
cat("\nVersion differences can affect results. This check does not validate Julia packages or numerical reproduction.\n")
quit(status = as.integer(length(failed) > 0 || missing))
