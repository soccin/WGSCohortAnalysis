#!/usr/bin/env Rscript
#
# wcaNewProject.R: scaffold a new analysis project from templates/project.
#
#   Rscript $WCA_HOME/bin/wcaNewProject.R <project_dir> [prefix]
#
# Copies the template, fills `toolkit:` and `prefix:` in 00.PARAMS.yml, and
# creates data/raw, cache and results. Refuses to overwrite an existing
# directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: wcaNewProject.R <project_dir> [prefix]")
}
target <- args[[1]]
prefix <- if (length(args) >= 2) args[[2]] else basename(target)

script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
toolkit <- normalizePath(file.path(dirname(script), ".."))
template <- file.path(toolkit, "templates", "project")

if (dir.exists(target)) stop("Directory exists: ", target)
dir.create(target, recursive = TRUE)
files <- list.files(template, all.files = TRUE, recursive = TRUE, no.. = TRUE)
for (f in files) {
  dir.create(dirname(file.path(target, f)), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(template, f), file.path(target, f))
}
dir.create(file.path(target, "data", "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(target, "cache"), showWarnings = FALSE)
dir.create(file.path(target, "results"), showWarnings = FALSE)

params_file <- file.path(target, "00.PARAMS.yml")
txt <- readLines(params_file)
txt <- sub("^toolkit: .*", sprintf("toolkit: %s", toolkit), txt)
txt <- sub("^prefix: .*", sprintf("prefix: %s", prefix), txt)
writeLines(txt, params_file)

cat(sprintf("Created %s\n  toolkit: %s\n  prefix:  %s\nNext: put the manifests in data/raw/, edit 00.PARAMS.yml, then Rscript run_all.R\n",
  normalizePath(target), toolkit, prefix))
