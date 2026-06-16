options(BioC_mirror = "https://bioconductor.org")

if (requireNamespace("BiocManager", quietly = TRUE)) {
  options(repos = BiocManager::repositories())
} else {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

cat("commaKit devcontainer setup\n")
cat("R version: ", as.character(getRversion()), "\n", sep = "")

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.8.tar.gz",
    repos = NULL,
    type = "source"
  )
}

renv_library <- renv::paths$library(project = getwd())
lock_dirs <- Sys.glob(file.path(renv_library, "00LOCK*"))
if (length(lock_dirs) > 0) {
  cat("Removing stale renv lock directories:\n")
  cat(paste0("  - ", lock_dirs, "\n"), sep = "")
  unlink(lock_dirs, recursive = TRUE, force = TRUE)
}

renv::restore(prompt = FALSE)
renv::load(project = getwd())

editor_packages <- c("languageserver", "httpgd", "lintr", "styler")
missing_editor_packages <- editor_packages[
  !vapply(editor_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_editor_packages) > 0) {
  renv::install(missing_editor_packages, prompt = FALSE)
}

cat("Devcontainer R environment is ready.\n")
