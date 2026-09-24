options(
  BioC_mirror = "https://bioconductor.org",
  repos = BiocManager::repositories()
)

library <- Sys.getenv(
  "COMMAKIT_R_LIBRARY",
  "/usr/local/lib/R/site-library"
)
lockfile <- "/opt/commakit/project/renv.lock"

dir.create(library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library, .libPaths()))

renv::restore(
  project = "/opt/commakit/project",
  library = library,
  lockfile = lockfile,
  prompt = FALSE,
  clean = FALSE
)

development <- c(
  "BiocCheck",
  "devtools",
  "pkgdown",
  "rcmdcheck",
  "renv",
  "roxygen2",
  "sessioninfo",
  "V8"
)

missing <- development[
  !vapply(development, requireNamespace, logical(1L), quietly = TRUE)
]

if (length(missing) > 0L) {
  stop(
    "Development image is missing required package(s): ",
    paste(missing, collapse = ", ")
  )
}

cat("commaKit development environment ready\n")
cat("R ", as.character(getRversion()), "\n", sep = "")
cat("Bioconductor ", as.character(BiocManager::version()), "\n", sep = "")
