options(
  BioC_mirror = "https://bioconductor.org",
  repos = BiocManager::repositories()
)

library <- Sys.getenv(
  "COMMAKIT_R_LIBRARY",
  "/usr/local/lib/R/site-library"
)
lockfile <- "/opt/commakit/project/renv.lock"
description <- "/opt/commakit/project/DESCRIPTION"

dir.create(library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library, .libPaths()))

renv::restore(
  project = "/opt/commakit/project",
  library = library,
  lockfile = lockfile,
  prompt = FALSE,
  clean = FALSE
)

dcf <- read.dcf(description)
fields <- intersect(c("Depends", "Imports", "Suggests"), colnames(dcf))
declared <- paste(dcf[1L, fields], collapse = ",")
declared <- trimws(unlist(strsplit(gsub("\n", "", declared), ",")))
declared <- sub(" \\(.*", "", declared)
declared <- setdiff(
  declared,
  c("", "R", "grDevices", "methods", "stats", "utils")
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

required <- unique(c(declared, development))
missing <- required[
  !vapply(required, requireNamespace, logical(1L), quietly = TRUE)
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
