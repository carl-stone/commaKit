# Check the current source checkout, not the DESCRIPTION used when the image
# was built. Version-only metadata edits need no rebuild, but new package
# dependencies must be available in the locked development image.
description <- read.dcf("/workspace/commaKit/DESCRIPTION")
fields <- intersect(c("Depends", "Imports", "Suggests"), colnames(description))
declared <- paste(description[1L, fields], collapse = ",")
declared <- trimws(unlist(strsplit(gsub("\n", "", declared), ",")))
declared <- sub(" \\(.*", "", declared)
declared <- setdiff(
  declared,
  c("", "R", "grDevices", "methods", "stats", "utils")
)

# A fast presence check, not an attempt to load every optional package on
# every command. Built-package checks still validate version requirements and
# whether those packages actually load.
available <- vapply(
  declared,
  function(package) length(find.package(package, quiet = TRUE)) > 0L,
  logical(1L)
)
missing <- declared[!available]
if (length(missing) > 0L) {
  stop(
    "Development image is missing package(s) declared in DESCRIPTION: ",
    paste(missing, collapse = ", "),
    ". Update renv.lock and rebuild the image.",
    call. = FALSE
  )
}
