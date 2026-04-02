---
paths:
  - "R/**/*.R"
  - "vignettes/**/*.Rmd"
  - "man/**/*.Rd"
---

# Documentation Standards

## Roxygen2

All documentation is generated via roxygen2. Run `devtools::document()` to rebuild `man/` and `NAMESPACE` after any roxygen2 change.

Every exported function requires:

```r
#' Short title (one line)
#'
#' Longer description explaining what the function does and why you'd use it.
#'
#' @param object A \code{commaData} object.
#' @param mod_type Character string specifying the modification type
#'   (e.g., \code{"6mA"}, \code{"5mC"}). If \code{NULL}, uses all
#'   types present in \code{object}.
#' @param mod_context Character string specifying the modification context
#'   (e.g., \code{"6mA_GATC"}). If \code{NULL}, uses all contexts present.
#' @return Description of what is returned, including structure and key columns.
#' @examples
#' data(comma_example_data)
#' result <- functionName(comma_example_data)
#' @export
```

Stub documentation (`"A dataframe."`, `"A string."`) is **not acceptable** for `@param` or `@return`.

Use `\donttest{}` for examples that require optional dependencies (e.g., clusterProfiler, methylKit) or are slow.

## Package-level Documentation

`R/comma-package.R` provides the `?comma` documentation page (Bioconductor requirement). It describes the five-step workflow: Load → QC → Annotate → Visualize → Differential methylation. Do not remove this file.

## Vignettes

Two vignettes in `vignettes/`:

- **`getting-started.Rmd`** (~214 lines) — end-to-end workflow using `comma_example_data`: construct → characterize → diff methylation → visualize
- **`multiple-modification-types.Rmd`** (~173 lines) — joint 6mA + 5mC analysis; demonstrates subsetting by `mod_type`, comparing patterns, visualizing both simultaneously

Vignettes use `BiocStyle` for formatting. Both must knit without errors before Bioconductor submission.
