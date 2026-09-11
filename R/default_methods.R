# Default methods for commaKit S4 generics.
#
# Every package-defined generic has a commaData method that does the real
# work. Without a default method, calling the generic on any other class
# fails with the raw S4 dispatch error:
#
#   Error: unable to find an inherited method for function 'methylation'
#   for signature 'object = "data.frame"'
#
# That error names the function and the wrong class but never says what
# the caller should have passed. For the package's target audience —
# scientists and coding agents assembling analyses from generated code —
# the actionable contract is: "X() expects a commaData object". These
# default methods provide it.
#
# The default signature is "ANY": it fires only when no more specific
# method exists, so future methods for other classes (or the existing
# commaData methods) always take precedence and nothing is shadowed.

# ── Accessors (R/accessors.R) ────────────────────────────────────────────────

#' @rdname methylation
setMethod("methylation", "ANY", function(object) {
  stop("methylation() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname siteCoverage
setMethod("siteCoverage", "ANY", function(object) {
  stop("siteCoverage() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname modCounts
setMethod("modCounts", "ANY", function(object) {
  stop("modCounts() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname canonicalCounts
setMethod("canonicalCounts", "ANY", function(object) {
  stop(
    "canonicalCounts() expects a commaData object, got ",
    class(object)[1],
    "."
  )
})

#' @rdname otherModCounts
setMethod("otherModCounts", "ANY", function(object) {
  stop(
    "otherModCounts() expects a commaData object, got ",
    class(object)[1],
    "."
  )
})

#' @rdname assayProvenance
setMethod("assayProvenance", "ANY", function(object) {
  stop(
    "assayProvenance() expects a commaData object, got ",
    class(object)[1],
    "."
  )
})

#' @rdname sampleInfo
setMethod("sampleInfo", "ANY", function(object) {
  stop("sampleInfo() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname siteInfo
setMethod("siteInfo", "ANY", function(object) {
  stop("siteInfo() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname modTypes
setMethod("modTypes", "ANY", function(object) {
  stop("modTypes() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname motifs
setMethod("motifs", "ANY", function(object) {
  stop("motifs() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname modContexts
setMethod("modContexts", "ANY", function(object) {
  stop("modContexts() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname genomeSizes
setMethod("genomeSizes", "ANY", function(object) {
  stop("genomeSizes() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname motifSites
setMethod("motifSites", "ANY", function(object) {
  stop("motifSites() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname minCoverage
setMethod("minCoverage", "ANY", function(object) {
  stop("minCoverage() expects a commaData object, got ", class(object)[1], ".")
})

# ── Assay and result layers (R/assay_layers.R, R/result_layers.R,
#    R/results_methods.R) ────────────────────────────────────────────────────

#' @rdname assayLayers
setMethod("assayLayers", "ANY", function(object) {
  stop("assayLayers() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname resultLayers
setMethod("resultLayers", "ANY", function(object) {
  stop(
    "resultLayers() expects a commaData object, got ",
    class(object)[1],
    "."
  )
})

#' @rdname results
setMethod("results", "ANY", function(object, ...) {
  stop("results() expects a commaData object, got ", class(object)[1], ".")
})

#' @rdname filterResults
setMethod("filterResults", "ANY", function(object, ...) {
  stop(
    "filterResults() expects a commaData object, got ",
    class(object)[1],
    "."
  )
})
