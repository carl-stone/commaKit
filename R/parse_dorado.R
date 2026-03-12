NULL

#' Parse a Dorado BAM file with MM/ML tags (not yet implemented)
#'
#' This function will read a Dorado BAM file containing MM/ML base modification
#' tags and aggregate per-read modification probabilities into per-site beta
#' values. It is an internal function called when \code{caller = "dorado"}.
#'
#' @details
#' **Not yet implemented.** Dorado BAM parsing requires reading MM/ML tags
#' from aligned BAM files using \code{\link[Rsamtools]{BamFile}} and
#' \code{\link[Rsamtools]{scanBam}}. This is lower priority than modkit
#' parsing because the recommended workflow is to convert Dorado BAM output
#' to modkit pileup BED using the \code{modkit pileup} command before loading
#' data into \code{comma}.
#'
#' **Recommended workflow:**
#' \preformatted{
#'   modkit pileup aligned.bam output.bed --ref genome.fa
#' }
#' Then load \code{output.bed} using \code{caller = "modkit"}.
#'
#' @param file Character string. Path to the Dorado-aligned BAM file.
#' @param sample_name Character string. Sample name.
#' @param mod_type Character vector or \code{NULL}. Modification type filter.
#' @param min_coverage Integer. Minimum read depth. Default \code{5}.
#'
#' @return Does not return; stops with an informative error message.
#'
#' @keywords internal
.parseDorado <- function(file, sample_name, mod_type = NULL, min_coverage = 5L) {
    stop(
        "Dorado BAM parsing is not yet implemented in comma v0.2.0.\n",
        "Recommended workflow: convert Dorado BAM to modkit pileup BED first:\n",
        "  modkit pileup ", file, " output.bed --ref genome.fa\n",
        "Then load with: commaData(..., caller = 'modkit')"
    )
}
