#' **co**mma: for **m**icrobial **m**ethylation **a**nalysis
#'
#' The \pkg{comma} package provides a complete toolkit for genome-wide analysis
#' of bacterial DNA methylation from Oxford Nanopore sequencing data. It
#' supports the three major modification types (6mA, 5mC, 4mC), handles input
#' from modkit pileup (primary), Dorado BAM, and Megalodon (legacy) callers,
#' and provides a unified data container, annotation utilities, differential
#' methylation testing, and a full set of visualization functions.
#'
#' @section Main workflow:
#' \enumerate{
#'   \item \strong{Load data} — \code{\link{commaData}()} constructs the
#'     central \code{commaData} S4 object from per-sample methylation files.
#'   \item \strong{QC} — \code{\link{methylomeSummary}()},
#'     \code{\link{coverageDepth}()}, \code{\link{plot_coverage}()},
#'     \code{\link{plot_methylation_distribution}()}, and
#'     \code{\link{plot_pca}()} provide sample-level quality assessment.
#'   \item \strong{Annotate} — \code{\link{loadAnnotation}()} imports a GFF3
#'     or BED file; \code{\link{annotateSites}()} maps methylation sites to
#'     genomic features.
#'   \item \strong{Visualize} — \code{\link{plot_genome_track}()} and
#'     \code{\link{plot_metagene}()} show methylation in a genomic context.
#'   \item \strong{Differential methylation} — \code{\link{diffMethyl}()}
#'     tests each site; \code{\link{results}()} and
#'     \code{\link{filterResults}()} extract and filter the results table;
#'     \code{\link{plot_volcano}()} and \code{\link{plot_heatmap}()} visualize
#'     the findings.
#' }
#'
#' @section Key classes and constructors:
#' \describe{
#'   \item{\code{\link{commaData}}}{The central S4 data container, extending
#'     \code{SummarizedExperiment}. Stores methylation (beta) and coverage
#'     matrices, per-site and per-sample metadata, genome information,
#'     genomic annotation, and motif site locations.}
#' }
#'
#' @section Package options:
#' None. All parameters are passed directly to individual functions.
#'
#' @return No return value. This page provides package-level documentation.
#'   See individual function pages for return values.
#'
#' @references
#' The modkit pileup format is documented at
#' \url{https://nanoporetech.github.io/modkit/}.
#'
#' @name comma-package
#' @aliases comma
#' @docType package
#' @keywords package
"_PACKAGE"

## Suppress R CMD check NOTEs for rlang .data pronoun used in ggplot2 aes()
## calls throughout the plot_* functions.
utils::globalVariables(".data")
