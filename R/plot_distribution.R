#' @importFrom ggplot2 ggplot aes geom_density scale_x_continuous
#' @importFrom ggplot2 facet_wrap labs theme_bw
NULL

#' Plot methylation beta value distributions
#'
#' Produces a density plot of methylation beta values (0--1) for each sample
#' in a \code{\link{commaData}} object. Useful for QC and for comparing
#' methylation level distributions across samples and modification types.
#'
#' @param object A \code{\link{commaData}} object.
#' @param mod_type Character vector of modification types to include
#'   (e.g., \code{"6mA"}, \code{c("6mA", "5mC")}). If \code{NULL} (default), all
#'   modification types are included and the plot is faceted by
#'   \code{mod_type}.
#' @param motif Character vector or \code{NULL}. If provided, only sites with
#'   matching sequence context motif(s) are included (e.g., \code{"GATC"}).
#'   If \code{NULL} (default), all motifs are included.
#' @param mod_context Character vector or \code{NULL}. If provided, only sites
#'   whose \code{mod_context} rowData column matches one of the supplied values
#'   are included. If \code{NULL} (default), all modification contexts are
#'   included.
#' @param per_sample Logical. If \code{TRUE} (default), a separate density
#'   curve is drawn for each sample. If \code{FALSE}, a single aggregate
#'   density curve is drawn per modification type.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The x-axis shows beta
#'   values (0 = unmethylated, 1 = fully methylated); the y-axis shows kernel
#'   density. When \code{per_sample = TRUE}, curves are colored by
#'   \code{sample_name}. When multiple modification types are present (and
#'   \code{mod_type = NULL}), the plot is faceted by \code{mod_type}. Sites
#'   with \code{NA} beta values (below coverage threshold) are silently
#'   excluded.
#'
#' @examples
#' data(comma_example_data)
#' plot_methylation_distribution(comma_example_data)
#'
#' # One modification type only
#' plot_methylation_distribution(comma_example_data, mod_type = "6mA")
#'
#' # Aggregate across samples
#' plot_methylation_distribution(comma_example_data, per_sample = FALSE)
#'
#' @seealso \code{\link{methylomeSummary}}, \code{\link{plot_coverage}}
#'
#' @export
plot_methylation_distribution <- function(object,
                                          mod_type = NULL,
                                          motif = NULL,
                                          mod_context = NULL,
                                          per_sample = TRUE) {
  ## --- Input validation ---------------------------------------------------
  if (!is(object, "commaData")) {
    stop("'object' must be a commaData object.")
  }
  invalid_per_sample <- !is.logical(per_sample) ||
    length(per_sample) != 1L ||
    is.na(per_sample)
  if (invalid_per_sample) {
    stop("'per_sample' must be TRUE or FALSE.")
  }

  ## --- Optional site filters ----------------------------------------------
  object <- .applySiteFilters(
    object,
    mod_type = mod_type,
    motif = motif,
    mod_context = mod_context,
    caller = "plot_methylation_distribution()"
  )

  ## --- Extract data -------------------------------------------------------
  methyl_mat <- methylation(object)
  df <- .simpleQcAssayLongData(
    object,
    assay_matrix = methyl_mat,
    value_name = "beta",
    site_columns = "mod_type",
    empty_message = paste0(
      "No non-NA methylation values found after filtering. ",
      "Check coverage thresholds."
    )
  )

  ## --- Build ggplot -------------------------------------------------------
  multi_mod <- length(unique(df$mod_type)) > 1L

  if (per_sample) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[["beta"]],
        color = .data[["sample_name"]],
        fill = .data[["sample_name"]]
      )
    ) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::labs(color = "Sample", fill = "Sample")
  } else {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[["beta"]])
    ) +
      ggplot2::geom_density(
        fill = "steelblue",
        alpha = 0.4,
        color = "steelblue4"
      )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      expand = c(0.01, 0.01),
      name   = "Methylation"
    ) +
    ggplot2::labs(
      y     = "Density",
      title = "Methylation Beta Distribution"
    ) +
    ggplot2::theme_bw()

  if (multi_mod) {
    p <- p + ggplot2::facet_wrap("mod_type")
  }

  p
}
