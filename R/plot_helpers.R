.MOD_TYPE_PALETTE <- c(
  "6mA" = "#e41a1c",
  "5mC" = "#377eb8",
  "4mC" = "#4daf4a"
)

.modTypePalette <- function(mod_type, warn_unknown = TRUE) {
  mod_type <- unique(as.character(mod_type))
  mod_type <- mod_type[!is.na(mod_type)]

  unknown <- setdiff(mod_type, names(.MOD_TYPE_PALETTE))
  if (length(unknown) > 0L) {
    if (isTRUE(warn_unknown)) {
      warning(
        "No palette color defined for mod_type value(s): ",
        paste(unknown, collapse = ", "),
        ". Using grey50 fallback.",
        call. = FALSE
      )
    }
    fallback <- stats::setNames(rep("grey50", length(unknown)), unknown)
    return(c(.MOD_TYPE_PALETTE, fallback))
  }

  .MOD_TYPE_PALETTE
}

.simpleQcAssayLongData <- function(
  object,
  assay_matrix,
  value_name,
  site_columns = character(),
  empty_message
) {
  sample_info <- sampleInfo(object)
  sample_names <- colnames(assay_matrix)
  n_sites <- nrow(assay_matrix)
  n_samples <- length(sample_names)

  if (length(site_columns) > 0L) {
    site_info <- siteInfo(object)
    missing_site_columns <- setdiff(site_columns, colnames(site_info))
    if (length(missing_site_columns) > 0L) {
      stop(
        "Site column(s) not found: ",
        paste(missing_site_columns, collapse = ", "),
        "."
      )
    }
  }

  df <- data.frame(
    value = as.vector(assay_matrix),
    sample_name = rep(sample_names, each = n_sites),
    stringsAsFactors = FALSE
  )
  names(df)[names(df) == "value"] <- value_name

  for (site_column in site_columns) {
    df[[site_column]] <- rep(site_info[[site_column]], times = n_samples)
  }

  df <- df[!is.na(df[[value_name]]), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop(empty_message)
  }

  ## Join optional condition from sampleInfo when present. condition is not a
  ## commaData invariant; QC plots must still work for import/QC-only objects.
  sample_info_cols <- intersect(
    c("sample_name", "condition"),
    colnames(sample_info)
  )
  sample_info_sub <- sample_info[, sample_info_cols, drop = FALSE]
  merge(df, sample_info_sub, by = "sample_name", all.x = TRUE)
}
