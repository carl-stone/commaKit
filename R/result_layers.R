#' @importFrom methods setGeneric setMethod is validObject
#' @importFrom SummarizedExperiment rowData rowRanges "rowRanges<-"
#' @importFrom S4Vectors DataFrame metadata "metadata<-"
#' @importFrom IRanges CharacterList
#' @importFrom GenomicRanges mcols "mcols<-"
NULL

.DIFFMETHYL_DEFAULT_RESULT_NAME <- "diffMethyl"

.DIFFMETHYL_CORE_RESULT_COLS <- c(
    "dm_pvalue",
    "dm_padj",
    "dm_delta_beta"
)

.validateResultLayerName <- function(result_name) {
    if (!is.character(result_name) || length(result_name) != 1L ||
            is.na(result_name) || !nzchar(result_name)) {
        stop("'result_name' must be a single non-empty character string.")
    }
    if (!grepl("^[A-Za-z][A-Za-z0-9_.:-]*$", result_name)) {
        stop(
            "'result_name' must start with a letter and contain only letters, ",
            "numbers, '_', '.', ':', or '-'."
        )
    }
    invisible(result_name)
}

.makeDiffMethylResultRecord <- function(result_name,
                                        result_cols,
                                        params,
                                        role = "diffMethyl",
                                        type = "differential_methylation",
                                        source = "diffMethyl",
                                        timestamp = Sys.time(),
                                        package_version = as.character(utils::packageVersion("commaKit"))) {
    .validateResultLayerName(result_name)
    if (!is.character(result_cols) || length(result_cols) == 0L ||
            any(is.na(result_cols) | !nzchar(result_cols))) {
        stop("'result_cols' must be a non-empty character vector.")
    }
    if (!is.list(params)) {
        stop("'params' must be a list.")
    }

    list(
        name = result_name,
        role = role,
        type = type,
        source = source,
        result_cols = as.character(result_cols),
        params = params,
        timestamp = timestamp,
        package_version = as.character(package_version)
    )
}

.diffMethylResultRegistry <- function(object) {
    md <- S4Vectors::metadata(object)
    registry <- md$diffMethyl_result_layers
    if (is.null(registry)) {
        registry <- list()
    }

    result_data <- md$diffMethyl_results
    if (!is.null(result_data)) {
        missing_records <- setdiff(names(result_data), names(registry))
        for (nm in missing_records) {
            registry[[nm]] <- .makeDiffMethylResultRecord(
                result_name = nm,
                result_cols = colnames(result_data[[nm]]),
                params = md$diffMethyl_params %||% list(),
                timestamp = NA
            )
        }
    }

    if (length(registry) == 0L && !is.null(md$diffMethyl_result_cols)) {
        registry[[.DIFFMETHYL_DEFAULT_RESULT_NAME]] <- .makeDiffMethylResultRecord(
            result_name = .DIFFMETHYL_DEFAULT_RESULT_NAME,
            result_cols = md$diffMethyl_result_cols,
            params = md$diffMethyl_params %||% list(),
            timestamp = NA
        )
    }

    registry
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.diffMethylDefaultResultName <- function(object) {
    md <- S4Vectors::metadata(object)
    default_name <- md$diffMethyl_default_result
    if (!is.null(default_name) && length(default_name) == 1L &&
            !is.na(default_name) && nzchar(default_name)) {
        return(as.character(default_name))
    }

    registry <- .diffMethylResultRegistry(object)
    if (!is.null(md$diffMethyl_result_cols)) {
        params_name <- md$diffMethyl_params$result_name
        if (!is.null(params_name) && length(params_name) == 1L &&
                !is.na(params_name) && params_name %in% names(registry)) {
            return(as.character(params_name))
        }
        if (.DIFFMETHYL_DEFAULT_RESULT_NAME %in% names(registry)) {
            return(.DIFFMETHYL_DEFAULT_RESULT_NAME)
        }
        if (length(registry) == 1L) {
            return(names(registry))
        }
    }

    NA_character_
}

.diffMethylResultNames <- function(object) {
    names(.diffMethylResultRegistry(object))
}

.hasDiffMethylResults <- function(object) {
    length(.diffMethylResultNames(object)) > 0L ||
        !is.null(S4Vectors::metadata(object)$diffMethyl_result_cols)
}

.resolveDiffMethylResultName <- function(object, result_name = NULL) {
    if (is.null(result_name)) {
        result_name <- .diffMethylDefaultResultName(object)
    }
    if (length(result_name) != 1L || is.na(result_name) || !nzchar(result_name)) {
        stop(
            "No default differential methylation result layer is set. ",
            "Run diffMethyl() first or provide 'result_name'."
        )
    }
    .validateResultLayerName(result_name)

    available <- .diffMethylResultNames(object)
    if (!result_name %in% available) {
        available_label <- if (length(available) > 0L) {
            paste(available, collapse = ", ")
        } else {
            "<none>"
        }
        stop(
            "Differential methylation result layer '", result_name,
            "' not found. Available result layers: ", available_label, "."
        )
    }

    result_name
}

.knownDiffMethylResultColsFromMetadata <- function(md) {
    cols <- character()
    if (!is.null(md$diffMethyl_result_cols)) {
        cols <- c(cols, md$diffMethyl_result_cols)
    }

    registry <- md$diffMethyl_result_layers
    if (!is.null(registry)) {
        for (record in registry) {
            cols <- c(cols, record$result_cols)
        }
    }

    result_data <- md$diffMethyl_results
    if (!is.null(result_data)) {
        for (df in result_data) {
            cols <- c(cols, colnames(df))
        }
    }

    cols <- unique(as.character(cols))
    cols[!is.na(cols) & nzchar(cols)]
}

.knownDiffMethylResultCols <- function(object) {
    .knownDiffMethylResultColsFromMetadata(S4Vectors::metadata(object))
}

.diffMethylResultData <- function(object, result_name) {
    md <- S4Vectors::metadata(object)
    result_data <- md$diffMethyl_results
    if (!is.null(result_data) && result_name %in% names(result_data)) {
        return(result_data[[result_name]])
    }

    if (result_name == .DIFFMETHYL_DEFAULT_RESULT_NAME &&
            !is.null(md$diffMethyl_result_cols)) {
        rd <- SummarizedExperiment::rowData(object)
        result_cols <- intersect(md$diffMethyl_result_cols, colnames(rd))
        if (length(result_cols) == length(md$diffMethyl_result_cols)) {
            return(S4Vectors::DataFrame(rd[, result_cols, drop = FALSE]))
        }
    }

    NULL
}

.setActiveDiffMethylResult <- function(object, result_data, result_cols) {
    rr <- SummarizedExperiment::rowRanges(object)
    mc <- GenomicRanges::mcols(rr)
    drop_cols <- intersect(
        .knownDiffMethylResultCols(object),
        colnames(mc)
    )
    if (length(drop_cols) > 0L) {
        keep_cols <- setdiff(colnames(mc), drop_cols)
        GenomicRanges::mcols(rr) <- mc[, keep_cols, drop = FALSE]
    }

    for (col_nm in result_cols) {
        GenomicRanges::mcols(rr)[[col_nm]] <- result_data[[col_nm]]
    }

    SummarizedExperiment::rowRanges(object) <- rr
    object
}

.addDiffMethylResultLayer <- function(object,
                                      result_name,
                                      result_data,
                                      params,
                                      result_cols,
                                      make_default = TRUE,
                                      overwrite = FALSE,
                                      timestamp = Sys.time()) {
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    .validateResultLayerName(result_name)
    if (!is(result_data, "DataFrame")) {
        result_data <- S4Vectors::DataFrame(result_data)
    }
    if (nrow(result_data) != nrow(object)) {
        stop("'result_data' must have one row per site in 'object'.")
    }
    missing_cols <- setdiff(result_cols, colnames(result_data))
    if (length(missing_cols) > 0L) {
        stop(
            "'result_data' is missing result column(s): ",
            paste(missing_cols, collapse = ", ")
        )
    }

    md <- S4Vectors::metadata(object)
    if (is.null(md$diffMethyl_results)) {
        md$diffMethyl_results <- list()
    }
    if (is.null(md$diffMethyl_result_layers)) {
        md$diffMethyl_result_layers <- list()
    }

    if (result_name %in% names(md$diffMethyl_results) && !isTRUE(overwrite)) {
        stop(
            "Differential methylation result layer '", result_name,
            "' already exists. Use a new 'result_name' or set overwrite = TRUE."
        )
    }

    md$diffMethyl_results[[result_name]] <- result_data[, result_cols, drop = FALSE]
    md$diffMethyl_result_layers[[result_name]] <- .makeDiffMethylResultRecord(
        result_name = result_name,
        result_cols = result_cols,
        params = params,
        timestamp = timestamp
    )

    S4Vectors::metadata(object) <- md

    if (isTRUE(make_default)) {
        object <- .setActiveDiffMethylResult(object, result_data, result_cols)
        md <- S4Vectors::metadata(object)
        md$diffMethyl_default_result <- result_name
        md$diffMethyl_result_cols <- result_cols
        md$diffMethyl_params <- params
        S4Vectors::metadata(object) <- md
    }

    methods::validObject(object)
    object
}

.subsetDiffMethylResultLayers <- function(md, row_index) {
    result_data <- md$diffMethyl_results
    if (is.null(result_data)) {
        return(md)
    }

    for (nm in names(result_data)) {
        df <- result_data[[nm]]
        if (nrow(df) >= max(row_index, 0L)) {
            result_data[[nm]] <- df[row_index, , drop = FALSE]
        }
    }
    md$diffMethyl_results <- result_data
    md
}

#' List differential methylation result layers
#'
#' Returns a registry of named differential methylation analyses stored in a
#' \code{\link{commaData}} object. Each call to \code{\link{diffMethyl}} can
#' write a named result layer, allowing multiple analysis runs to coexist while
#' the default layer remains available through \code{\link{results}}.
#'
#' @param object A \code{commaData} object.
#'
#' @return A \code{\link[S4Vectors]{DataFrame}} with one row per result layer
#'   and columns describing the layer name, default status, statistical method,
#'   result columns, filters, and provenance.
#'
#' @seealso \code{\link{diffMethyl}}, \code{\link{results}},
#'   \code{\link{filterResults}}
#'
#' @examples
#' data(comma_example_data)
#' dm <- diffMethyl(comma_example_data, formula = ~ condition,
#'                  method = "quasi_f", result_name = "quasi_f.v1")
#' resultLayers(dm)
#'
#' @export
setGeneric("resultLayers", function(object) standardGeneric("resultLayers"))

#' @rdname resultLayers
setMethod("resultLayers", "commaData", function(object) {
    registry <- .diffMethylResultRegistry(object)
    default_name <- .diffMethylDefaultResultName(object)

    records <- lapply(names(registry), function(nm) {
        record <- registry[[nm]]
        params <- record$params %||% list()
        list(
            name = nm,
            role = record$role %||% "diffMethyl",
            type = record$type %||% "differential_methylation",
            source = record$source %||% "diffMethyl",
            is_default = identical(nm, default_name),
            method = params$method %||% NA_character_,
            formula = paste(as.character(params$formula %||% NA_character_), collapse = " "),
            reference = params$reference %||% NA_character_,
            treatment = params$treatment %||% NA_character_,
            mod_context = as.character(params$mod_context %||% character()),
            mod_type = as.character(params$mod_type %||% character()),
            motif = as.character(params$motif %||% character()),
            p_adjust_method = params$p_adjust_method %||% NA_character_,
            min_coverage = params$min_coverage %||% NA_integer_,
            alpha = params$alpha %||% NA_real_,
            result_cols = as.character(record$result_cols %||% character()),
            timestamp = record$timestamp %||% as.POSIXct(NA_real_, origin = "1970-01-01"),
            package_version = record$package_version %||% NA_character_
        )
    })

    S4Vectors::DataFrame(
        name = unname(vapply(records, `[[`, character(1L), "name")),
        role = unname(vapply(records, `[[`, character(1L), "role")),
        type = unname(vapply(records, `[[`, character(1L), "type")),
        source = unname(vapply(records, `[[`, character(1L), "source")),
        is_default = unname(vapply(records, `[[`, logical(1L), "is_default")),
        method = unname(vapply(records, `[[`, character(1L), "method")),
        formula = unname(vapply(records, `[[`, character(1L), "formula")),
        reference = unname(vapply(records, `[[`, character(1L), "reference")),
        treatment = unname(vapply(records, `[[`, character(1L), "treatment")),
        mod_context = IRanges::CharacterList(
            lapply(records, `[[`, "mod_context")
        ),
        mod_type = IRanges::CharacterList(
            lapply(records, `[[`, "mod_type")
        ),
        motif = IRanges::CharacterList(
            lapply(records, `[[`, "motif")
        ),
        p_adjust_method = unname(vapply(records, `[[`, character(1L), "p_adjust_method")),
        min_coverage = unname(vapply(records, `[[`, integer(1L), "min_coverage")),
        alpha = unname(vapply(records, `[[`, numeric(1L), "alpha")),
        result_cols = IRanges::CharacterList(
            lapply(records, `[[`, "result_cols")
        ),
        timestamp = unname(vapply(records, function(record) {
            ts <- record$timestamp
            if (length(ts) == 0L || all(is.na(ts))) {
                NA_character_
            } else {
                format(ts[[1L]], usetz = TRUE)
            }
        }, character(1L))),
        package_version = unname(vapply(records, `[[`, character(1L), "package_version"))
    )
})
