#' @importFrom methods setGeneric setMethod is
#' @importFrom SummarizedExperiment assay assayNames "assay<-"
#' @importFrom S4Vectors DataFrame metadata "metadata<-"
#' @importFrom IRanges CharacterList
NULL

.ASSAY_LAYER_DEFAULT_ROLES <- c(
    methylation = "methylation",
    coverage = "coverage",
    mod_counts = "mod_counts",
    canonical_counts = "canonical_counts"
)

.ASSAY_LAYER_DEFAULT_TYPES <- c(
    methylation = "filtered_beta",
    coverage = "observed_total_coverage",
    mod_counts = "observed_counts",
    canonical_counts = "observed_counts"
)

.ASSAY_LAYER_DEFAULT_SOURCES <- c(
    methylation = "unknown",
    coverage = "unknown",
    mod_counts = "unknown",
    canonical_counts = "unknown"
)

.validateAssayLayerName <- function(assay_name) {
    if (!is.character(assay_name) || length(assay_name) != 1L ||
            is.na(assay_name) || !nzchar(assay_name)) {
        stop("'assay_name' must be a single non-empty character string.")
    }
    if (!grepl("^[A-Za-z][A-Za-z0-9_.:-]*$", assay_name)) {
        stop(
            "'assay_name' must start with a letter and contain only letters, ",
            "numbers, '_', '.', ':', or '-'."
        )
    }
    invisible(assay_name)
}

.inferAssayLayerRole <- function(assay_name, record = NULL) {
    role <- record$role
    if (is.null(role) || length(role) == 0L || is.na(role) || !nzchar(role)) {
        role <- .ASSAY_LAYER_DEFAULT_ROLES[[assay_name]]
    }
    if (is.null(role) || length(role) == 0L || is.na(role) || !nzchar(role)) {
        role <- "derived"
    }
    as.character(role[[1L]])
}

.makeAssayLayerRecord <- function(type,
                                  source,
                                  role = "derived",
                                  parent_assays = character(),
                                  method = NA_character_,
                                  params = list(),
                                  default_for = character(),
                                  timestamp = NA,
                                  package_version = as.character(utils::packageVersion("commaKit"))) {
    if (is.null(type) || length(type) != 1L || is.na(type) || !nzchar(type)) {
        stop("'type' must be a single non-empty character string.")
    }
    if (is.null(source) || length(source) != 1L || is.na(source) || !nzchar(source)) {
        stop("'source' must be a single non-empty character string.")
    }
    if (!is.list(params)) {
        stop("'params' must be a list.")
    }

    list(
        type = as.character(type),
        source = as.character(source),
        role = as.character(role[[1L]]),
        parent_assays = as.character(parent_assays),
        method = if (length(method) == 0L) NA_character_ else as.character(method[[1L]]),
        params = params,
        default_for = as.character(default_for),
        timestamp = timestamp,
        package_version = as.character(package_version)
    )
}

.normalizeAssayLayerRecord <- function(assay_name, record = NULL, defaults = character()) {
    if (is.null(record)) {
        record <- list()
    }

    type <- record$type
    if (is.null(type) || length(type) == 0L || is.na(type) || !nzchar(type)) {
        type <- .ASSAY_LAYER_DEFAULT_TYPES[[assay_name]]
    }
    if (is.null(type) || length(type) == 0L || is.na(type) || !nzchar(type)) {
        type <- "derived"
    }

    source <- record$source
    if (is.null(source) || length(source) == 0L || is.na(source) || !nzchar(source)) {
        source <- .ASSAY_LAYER_DEFAULT_SOURCES[[assay_name]]
    }
    if (is.null(source) || length(source) == 0L || is.na(source) || !nzchar(source)) {
        source <- "unknown"
    }

    role <- .inferAssayLayerRole(assay_name, record)

    default_for <- names(defaults)[defaults == assay_name]

    parent_assays <- record$parent_assays
    if (is.null(parent_assays)) {
        parent_assays <- character()
    }

    method <- record$method
    if (is.null(method) || length(method) == 0L) {
        method <- NA_character_
    }

    timestamp <- record$timestamp
    if (is.null(timestamp) || length(timestamp) == 0L) {
        timestamp <- as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC")
    }

    package_version <- record$package_version
    if (is.null(package_version) || length(package_version) == 0L) {
        package_version <- NA_character_
    }

    params <- record$params
    if (is.null(params)) {
        params <- list()
    }

    list(
        assay = assay_name,
        type = as.character(type[[1L]]),
        source = as.character(source[[1L]]),
        role = role,
        is_default = assay_name %in% defaults,
        default_for = as.character(default_for),
        parent_assays = as.character(parent_assays),
        method = as.character(method[[1L]]),
        timestamp = timestamp,
        package_version = as.character(package_version[[1L]]),
        params = params
    )
}

.assayLayerDefaults <- function(object) {
    md <- S4Vectors::metadata(object)
    defaults <- md$assay_defaults
    if (is.null(defaults)) {
        defaults <- list()
    }
    if (is.list(defaults)) {
        defaults <- unlist(defaults, use.names = TRUE)
    } else {
        defaults <- as.character(defaults)
    }
    defaults <- defaults[nzchar(names(defaults)) & nzchar(defaults)]

    provenance <- md$assay_provenance
    if (!is.null(provenance)) {
        for (assay_name in names(provenance)) {
            default_for <- provenance[[assay_name]]$default_for
            default_for <- as.character(default_for)
            default_for <- default_for[nzchar(default_for)]
            missing_roles <- setdiff(default_for, names(defaults))
            defaults[missing_roles] <- assay_name
        }
    }

    present <- SummarizedExperiment::assayNames(object)
    inferred_roles <- .ASSAY_LAYER_DEFAULT_ROLES[
        names(.ASSAY_LAYER_DEFAULT_ROLES) %in% present
    ]
    missing_roles <- setdiff(names(inferred_roles), names(defaults))
    defaults <- c(defaults, inferred_roles[missing_roles])
    defaults[defaults %in% present]
}

#' List assay layers and defaults
#'
#' Returns a compact registry of assay matrices stored in a
#' \code{\link{commaData}} object. The registry is derived from
#' \code{assayNames(object)}, \code{metadata(object)$assay_provenance}, and
#' \code{metadata(object)$assay_defaults}. This mirrors layered-assay behavior
#' used by Bioconductor and single-cell workflows: raw layers remain present,
#' while derived layers are named explicitly and can be marked as defaults for a
#' role without overwriting the raw data.
#'
#' @param object A \code{commaData} object.
#'
#' @return A \code{\link[S4Vectors]{DataFrame}} with one row per assay and
#'   columns \code{assay}, \code{role}, \code{type}, \code{source},
#'   \code{is_default}, \code{default_for}, \code{parent_assays},
#'   \code{method}, \code{timestamp}, and \code{package_version}.
#'
#' @seealso \code{\link{assayProvenance}}, \code{\link{methylation}},
#'   \code{\link{siteCoverage}}, \code{\link{modCounts}},
#'   \code{\link{canonicalCounts}}
#'
#' @examples
#' data(comma_example_data)
#' assayLayers(comma_example_data)
#'
#' @export
setGeneric("assayLayers", function(object) standardGeneric("assayLayers"))

#' @rdname assayLayers
setMethod("assayLayers", "commaData", function(object) {
    assay_names <- SummarizedExperiment::assayNames(object)
    provenance <- S4Vectors::metadata(object)$assay_provenance
    if (is.null(provenance)) {
        provenance <- list()
    }
    defaults <- .assayLayerDefaults(object)

    records <- lapply(assay_names, function(nm) {
        .normalizeAssayLayerRecord(nm, provenance[[nm]], defaults)
    })
    names(records) <- assay_names

    S4Vectors::DataFrame(
        assay = unname(vapply(records, `[[`, character(1L), "assay")),
        role = unname(vapply(records, `[[`, character(1L), "role")),
        type = unname(vapply(records, `[[`, character(1L), "type")),
        source = unname(vapply(records, `[[`, character(1L), "source")),
        is_default = unname(vapply(records, `[[`, logical(1L), "is_default")),
        default_for = IRanges::CharacterList(
            lapply(records, `[[`, "default_for")
        ),
        parent_assays = IRanges::CharacterList(
            lapply(records, `[[`, "parent_assays")
        ),
        method = unname(vapply(records, `[[`, character(1L), "method")),
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

.addAssayLayer <- function(object,
                           assay_name,
                           value,
                           type,
                           source,
                           role = "derived",
                           parent_assays = character(),
                           method = NA_character_,
                           params = list(),
                           default_for = character(),
                           make_default = length(default_for) > 0L,
                           overwrite = FALSE,
                           timestamp = Sys.time()) {
    if (!is(object, "commaData")) {
        stop("'object' must be a commaData object.")
    }
    .validateAssayLayerName(assay_name)
    if (!is.matrix(value)) {
        stop("'value' must be a matrix.")
    }
    if (!identical(dim(value), dim(object))) {
        stop("'value' must have the same dimensions as 'object'.")
    }
    if (assay_name %in% SummarizedExperiment::assayNames(object) && !isTRUE(overwrite)) {
        stop(
            "Assay layer '", assay_name, "' already exists. Use a new assay ",
            "name for another version, or set overwrite = TRUE explicitly."
        )
    }

    parent_assays <- as.character(parent_assays)
    missing_parent <- setdiff(parent_assays, SummarizedExperiment::assayNames(object))
    if (length(missing_parent) > 0L) {
        stop(
            "'parent_assays' not found in object: ",
            paste(missing_parent, collapse = ", ")
        )
    }

    assay(object, assay_name, withDimnames = FALSE) <- value

    md <- S4Vectors::metadata(object)
    if (is.null(md$assay_provenance)) {
        md$assay_provenance <- list()
    }
    md$assay_provenance[[assay_name]] <- .makeAssayLayerRecord(
        type = type,
        source = source,
        role = role,
        parent_assays = parent_assays,
        method = method,
        params = params,
        default_for = if (isTRUE(make_default)) default_for else character(),
        timestamp = timestamp
    )

    if (isTRUE(make_default)) {
        default_for <- as.character(default_for)
        default_for <- default_for[nzchar(default_for)]
        if (length(default_for) > 0L) {
            defaults <- md$assay_defaults
            if (is.null(defaults)) {
                defaults <- list()
            }
            if (!is.list(defaults)) {
                defaults <- as.list(.assayLayerDefaults(object))
            }
            for (role_name in default_for) {
                defaults[[role_name]] <- assay_name
            }
            md$assay_defaults <- defaults
        }
    }

    S4Vectors::metadata(object) <- md
    methods::validObject(object)
    object
}
