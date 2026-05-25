NULL

# ─── diffMethyl design helpers ────────────────────────────────────────────────

#' Resolve the two-level differential methylation design contract
#'
#' Internal helper shared by diffMethyl() and all statistical backends. comma
#' currently supports one two-level contrast per diffMethyl() call. Multi-level
#' primary variables must be modeled in a future explicit-contrast API.
#'
#' @param coldata Sample metadata as a data.frame-like object.
#' @param formula One-sided design formula.
#' @param ref_level Optional reference level.
#'
#' @return A list containing primary_var, ref_level, treat_level, cond_levels,
#'   cond, and group_idx.
#' @keywords internal
.resolveDiffMethylDesign <- function(coldata, formula, ref_level = NULL) {
    if (!inherits(formula, "formula")) {
        stop("'formula' must be a formula object (e.g., ~ condition).")
    }

    rhs_vars <- all.vars(formula)
    if (length(rhs_vars) == 0L) {
        stop("'formula' must contain at least one RHS variable (e.g., ~ condition).")
    }
    primary_var <- rhs_vars[[1L]]

    coldata <- as.data.frame(coldata)
    if (!primary_var %in% colnames(coldata)) {
        stop(
            "Variable '", primary_var, "' from formula not found in sample metadata. ",
            "Available columns: ", paste(colnames(coldata), collapse = ", ")
        )
    }

    cond <- as.character(coldata[[primary_var]])
    if (anyNA(cond)) {
        stop("Column '", primary_var, "' contains NA values; diffMethyl() requires complete group labels.")
    }

    all_levels <- if (is.factor(coldata[[primary_var]])) {
        levels(coldata[[primary_var]])[levels(coldata[[primary_var]]) %in% cond]
    } else {
        sort(unique(cond))
    }

    if (length(all_levels) < 2L) {
        stop(
            "Differential methylation requires exactly 2 distinct levels of '",
            primary_var, "'. Found only: '", all_levels[[1L]], "'."
        )
    }
    if (length(all_levels) > 2L) {
        stop(
            "diffMethyl() currently supports exactly 2 levels for '", primary_var,
            "' per call. Found ", length(all_levels), " levels: ",
            paste(all_levels, collapse = ", "),
            ". For now, subset the object to the two groups you want to compare; ",
            "explicit multi-level contrasts are planned for a future API."
        )
    }

    if (!is.null(ref_level)) {
        if (!is.character(ref_level) || length(ref_level) != 1L || is.na(ref_level)) {
            stop("'reference' must be a single non-NA character string or NULL.")
        }
        if (!ref_level %in% all_levels) {
            stop(
                "'reference' value '", ref_level, "' not found in column '",
                primary_var, "'. Available values: ",
                paste(all_levels, collapse = ", ")
            )
        }
    } else {
        ref_level <- all_levels[[1L]]
    }

    cond_levels <- c(ref_level, setdiff(all_levels, ref_level))
    treat_level <- cond_levels[[2L]]
    group_idx <- lapply(cond_levels, function(lv) which(cond == lv))
    names(group_idx) <- cond_levels

    list(
        primary_var = primary_var,
        ref_level   = ref_level,
        treat_level = treat_level,
        cond_levels = cond_levels,
        cond        = cond,
        group_idx   = group_idx
    )
}

#' Compute per-group beta means and treatment-reference delta beta
#'
#' @param methyl_mat Sites x samples methylation matrix.
#' @param design Output from .resolveDiffMethylDesign().
#'
#' @return A list with group_means matrix and delta_beta vector.
#' @keywords internal
.computeDiffMethylGroupStats <- function(methyl_mat, design) {
    n_sites <- nrow(methyl_mat)
    group_means <- vapply(design$cond_levels, function(lv) {
        idx <- design$group_idx[[lv]]
        if (length(idx) == 1L) {
            methyl_mat[, idx]
        } else {
            rowMeans(methyl_mat[, idx, drop = FALSE], na.rm = TRUE)
        }
    }, numeric(n_sites))

    if (is.null(dim(group_means))) {
        group_means <- matrix(group_means, nrow = 1L,
                              dimnames = list(NULL, design$cond_levels))
    }
    group_means[is.nan(group_means)] <- NA_real_

    list(
        group_means = group_means,
        delta_beta  = group_means[, design$treat_level] -
                      group_means[, design$ref_level]
    )
}
