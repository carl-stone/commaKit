test_that("resultLayers() lists the default diffMethyl result layer", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 12L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")

    layers <- resultLayers(dm)

    expect_s4_class(layers, "DFrame")
    expect_equal(layers$name, "diffMethyl")
    expect_true(layers$is_default)
    expect_equal(layers$method, "quasi_f")
    expect_true("dm_pvalue" %in% as.character(layers$result_cols[[1L]]))
})

test_that("diffMethyl() stores multiple named result layers", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 12L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
        obj,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.min0",
        min_coverage = 0L
    )
    dm <- diffMethyl(
        dm,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.min1000",
        min_coverage = 1000L
    )

    layers <- resultLayers(dm)
    expect_equal(layers$name, c("quasi_f.min0", "quasi_f.min1000"))
    expect_equal(layers$name[layers$is_default], "quasi_f.min1000")

    res_low <- results(dm, result_name = "quasi_f.min0")
    res_high <- results(dm, result_name = "quasi_f.min1000")
    res_default <- results(dm)

    expect_true(any(!is.na(res_low$dm_pvalue)))
    expect_true(all(is.na(res_high$dm_pvalue)))
    expect_equal(res_default$dm_pvalue, res_high$dm_pvalue)
})

test_that("diffMethyl() protects explicit result layer names unless overwrite is set", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 8L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
        obj,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.v1"
    )

    expect_error(
        diffMethyl(
            dm,
            formula = ~ condition,
            method = "quasi_f",
            result_name = "quasi_f.v1"
        ),
        "already exists"
    )

    dm2 <- diffMethyl(
        dm,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.v1",
        min_coverage = 1000L,
        overwrite = TRUE
    )
    expect_true(all(is.na(results(dm2, result_name = "quasi_f.v1")$dm_pvalue)))
})



test_that("cannot overwrite active result while keeping stale default mirror", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 8L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
        obj,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.active"
    )

    expect_error(
        diffMethyl(
            dm,
            formula = ~ condition,
            method = "quasi_f",
            result_name = "quasi_f.active",
            min_coverage = 1000L,
            overwrite = TRUE,
            make_default = FALSE
        ),
        "Cannot overwrite the active"
    )
})

test_that("make_default = FALSE keeps the active result unchanged", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 10L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
        obj,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.active",
        min_coverage = 0L
    )
    active_before <- results(dm)$dm_pvalue
    dm <- diffMethyl(
        dm,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.archived",
        min_coverage = 1000L,
        make_default = FALSE
    )

    expect_equal(results(dm)$dm_pvalue, active_before)
    expect_true(all(is.na(results(dm, result = "quasi_f.archived")$dm_pvalue)))

    layers <- resultLayers(dm)
    expect_equal(layers$name[layers$is_default], "quasi_f.active")
})

test_that("make_default = FALSE can archive a run without creating a default", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 8L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
        obj,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.archived",
        make_default = FALSE
    )

    layers <- resultLayers(dm)
    expect_false(any(layers$is_default))
    expect_error(results(dm), "No default")
    expect_s3_class(results(dm, result_name = "quasi_f.archived"), "data.frame")
})

test_that("filterResults() can filter a selected result layer", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 12L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
        obj,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.loose",
        min_coverage = 0L
    )
    dm <- diffMethyl(
        dm,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.empty",
        min_coverage = 1000L
    )

    sig_loose <- filterResults(
        dm,
        padj = 1,
        delta_beta = 0,
        result_name = "quasi_f.loose"
    )
    sig_empty <- filterResults(
        dm,
        padj = 1,
        delta_beta = 0,
        name = "quasi_f.empty"
    )

    expect_gt(nrow(sig_loose), 0L)
    expect_equal(nrow(sig_empty), 0L)
})

test_that("row subsetting keeps named result layers aligned", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 12L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
        obj,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.v1"
    )
    dm <- diffMethyl(
        dm,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.empty",
        min_coverage = 1000L
    )

    keep <- c(2L, 4L, 6L)
    res_before <- results(dm, result_name = "quasi_f.v1")
    sub <- dm[keep, ]
    res_after <- results(sub, result_name = "quasi_f.v1")

    expect_equal(nrow(res_after), length(keep))
    expect_equal(res_after$dm_pvalue, res_before$dm_pvalue[keep])
    expect_true(all(is.na(results(sub, result_name = "quasi_f.empty")$dm_pvalue)))
})

test_that("results() validates selected result layer arguments", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 8L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(
        obj,
        formula = ~ condition,
        method = "quasi_f",
        result_name = "quasi_f.v1"
    )

    expect_error(results(dm, result_name = "missing"), "not found")
    expect_equal(
        results(dm, name = "quasi_f.v1")$dm_pvalue,
        results(dm, result_name = "quasi_f.v1")$dm_pvalue
    )
    expect_error(
        results(dm, result = "quasi_f.v1", result_name = "other"),
        "only one"
    )
    expect_error(
        results(dm, result = "quasi_f.v1", name = "other"),
        "only one"
    )
})



test_that("resultLayers() handles formulas whose deparse spans multiple lines", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 8L, n_ctrl = 2L, n_treat = 2L)
    long_var <- "very_long_batch_covariate_name_used_to_force_formula_deparse_wrapping"
    SummarizedExperiment::colData(obj)[[long_var]] <- c("a", "a", "b", "b")
    form <- stats::as.formula(paste("~ condition +", long_var))

    dm <- diffMethyl(
        obj,
        formula = form,
        method = "quasi_f",
        result_name = "quasi_f.long_formula"
    )

    layers <- resultLayers(dm)
    expect_length(layers$formula, 1L)
    expect_equal(layers$formula, paste(deparse(form), collapse = " "))
})

test_that("resultLayers() infers a legacy diffMethyl result row", {
    skip_if_not_installed("limma")
    obj <- .make_diff_methyl_fixture(n_sites = 8L, n_ctrl = 2L, n_treat = 2L)
    dm <- diffMethyl(obj, formula = ~ condition, method = "quasi_f")
    S4Vectors::metadata(dm)$diffMethyl_results <- NULL
    S4Vectors::metadata(dm)$diffMethyl_result_layers <- NULL
    S4Vectors::metadata(dm)$diffMethyl_default_result <- NULL

    layers <- resultLayers(dm)
    expect_equal(layers$name, "diffMethyl")
    expect_true(layers$is_default)
    expect_true("dm_padj" %in% as.character(layers$result_cols[[1L]]))
    expect_s3_class(results(dm), "data.frame")
})
