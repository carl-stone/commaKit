test_that("assayLayers() lists core assays and inferred defaults", {
    obj <- .make_two_modtype_fixture()
    layers <- assayLayers(obj)

    expect_s4_class(layers, "DFrame")
    expect_equal(layers$assay, SummarizedExperiment::assayNames(obj))
    expect_true(all(c("assay", "role", "type", "source", "is_default",
                      "default_for", "parent_assays", "method") %in%
                        colnames(layers)))
    expect_true(all(c("methylation", "coverage", "mod_counts",
                      "canonical_counts") %in% layers$assay))
    expect_true(layers$is_default[layers$assay == "methylation"])
    expect_equal(
        as.character(layers$default_for[layers$assay == "methylation"][[1L]]),
        "methylation"
    )
    expect_equal(layers$type[layers$assay == "mod_counts"], "reconstructed_counts")
})

test_that("assayLayers() infers registry rows for legacy objects", {
    obj <- .make_two_modtype_fixture()
    S4Vectors::metadata(obj)$assay_provenance <- NULL
    S4Vectors::metadata(obj)$assay_defaults <- NULL

    layers <- assayLayers(obj)

    expect_equal(layers$assay, SummarizedExperiment::assayNames(obj))
    expect_equal(layers$type[layers$assay == "methylation"], "filtered_beta")
    expect_true(layers$is_default[layers$assay == "coverage"])
})

test_that("commaKit:::.addAssayLayer() adds named derived layers without replacing raw assays", {
    obj <- .make_two_modtype_fixture()
    original_methyl <- methylation(obj)
    scaled <- log2((modCounts(obj) + 0.5) / (canonicalCounts(obj) + 0.5))

    obj2 <- commaKit:::.addAssayLayer(
        obj,
        assay_name = "methylation_mvalue.alpha05",
        value = scaled,
        type = "m_value",
        source = "mValues",
        role = "methylation_transform",
        parent_assays = c("mod_counts", "canonical_counts"),
        method = "mValues",
        params = list(alpha = 0.5),
        default_for = "methylation_transform"
    )

    expect_true("methylation_mvalue.alpha05" %in% SummarizedExperiment::assayNames(obj2))
    expect_equal(methylation(obj2), original_methyl)
    expect_equal(
        SummarizedExperiment::assay(obj2, "methylation_mvalue.alpha05"),
        scaled
    )

    layers <- assayLayers(obj2)
    added <- layers[layers$assay == "methylation_mvalue.alpha05", ]
    expect_equal(added$type, "m_value")
    expect_true(added$is_default)
    expect_equal(as.character(added$default_for[[1L]]), "methylation_transform")
    expect_equal(as.character(added$parent_assays[[1L]]),
                 c("mod_counts", "canonical_counts"))
})

test_that("commaKit:::.addAssayLayer() allows multiple explicit versions", {
    obj <- .make_two_modtype_fixture()
    layer1 <- methylation(obj) + 0.01
    layer2 <- methylation(obj) + 0.02

    obj <- commaKit:::.addAssayLayer(
        obj,
        assay_name = "methylation_norm.v1",
        value = layer1,
        type = "normalized_beta",
        source = "test",
        parent_assays = "methylation",
        method = "normalize_test"
    )
    obj <- commaKit:::.addAssayLayer(
        obj,
        assay_name = "methylation_norm.v2",
        value = layer2,
        type = "normalized_beta",
        source = "test",
        parent_assays = "methylation",
        method = "normalize_test"
    )

    expect_true(all(c("methylation_norm.v1", "methylation_norm.v2") %in%
                    SummarizedExperiment::assayNames(obj)))
    expect_false(any(duplicated(SummarizedExperiment::assayNames(obj))))
})

test_that("commaKit:::.addAssayLayer() protects existing layers unless overwrite is explicit", {
    obj <- .make_two_modtype_fixture()
    expect_error(
        commaKit:::.addAssayLayer(
            obj,
            assay_name = "methylation",
            value = methylation(obj),
            type = "filtered_beta",
            source = "test"
        ),
        "already exists"
    )

    obj2 <- commaKit:::.addAssayLayer(
        obj,
        assay_name = "methylation",
        value = methylation(obj),
        type = "filtered_beta",
        source = "test",
        overwrite = TRUE
    )
    expect_s4_class(obj2, "commaData")
})

test_that("commaKit:::.addAssayLayer() validates names, dimensions, and parents", {
    obj <- .make_two_modtype_fixture()
    expect_error(
        commaKit:::.addAssayLayer(obj, "bad name", methylation(obj), "x", "test"),
        "assay_name"
    )
    expect_error(
        commaKit:::.addAssayLayer(obj, "good_name", methylation(obj)[-1, ], "x", "test"),
        "same dimensions"
    )
    expect_error(
        commaKit:::.addAssayLayer(
            obj,
            "good_name",
            methylation(obj),
            "x",
            "test",
            parent_assays = "missing_layer"
        ),
        "parent_assays"
    )
})
