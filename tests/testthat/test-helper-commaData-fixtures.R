test_that("shared commaData fixture helper constructs valid objects", {
    beta <- matrix(
        c(0.1, 0.2, 0.3, 0.4),
        nrow = 2L,
        dimnames = list(NULL, c("s1", "s2"))
    )
    sample_info <- data.frame(
        sample_name = c("s1", "s2"),
        condition = c("ctrl", "treat"),
        replicate = 1:2,
        stringsAsFactors = FALSE
    )

    obj <- .make_commaData_fixture(
        beta,
        sample_info = sample_info,
        positions = c(10L, 20L),
        mod_type = c("6mA", "5mC"),
        motif = c("GATC", "CCWGG")
    )

    expect_s4_class(obj, "commaData")
    expect_equal(dim(methylation(obj)), c(2L, 2L))
    expect_setequal(as.character(siteInfo(obj)$mod_type), c("6mA", "5mC"))
    expect_setequal(motifs(obj), c("CCWGG", "GATC"))
})

test_that("shared two-mod-type fixture exposes both modification contexts", {
    obj <- .make_two_modtype_fixture(n_6ma = 3L, n_5mc = 2L)

    expect_s4_class(obj, "commaData")
    expect_equal(nrow(obj), 5L)
    expect_setequal(modTypes(obj), c("5mC", "6mA"))
    expect_setequal(modContexts(obj), c("5mC_CCWGG", "6mA_GATC"))
})

test_that("shared diffMethyl fixture includes ground-truth site metadata", {
    obj <- .make_diff_methyl_fixture(n_sites = 6L, n_ctrl = 2L, n_treat = 2L)
    si <- siteInfo(obj)

    expect_s4_class(obj, "commaData")
    expect_equal(dim(obj), c(6L, 4L))
    expect_true("is_diff" %in% colnames(si))
    expect_equal(sum(si$is_diff), 3L)
    expect_setequal(sampleInfo(obj)$condition, c("control", "treatment"))
})
