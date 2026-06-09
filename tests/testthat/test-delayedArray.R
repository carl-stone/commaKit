test_that("commaData supports basic DelayedArray-backed assay operations", {
    skip_if_not_installed("DelayedArray")

    obj <- .make_two_modtype_fixture()
    for (nm in SummarizedExperiment::assayNames(obj)) {
        SummarizedExperiment::assay(obj, nm) <-
            DelayedArray::DelayedArray(SummarizedExperiment::assay(obj, nm))
    }

    expect_true(validObject(obj))
    expect_s4_class(methylation(obj), "DelayedMatrix")
    expect_s4_class(siteCoverage(obj), "DelayedMatrix")
    expect_s4_class(modCounts(obj), "DelayedMatrix")
    expect_s4_class(canonicalCounts(obj), "DelayedMatrix")

    sub <- obj[seq_len(3L), seq_len(2L)]
    expect_s4_class(sub, "commaData")
    expect_equal(dim(sub), c(3L, 2L))
    expect_s4_class(SummarizedExperiment::assay(sub, "methylation"),
                    "DelayedMatrix")
})
