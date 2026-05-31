test_that(".applySiteFilters filters by mod_type, motif, and mod_context", {
    data(comma_example_data)

    by_type <- commaKit:::.applySiteFilters(comma_example_data, mod_type = "6mA")
    expect_equal(nrow(by_type), 393L)
    expect_true(all(siteInfo(by_type)$mod_type == "6mA"))

    by_motif <- commaKit:::.applySiteFilters(comma_example_data, motif = "CCWGG")
    expect_equal(nrow(by_motif), 195L)
    expect_true(all(siteInfo(by_motif)$motif == "CCWGG"))

    by_context <- commaKit:::.applySiteFilters(
        comma_example_data,
        mod_context = "5mC_CCWGG"
    )
    expect_equal(nrow(by_context), 195L)
    expect_true(all(siteInfo(by_context)$mod_context == "5mC_CCWGG"))
})

test_that(".applySiteFilters validates filters sequentially", {
    data(comma_example_data)

    expect_error(
        commaKit:::.applySiteFilters(
            comma_example_data,
            mod_type = "6mA",
            motif = "CCWGG"
        ),
        "'motif' value\\(s\\) not found"
    )
    expect_error(
        commaKit:::.applySiteFilters(
            comma_example_data,
            motif = "GATC",
            mod_context = "5mC_CCWGG"
        ),
        "'mod_context' value\\(s\\) not found"
    )
})

test_that("exported site-filtering callers share validation behavior", {
    data(comma_example_data)

    expect_error(
        methylomeSummary(comma_example_data, mod_type = "6mA", motif = "CCWGG"),
        "'motif' value\\(s\\) not found"
    )
    expect_error(
        plot_coverage(comma_example_data, motif = "missing_motif"),
        "'motif' value\\(s\\) not found"
    )
    expect_error(
        mValues(comma_example_data, mod_context = "missing_context"),
        "'mod_context' value\\(s\\) not found"
    )
})
