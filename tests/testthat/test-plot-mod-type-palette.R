test_that("mod_type palette uses consistent colors for known modification types", {
    pal <- commaKit:::.modTypePalette(c("5mC", "6mA", "4mC", "6mA"))

    expect_named(pal, c("6mA", "5mC", "4mC"))
    expect_equal(pal[["6mA"]], "#e41a1c")
    expect_equal(pal[["5mC"]], "#377eb8")
    expect_equal(pal[["4mC"]], "#4daf4a")
})

test_that("mod_type palette warns and falls back for unknown modification types", {
    expect_warning(
        pal <- commaKit:::.modTypePalette(c("6mA", "9mX")),
        "No palette color defined for mod_type"
    )

    expect_equal(pal[["6mA"]], "#e41a1c")
    expect_equal(pal[["9mX"]], "grey50")
})

test_that("plot_genome_track applies shared mod_type colors", {
    data(comma_example_data)

    p <- plot_genome_track(
        comma_example_data,
        chromosome = "chr_sim",
        annotation = FALSE
    )
    color_scale <- p$scales$get_scales("colour")
    pal <- color_scale$palette(3L)

    expect_s3_class(p, "ggplot")
    expect_equal(pal[["6mA"]], "#e41a1c")
    expect_equal(pal[["5mC"]], "#377eb8")
})

test_that("plot_tss_profile applies shared mod_type colors", {
    data(comma_example_data)

    p <- suppressWarnings(plot_tss_profile(
        comma_example_data,
        feature_type = "gene",
        color_by = "mod_type",
        facet_by = "none"
    ))
    color_scale <- p$scales$get_scales("colour")
    pal <- color_scale$palette(3L)

    expect_s3_class(p, "ggplot")
    expect_equal(pal[["6mA"]], "#e41a1c")
    expect_equal(pal[["5mC"]], "#377eb8")
})
