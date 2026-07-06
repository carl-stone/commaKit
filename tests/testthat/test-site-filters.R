test_that("filterSites() applies combined site filters in object row order", {
  obj <- .make_two_modtype_fixture(n_6ma = 4L, n_5mc = 3L)
  info <- as.data.frame(siteInfo(obj))
  expected <- which(
    info$mod_type == "5mC" &
      info$motif == "CCWGG" &
      info$mod_context == "5mC_CCWGG"
  )

  filtered <- filterSites(
    obj,
    mod_type = "5mC",
    motif = "CCWGG",
    mod_context = "5mC_CCWGG"
  )

  expect_equal(nrow(filtered), length(expected))
  expect_equal(siteInfo(filtered)$site_key, info$site_key[expected])
  expect_true(all(siteInfo(filtered)$mod_context == "5mC_CCWGG"))
})

test_that("filterSites() keeps current empty-result behavior", {
  obj <- .make_two_modtype_fixture(n_6ma = 4L, n_5mc = 3L)

  filtered <- filterSites(obj, mod_type = "6mA", motif = "CCWGG")

  expect_s4_class(filtered, "commaData")
  expect_equal(nrow(filtered), 0L)
})

test_that("exported callers validate filters sequentially", {
  obj <- .make_two_modtype_fixture(n_6ma = 4L, n_5mc = 3L)

  expect_error(
    methylomeSummary(
      obj,
      mod_type = "6mA",
      motif = "CCWGG"
    ),
    "'motif' value\\(s\\) not found"
  )
  expect_error(
    mValues(
      obj,
      motif = "GATC",
      mod_context = "5mC_CCWGG"
    ),
    "'mod_context' value\\(s\\) not found"
  )
})

test_that("exported site-filtering callers share validation behavior", {
  obj <- .make_two_modtype_fixture(n_6ma = 4L, n_5mc = 3L)

  expect_error(
    methylomeSummary(obj, mod_type = "6mA", motif = "CCWGG"),
    "'motif' value\\(s\\) not found"
  )
  expect_error(
    plot_coverage(obj, motif = "missing_motif"),
    "'motif' value\\(s\\) not found"
  )
  expect_error(
    mValues(obj, mod_context = "missing_context"),
    "'mod_context' value\\(s\\) not found"
  )
})
