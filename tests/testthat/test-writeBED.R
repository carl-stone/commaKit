# Tests for writeBED() — Phase 3 BED export function

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_bed_data <- function() {
    # 5 sites at known positions with beta values spanning all 5 RGB color bands:
    #   0.1 → score 100 ≤ 200 → blue     (0,0,255)
    #   0.3 → score 300 ≤ 400 → blue-pur (83,0,172)
    #   0.5 → score 500 ≤ 600 → purple   (167,0,85)
    #   0.7 → score 700 ≤ 800 → red-pur  (222,0,28)
    #   0.9 → score 900 ≤ 1000→ red      (250,0,0)
    n_sites   <- 5L
    positions <- c(1000L, 2000L, 3000L, 4000L, 5000L)
    site_keys <- paste0("chr_sim:", positions, ":+:6mA:GATC")
    betas     <- c(0.1, 0.3, 0.5, 0.7, 0.9)

    methyl_mat <- matrix(
        c(betas, betas + 0.05),
        nrow = n_sites, ncol = 2L,
        dimnames = list(site_keys, c("samp1", "samp2"))
    )
    # Clamp samp2 to [0,1]
    methyl_mat[, "samp2"] <- pmin(1, pmax(0, methyl_mat[, "samp2"]))

    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2L,
                      dimnames = dimnames(methyl_mat))

    rd <- S4Vectors::DataFrame(
        chrom    = rep("chr_sim", n_sites),
        position = positions,
        strand   = rep("+", n_sites),
        mod_type = rep("6mA", n_sites),
        motif    = rep("GATC", n_sites),
        row.names = site_keys
    )
    cd <- S4Vectors::DataFrame(
        sample_name = c("samp1", "samp2"),
        condition   = c("ctrl", "treat"),
        replicate   = 1:2,
        row.names   = c("samp1", "samp2")
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(methylation = methyl_mat, coverage = cov_mat),
        rowData = rd,
        colData = cd
    )
    new("commaData", se,
        genomeInfo = c(chr_sim = 10000L),
        annotation = GenomicRanges::GRanges(),
        motifSites = GenomicRanges::GRanges())
}

# ─── File creation & return value ─────────────────────────────────────────────

test_that("writeBED: creates a file at the specified path", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    expect_true(file.exists(f))
})

test_that("writeBED: returns the file path invisibly", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    ret <- writeBED(obj, file = f, sample = "samp1")
    expect_equal(ret, f)
})

# ─── Track header ─────────────────────────────────────────────────────────────

test_that("writeBED: first line is a track header", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    lines <- readLines(f)
    expect_true(grepl("^track ", lines[1]))
})

test_that("writeBED: track header contains sample name by default", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    header <- readLines(f)[1]
    expect_true(grepl('name="samp1"', header))
})

test_that("writeBED: custom track_name appears in header", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1", track_name = "my_track")
    header <- readLines(f)[1]
    expect_true(grepl('name="my_track"', header))
})

test_that("writeBED: custom track_description appears in header", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1", track_description = "my desc")
    header <- readLines(f)[1]
    expect_true(grepl('description="my desc"', header))
})

test_that("writeBED: track header contains itemRgb='On'", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    header <- readLines(f)[1]
    expect_true(grepl('itemRgb="On"', header, ignore.case = TRUE))
})

# ─── BED line format ──────────────────────────────────────────────────────────

test_that("writeBED: each BED data line has 9 tab-separated columns", {
    obj   <- .make_bed_data()
    f     <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    lines  <- readLines(f)[-1]  # skip track header
    n_tabs <- vapply(lines, function(l) nchar(gsub("[^\t]", "", l)), integer(1L),
                     USE.NAMES = FALSE)
    expect_true(all(n_tabs == 8L))  # 9 columns = 8 tabs
})

test_that("writeBED: number of BED data lines equals non-NA site count", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    lines    <- readLines(f)[-1]
    n_non_na <- sum(!is.na(methylation(obj)[, "samp1"]))
    expect_equal(length(lines), n_non_na)
})

# ─── Coordinate conversion ────────────────────────────────────────────────────

test_that("writeBED: chromStart is 0-based (position - 1) and chromEnd is 1-based", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    lines     <- readLines(f)[-1]
    first_row <- strsplit(lines[1], "\t")[[1]]
    # First site has position = 1000; BED chromStart = 999, chromEnd = 1000
    expect_equal(as.integer(first_row[2]), 999L)   # chromStart (0-based)
    expect_equal(as.integer(first_row[3]), 1000L)  # chromEnd
})

# ─── Score field ──────────────────────────────────────────────────────────────

test_that("writeBED: score equals round(beta * 1000) for each site", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    lines     <- readLines(f)[-1]
    score_col <- vapply(strsplit(lines, "\t"),
                        function(x) as.integer(x[5]), integer(1L),
                        USE.NAMES = FALSE)
    betas    <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    expected <- as.integer(round(betas * 1000))
    expect_equal(score_col, expected)
})

# ─── RGB color scale ──────────────────────────────────────────────────────────

test_that("writeBED: rgb_scale=TRUE assigns correct colors for each score band", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1", rgb_scale = TRUE)
    lines   <- readLines(f)[-1]
    rgb_col <- vapply(strsplit(lines, "\t"),
                      function(x) x[9], character(1L),
                      USE.NAMES = FALSE)
    # score 100 ≤ 200  → blue
    expect_equal(rgb_col[1], "0,0,255")
    # score 300 ≤ 400  → blue-purple
    expect_equal(rgb_col[2], "83,0,172")
    # score 500 ≤ 600  → purple
    expect_equal(rgb_col[3], "167,0,85")
    # score 700 ≤ 800  → red-purple
    expect_equal(rgb_col[4], "222,0,28")
    # score 900 ≤ 1000 → red
    expect_equal(rgb_col[5], "250,0,0")
})

test_that("writeBED: rgb_scale=FALSE gives '0,0,0' for all sites", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1", rgb_scale = FALSE)
    lines   <- readLines(f)[-1]
    rgb_col <- vapply(strsplit(lines, "\t"),
                      function(x) x[9], character(1L),
                      USE.NAMES = FALSE)
    expect_true(all(rgb_col == "0,0,0"))
})

# ─── NA site exclusion ────────────────────────────────────────────────────────

test_that("writeBED: NA methylation sites are excluded from output", {
    obj   <- .make_bed_data()
    m     <- methylation(obj)
    m[1L, "samp1"] <- NA_real_
    SummarizedExperiment::assay(obj, "methylation") <- m
    f <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(obj, file = f, sample = "samp1")
    lines <- readLines(f)[-1]
    expect_equal(length(lines), 4L)  # 5 sites minus 1 NA
})

test_that("writeBED: all-NA sample produces only track header line with warning", {
    obj <- .make_bed_data()
    m   <- methylation(obj)
    m[, "samp1"] <- NA_real_
    SummarizedExperiment::assay(obj, "methylation") <- m
    f <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    expect_warning(writeBED(obj, file = f, sample = "samp1"), "No sites")
    lines <- readLines(f)
    expect_equal(length(lines), 1L)  # only the track header
})

# ─── mod_type filtering ───────────────────────────────────────────────────────

test_that("writeBED: mod_type filtering writes only matching sites", {
    data(comma_example_data)
    f <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    writeBED(comma_example_data, file = f, sample = "ctrl_1", mod_type = "6mA")
    lines  <- readLines(f)[-1]
    n_6mA  <- sum(!is.na(
        methylation(comma_example_data)[
            siteInfo(comma_example_data)$mod_type == "6mA", "ctrl_1"
        ]
    ))
    expect_equal(length(lines), n_6mA)
})

test_that("writeBED: mod_type=NULL writes all sites", {
    data(comma_example_data)
    f_all  <- tempfile(fileext = ".bed")
    f_6mA  <- tempfile(fileext = ".bed")
    on.exit({ unlink(f_all); unlink(f_6mA) })
    writeBED(comma_example_data, file = f_all, sample = "ctrl_1")
    writeBED(comma_example_data, file = f_6mA, sample = "ctrl_1", mod_type = "6mA")
    n_all <- length(readLines(f_all)) - 1L   # subtract header
    n_6mA <- length(readLines(f_6mA)) - 1L
    expect_gt(n_all, n_6mA)
})

# ─── Error handling ───────────────────────────────────────────────────────────

test_that("writeBED: error on non-commaData input", {
    f <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    expect_error(
        writeBED(data.frame(x = 1), file = f, sample = "s1"),
        "'object' must be a commaData"
    )
})

test_that("writeBED: error on invalid sample name", {
    obj <- .make_bed_data()
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    expect_error(
        writeBED(obj, file = f, sample = "no_such_sample"),
        "not found in object"
    )
})

test_that("writeBED: error when mod_type has no matching sites", {
    obj <- .make_bed_data()   # only has 6mA sites
    f   <- tempfile(fileext = ".bed")
    on.exit(unlink(f))
    expect_error(
        writeBED(obj, file = f, sample = "samp1", mod_type = "5mC"),
        "No sites remain"
    )
})
