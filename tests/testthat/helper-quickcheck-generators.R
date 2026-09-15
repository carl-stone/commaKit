# quickcheck generators for commaKit property-based tests.
#
# These are NOT quickcheck generators themselves (quickcheck/hedgehog
# generators are opaque hedgehog gen objects). Instead they are ordinary
# functions that, given randomly drawn scalars produced by quickcheck
# generators, build valid or deliberately invalid commaData components.
# They are used inside for_all(property = ...) bodies via integer_(),
# character_(), double_(), etc.
#
# Invariant-first design: every *_valid_* helper builds a component that
# MUST be accepted; every *_invalid_* helper builds a component that MUST
# be rejected with an actionable error (naming the offending argument).

# Valid mod types (mirrors .VALID_MOD_TYPES in R/commaData_class.R)
.qc_valid_mod_types <- c("6mA", "5mC", "4mC")

.qc_mod_type_gen <- function() {
  # A generator for a single valid mod_type string.
  quickcheck::constant(sample(.qc_valid_mod_types, 1L))
}

.qc_sample_names_gen <- function(n) {
  # Distinct sample names of length n. quickcheck's character_() can
  # produce duplicates and empty strings, so we derive names from an
  # index integer to guarantee uniqueness.
  quickcheck::constant(paste0("s", seq_len(n)))
}

.qc_positions_gen <- function(n, span = 100000L) {
  # n distinct, sorted, positive positions.
  quickcheck::constant(sort(sample.int(span, n)))
}

.qc_strand_gen <- function() {
  quickcheck::constant(sample(c("+", "-"), 1L))
}

.qc_motif_gen <- function() {
  quickcheck::constant(sample(c("GATC", "CCWGG", "CCWGG", "GATC"), 1L))
}

# ── bedMethyl row construction ──────────────────────────────────────────────
# 18-column modkit pileup bedMethyl row (tab-separated):
# chrom start end mod_code score strand thickStart thickEnd itemRgb
# Nvalid_cov fraction_modified Nmod Ncanonical Nother_mod Ndelete Nfail Ndiff Nnocall

.qc_make_bedmethyl_row <- function(
  chrom,
  position, # 1-based
  strand,
  mod_type,
  motif,
  coverage,
  beta
) {
  mod_code <- switch(
    mod_type,
    "6mA" = "a",
    "5mC" = "m",
    "4mC" = "21839",
    stop("unknown mod_type: ", mod_type)
  )
  mod_code <- paste0(mod_code, ",", motif, ",1")
  nmod <- round(coverage * beta)
  paste(
    chrom,
    position - 1L, # BED is 0-based
    position,
    mod_code,
    0L,
    strand,
    position - 1L,
    position,
    0,
    coverage,
    round(beta * 100, 1L),
    nmod,
    coverage - nmod,
    0L,
    0L,
    0L,
    0L,
    0L,
    sep = "\t"
  )
}

.qc_write_bedmethyl <- function(path, lines) {
  writeLines(lines, path)
  path
}

# ── valid commaData components ──────────────────────────────────────────────

.qc_valid_files_coldata <- function(n_samples) {
  # Returns list(files = named chr, colData = data.frame) for n_samples
  # valid samples backed by real (written-to-disk) bedMethyl files.
  # Each sample gets one site at a fixed position so parsing succeeds.
  sample_names <- paste0("s", seq_len(n_samples))
  paths <- file.path(
    tempdir(),
    paste0("qc-bed-", sample_names, "-", sample.int(2^30, n_samples), ".bed")
  )
  for (i in seq_len(n_samples)) {
    lines <- .qc_make_bedmethyl_row(
      chrom = "chr_qc",
      position = 1000L,
      strand = "+",
      mod_type = "6mA",
      motif = "GATC",
      coverage = 20L,
      beta = 0.5
    )
    .qc_write_bedmethyl(paths[[i]], lines)
  }
  list(
    files = stats::setNames(paths, sample_names),
    colData = data.frame(
      sample_name = sample_names,
      replicate = seq_len(n_samples),
      stringsAsFactors = FALSE
    )
  )
}

# ── invalid commaData components ────────────────────────────────────────────

.qc_bad_coldata_missing_col <- function() {
  data.frame(
    not_sample_name = "s1",
    replicate = 1L,
    stringsAsFactors = FALSE
  )
}

.qc_invalid_files_unnamed <- function() {
  c("/tmp/does-not-exist-a.bed", "/tmp/does-not-exist-b.bed")
}

.qc_bad_files_missing_sample <- function() {
  stats::setNames(
    c("/tmp/does-not-exist-a.bed", "/tmp/does-not-exist-b.bed"),
    c("s1", "ghost")
  )
}

.qc_invalid_files_extra_sample <- function(coldata) {
  stats::setNames(
    "/tmp/does-not-exist-a.bed",
    "s1"
  )
}

# ── fixture-based misuse harness ─────────────────────────────────────────────
# These wrap the existing .make_commaData_fixture() to generate random-but-
# valid commaData objects, then call the API with type-confused arguments.

.qc_random_fixture <- function(n_sites, n_samples, mod_type = "6mA") {
  sample_names <- paste0("s", seq_len(n_samples))
  beta <- matrix(
    stats::runif(n_sites * n_samples, 0.05, 0.95),
    nrow = n_sites,
    dimnames = list(NULL, sample_names)
  )
  coverage <- matrix(
    as.integer(20L + sample.int(30L, n_sites * n_samples, replace = TRUE) - 1L),
    nrow = n_sites,
    ncol = n_samples,
    dimnames = list(NULL, sample_names)
  )
  sample_info <- data.frame(
    sample_name = sample_names,
    condition = rep(c("control", "treatment"), length.out = n_samples),
    replicate = seq_len(n_samples),
    stringsAsFactors = FALSE
  )
  .make_commaData_fixture(
    beta = beta,
    coverage = coverage,
    sample_info = sample_info,
    positions = sort(sample.int(100000L, n_sites)),
    mod_type = mod_type,
    motif = "GATC"
  )
}
