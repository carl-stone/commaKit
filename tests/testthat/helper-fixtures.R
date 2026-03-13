# Shared test fixtures for CoMMA test suite.
#
# These objects are loaded automatically by testthat before every test file.
# They represent a minimal synthetic dataset — a small 10 kb "genome" with
# three methylation sites, two genomic features, and one Transcription-Unit —
# designed to make expected outputs easy to reason about by hand.

# ---------------------------------------------------------------------------
# Minimal methylation dataframe
# Three GATC sites at positions 500, 1500, and 3000.
# ---------------------------------------------------------------------------
fixture_methyl_df <- data.frame(
  Position = c(500L, 1500L, 3000L),
  beta     = c(0.2, 0.8, 0.5),
  Strand   = c("+", "+", "-"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Minimal feature annotation dataframe
# Two features: a Gene spanning 400–700 and a Promoter spanning 200–450.
# Position 500 falls inside Gene (and NOT Promoter: 500 > 450).
# Position 1500 falls inside neither → No_Feature.
# Position 3000 falls inside neither → No_Feature.
# ---------------------------------------------------------------------------
fixture_meta_df <- data.frame(
  Type   = c("Gene",    "Promoter"),
  Site   = c("geneA",   "promA"),
  Left   = c(400L,      200L),
  Right  = c(700L,      450L),
  Strand = c("+",       "+"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Transcription-Unit annotation for TSS / TTS tests
# Two TUs:
#   tuA  sense (+):     Left=1000, Right=2000
#   tuB  antisense (-): Left=2500, Right=3500
# ---------------------------------------------------------------------------
fixture_tu_df <- data.frame(
  Type   = c("Transcription-Units", "Transcription-Units"),
  Site   = c("tuA",                 "tuB"),
  Left   = c(1000L,                 2500L),
  Right  = c(2000L,                 3500L),
  Strand = c("+",                   "-"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Terminator annotation for annotateTTS(annotatedOnly = TRUE) tests
# One Rho-Independent terminator on sense strand at 2100.
# ---------------------------------------------------------------------------
fixture_term_df <- data.frame(
  Type   = c("Rho-Independent-Terminators", "Transcription-Units"),
  Site   = c("termA",                       "tuA"),
  Left   = c(2050L,                         1000L),
  Right  = c(2150L,                         2000L),
  Strand = c("+",                           "+"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Coverage dataframe for calculateMethylSiteDepth tests
# Positions 1–100 with uniform coverage of 10.
# ---------------------------------------------------------------------------
fixture_cov_df <- data.frame(
  Position = 1:100,
  coverage = rep(10L, 100),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Dataset for varByCoverage tests
# Three sites each at coverage 5 and 10, with known methylation deltas.
# ---------------------------------------------------------------------------
fixture_var_df <- data.frame(
  Coverage_Sample      = c(5L, 5L, 5L, 10L, 10L, 10L),
  Percent_Methyl_Sample = c(0.6, 0.7, 0.8,  0.4,  0.5,  0.6),
  Ancestor_Mean         = c(0.5, 0.5, 0.5,  0.5,  0.5,  0.5),
  stringsAsFactors = FALSE
)
