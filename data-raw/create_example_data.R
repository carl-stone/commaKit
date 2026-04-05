## Script to generate comma_example_data — the synthetic example dataset
## bundled with the comma package for use in tests, vignettes, and examples.
##
## Run this script from the package root with:
##   source("data-raw/create_example_data.R")
##
## Output: data/comma_example_data.rda
##
## Design:
##   - Simulated 100 kb single-chromosome genome ("chr_sim")
##   - 6 samples: ctrl_1/2/3 (control), treat_1/2/3 (treatment)
##   - Two modification types: 6mA (~393 GATC sites) and 5mC (~195 CCWGG sites)
##   - ~30 differentially methylated 6mA sites between control and treatment
##   - Site positions: deterministic evenly-spaced backbone + random complement +
##     anchor positions guaranteeing coverage of all position-sensitive test windows
##   - Random coverage 10–150 per site per sample
##   - Object constructed directly (no file I/O) to avoid round-trip artifacts

set.seed(1312)

# ── Parameters ────────────────────────────────────────────────────────────────
GENOME_SIZE   <- 100000L
CHR_NAME      <- "chr_sim"
SAMPLES       <- c("ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3")
CONDITIONS    <- c("control", "control", "control", "treatment", "treatment", "treatment")
REPLICATES    <- c(1L, 2L, 3L, 1L, 2L, 3L)
GATC_SPACING  <- 512L    # deterministic backbone interval for 6mA sites (~195 sites)
CCWGG_SPACING <- 1024L   # deterministic backbone interval for 5mC sites (~97 sites)
N_RANDOM_6MA  <- 195L    # random sites added on top of backbone
N_RANDOM_5MC  <- 98L     # random sites added on top of backbone
N_DIFF_SITES  <- 30L     # 6mA sites that are differentially methylated
MOTIF_6MA     <- "GATC"  # simulated Dam methyltransferase motif for 6mA sites
MOTIF_5MC     <- "CCWGG" # simulated Dcm methyltransferase motif for 5mC sites
# Anchor positions that MUST be present so position-sensitive tests never skip.
# 1024 (= 512x2) is already on the backbone and covers the [1000, 1050] window.
GATC_ANCHORS  <- c(1984L, 49984L, 60032L)

# ── Simulate site positions ───────────────────────────────────────────────────
# Strategy: deterministic evenly-spaced backbone + random complement + anchors.
# Combined, backbone + random gives ~GENOME_SIZE/256 GATC sites and
# ~GENOME_SIZE/512 5mC sites, matching expected biological motif frequencies.

# 6mA (GATC): backbone every 512 bp → ~195 sites
gatc_backbone  <- seq(GATC_SPACING, GENOME_SIZE - 4L, by = GATC_SPACING)
gatc_taken     <- sort(unique(c(gatc_backbone, GATC_ANCHORS)))
gatc_pool      <- setdiff(seq_len(GENOME_SIZE - 4L), gatc_taken)
gatc_random    <- sample(gatc_pool, N_RANDOM_6MA)
gatc_positions <- sort(unique(c(gatc_backbone, GATC_ANCHORS, gatc_random)))

# 5mC (CCWGG): backbone every 1024 bp → ~97 sites
ccgg_backbone  <- seq(CCWGG_SPACING, GENOME_SIZE - 4L, by = CCWGG_SPACING)
ccgg_taken     <- sort(unique(c(ccgg_backbone, gatc_positions)))
ccgg_pool      <- setdiff(seq_len(GENOME_SIZE - 4L), ccgg_taken)
ccgg_random    <- sample(ccgg_pool, N_RANDOM_5MC)
ccgg_positions <- sort(unique(c(ccgg_backbone, ccgg_random)))

# Derive final site counts used everywhere downstream
N_6MA_SITES <- length(gatc_positions)
N_5MC_SITES <- length(ccgg_positions)

# Assign strands
gatc_strands  <- sample(c("+", "-"), N_6MA_SITES, replace = TRUE)
ccgg_strands  <- sample(c("+", "-"), N_5MC_SITES, replace = TRUE)

# ── Simulate methylation beta values ──────────────────────────────────────────
# 6mA: control samples highly methylated (~0.9), treatment mostly methylated
# but ~30 sites are significantly hypomethylated (~0.3)
diff_idx_6ma   <- sample.int(N_6MA_SITES, N_DIFF_SITES)

sim_6ma_beta <- function(n, is_diff, is_treatment) {
    betas <- numeric(n)
    for (i in seq_len(n)) {
        if (is_diff[i] && is_treatment) {
            # Differentially methylated: low in treatment
            betas[i] <- rbeta(1, 2, 6)   # mean ~0.25
        } else {
            # Constitutively methylated: high
            betas[i] <- rbeta(1, 18, 2)  # mean ~0.90
        }
    }
    pmin(pmax(betas, 0.01), 0.99)
}

is_diff_6ma <- seq_len(N_6MA_SITES) %in% diff_idx_6ma

beta_6ma_ctrl1  <- sim_6ma_beta(N_6MA_SITES, is_diff_6ma, is_treatment = FALSE)
beta_6ma_ctrl2  <- sim_6ma_beta(N_6MA_SITES, is_diff_6ma, is_treatment = FALSE)
beta_6ma_ctrl3  <- sim_6ma_beta(N_6MA_SITES, is_diff_6ma, is_treatment = FALSE)
beta_6ma_treat1 <- sim_6ma_beta(N_6MA_SITES, is_diff_6ma, is_treatment = TRUE)
beta_6ma_treat2 <- sim_6ma_beta(N_6MA_SITES, is_diff_6ma, is_treatment = TRUE)
beta_6ma_treat3 <- sim_6ma_beta(N_6MA_SITES, is_diff_6ma, is_treatment = TRUE)

# 5mC: no differential methylation in this example
beta_5mc_ctrl1  <- rbeta(N_5MC_SITES, 8, 2)  # mean ~0.80
beta_5mc_ctrl2  <- rbeta(N_5MC_SITES, 8, 2)
beta_5mc_ctrl3  <- rbeta(N_5MC_SITES, 8, 2)
beta_5mc_treat1 <- rbeta(N_5MC_SITES, 8, 2)
beta_5mc_treat2 <- rbeta(N_5MC_SITES, 8, 2)
beta_5mc_treat3 <- rbeta(N_5MC_SITES, 8, 2)

# ── Simulate coverage ─────────────────────────────────────────────────────────
cov_6ma_ctrl1  <- sample(10L:150L, N_6MA_SITES, replace = TRUE)
cov_6ma_ctrl2  <- sample(10L:150L, N_6MA_SITES, replace = TRUE)
cov_6ma_ctrl3  <- sample(10L:150L, N_6MA_SITES, replace = TRUE)
cov_6ma_treat1 <- sample(10L:150L, N_6MA_SITES, replace = TRUE)
cov_6ma_treat2 <- sample(10L:150L, N_6MA_SITES, replace = TRUE)
cov_6ma_treat3 <- sample(10L:150L, N_6MA_SITES, replace = TRUE)

cov_5mc_ctrl1  <- sample(10L:150L, N_5MC_SITES, replace = TRUE)
cov_5mc_ctrl2  <- sample(10L:150L, N_5MC_SITES, replace = TRUE)
cov_5mc_ctrl3  <- sample(10L:150L, N_5MC_SITES, replace = TRUE)
cov_5mc_treat1 <- sample(10L:150L, N_5MC_SITES, replace = TRUE)
cov_5mc_treat2 <- sample(10L:150L, N_5MC_SITES, replace = TRUE)
cov_5mc_treat3 <- sample(10L:150L, N_5MC_SITES, replace = TRUE)

# ── Build site keys ───────────────────────────────────────────────────────────
keys_6ma <- paste(CHR_NAME, gatc_positions, gatc_strands, "6mA", MOTIF_6MA, sep = ":")
keys_5mc <- paste(CHR_NAME, ccgg_positions, ccgg_strands, "5mC", MOTIF_5MC, sep = ":")
all_keys <- c(keys_6ma, keys_5mc)
n_total  <- length(all_keys)

# ── Build assay matrices ──────────────────────────────────────────────────────
methyl_mat <- matrix(
    c(c(beta_6ma_ctrl1,  beta_5mc_ctrl1),
      c(beta_6ma_ctrl2,  beta_5mc_ctrl2),
      c(beta_6ma_ctrl3,  beta_5mc_ctrl3),
      c(beta_6ma_treat1, beta_5mc_treat1),
      c(beta_6ma_treat2, beta_5mc_treat2),
      c(beta_6ma_treat3, beta_5mc_treat3)),
    nrow = n_total, ncol = 6L,
    dimnames = list(all_keys, SAMPLES)
)

coverage_mat <- matrix(
    c(c(cov_6ma_ctrl1,  cov_5mc_ctrl1),
      c(cov_6ma_ctrl2,  cov_5mc_ctrl2),
      c(cov_6ma_ctrl3,  cov_5mc_ctrl3),
      c(cov_6ma_treat1, cov_5mc_treat1),
      c(cov_6ma_treat2, cov_5mc_treat2),
      c(cov_6ma_treat3, cov_5mc_treat3)),
    nrow = n_total, ncol = 6L,
    dimnames = list(all_keys, SAMPLES)
)
storage.mode(coverage_mat) <- "integer"

# ── Build rowData ─────────────────────────────────────────────────────────────
row_df <- S4Vectors::DataFrame(
    chrom       = rep(CHR_NAME, n_total),
    position    = c(gatc_positions, ccgg_positions),
    strand      = c(gatc_strands, ccgg_strands),
    mod_type    = c(rep("6mA", N_6MA_SITES), rep("5mC", N_5MC_SITES)),
    motif       = c(rep(MOTIF_6MA, N_6MA_SITES), rep(MOTIF_5MC, N_5MC_SITES)),
    mod_context = c(rep(paste0("6mA_", MOTIF_6MA), N_6MA_SITES),
                    rep(paste0("5mC_", MOTIF_5MC), N_5MC_SITES)),
    is_diff     = c(is_diff_6ma, rep(FALSE, N_5MC_SITES)),  # ground truth for testing
    row.names   = all_keys
)

# ── Build colData ─────────────────────────────────────────────────────────────
col_df <- S4Vectors::DataFrame(
    sample_name = SAMPLES,
    condition   = CONDITIONS,
    replicate   = REPLICATES,
    caller      = rep("modkit", 6L),
    row.names   = SAMPLES
)

# ── Build annotation GRanges ──────────────────────────────────────────────────
ann_gr <- GenomicRanges::GRanges(
    seqnames = rep(CHR_NAME, 5L),
    ranges   = IRanges::IRanges(
        start = c(1L,    600L,  1400L, 2500L, 4000L),
        end   = c(500L,  1200L, 2000L, 3500L, 5000L)
    ),
    strand   = c("+", "+", "-", "+", "-")
)
GenomicRanges::mcols(ann_gr)$feature_type <- c("gene", "gene", "gene", "rRNA", "tRNA")
GenomicRanges::mcols(ann_gr)$name         <- c("geneA", "geneB", "geneC", "geneD", "geneE")

# ── Assemble commaData object ─────────────────────────────────────────────────
library(SummarizedExperiment)

se <- SummarizedExperiment(
    assays  = list(methylation = methyl_mat, coverage = coverage_mat),
    rowData = row_df,
    colData = col_df
)

comma_example_data <- new("commaData",
    se,
    genomeInfo = c(chr_sim = GENOME_SIZE),
    annotation = ann_gr,
    motifSites = GenomicRanges::GRanges()
)

# ── Save ──────────────────────────────────────────────────────────────────────
usethis::use_data(comma_example_data, overwrite = TRUE, compress = "xz")

message("comma_example_data saved. Object size: ",
        format(object.size(comma_example_data), units = "KB"))
message("Sites: ", nrow(comma_example_data), " (", N_6MA_SITES, " 6mA + ",
        N_5MC_SITES, " 5mC)")
message("Samples: ", ncol(comma_example_data))
message("Differentially methylated 6mA sites (ground truth): ", N_DIFF_SITES)
