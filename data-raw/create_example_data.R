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
##   - 3 samples: ctrl_1, ctrl_2 (condition "control"), treat_1 ("treatment")
##   - Two modification types: 6mA (GATC-like, ~200 sites) and 5mC (~100 sites)
##   - ~30 differentially methylated 6mA sites between control and treatment
##   - Random coverage 10–50 per site per sample
##   - Object constructed directly (no file I/O) to avoid round-trip artifacts

set.seed(1312)

# ── Parameters ────────────────────────────────────────────────────────────────
GENOME_SIZE  <- 100000L
CHR_NAME     <- "chr_sim"
SAMPLES      <- c("ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3")
CONDITIONS   <- c("control", "control", "control", "treatment", "treatment", "treatment")
REPLICATES   <- c(1L, 2L, 3L, 1L, 2L, 3L)
N_6MA_SITES  <- 200L
N_5MC_SITES  <- 100L
N_DIFF_SITES <- 30L   # 6mA sites that are differentially methylated
MOTIF_6MA    <- "GATC"   # simulated Dam methyltransferase motif for 6mA sites
MOTIF_5MC    <- "CCWGG"  # simulated Dcm methyltransferase motif for 5mC sites

# ── Simulate site positions ───────────────────────────────────────────────────
# 6mA sites (mimic GATC motif spacing — roughly every 256 bp on average)
gatc_positions <- sort(sample.int(GENOME_SIZE - 4L, N_6MA_SITES)) + 1L
# 5mC sites (mimic CCGG motif — slightly less frequent)
ccgg_positions <- sort(sample.int(GENOME_SIZE - 4L, N_5MC_SITES)) + 1L
# Assign strands
gatc_strands   <- sample(c("+", "-"), N_6MA_SITES, replace = TRUE)
ccgg_strands   <- sample(c("+", "-"), N_5MC_SITES, replace = TRUE)

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
