.make_commaData_fixture <- function(beta, coverage = NULL, sample_info = NULL,
                                    positions = NULL, mod_type = "6mA",
                                    motif = "GATC", chrom = "chr_sim",
                                    seqlength = 100000L, strand = "+",
                                    site_metadata = NULL,
                                    mod_counts = NULL,
                                    canonical_counts = NULL,
                                    other_mod_counts = NULL) {
    stopifnot(is.matrix(beta))

    n_sites <- nrow(beta)
    sample_names <- colnames(beta)
    if (is.null(sample_names)) {
        sample_names <- paste0("sample_", seq_len(ncol(beta)))
        colnames(beta) <- sample_names
    }

    if (is.null(coverage)) {
        coverage <- matrix(20L, nrow = n_sites, ncol = ncol(beta),
                           dimnames = dimnames(beta))
    }
    stopifnot(is.matrix(coverage), identical(dim(coverage), dim(beta)))
    if (is.null(dimnames(coverage))) {
        dimnames(coverage) <- dimnames(beta)
    }
    if (is.null(mod_counts)) {
        mod_counts <- round(beta * coverage)
    }
    mod_counts[is.na(mod_counts)] <- NA_integer_
    mod_counts <- pmax(0L, pmin(as.integer(mod_counts), coverage))
    dim(mod_counts) <- dim(beta)
    dimnames(mod_counts) <- dimnames(beta)

    if (is.null(canonical_counts)) {
        canonical_counts <- coverage - mod_counts
    }
    canonical_counts[is.na(canonical_counts)] <- NA_integer_
    canonical_counts <- pmax(0L, as.integer(canonical_counts))
    dim(canonical_counts) <- dim(beta)
    dimnames(canonical_counts) <- dimnames(beta)

    if (is.null(other_mod_counts)) {
        other_mod_counts <- matrix(0L, nrow = n_sites, ncol = ncol(beta),
                                   dimnames = dimnames(beta))
    }
    other_mod_counts[is.na(other_mod_counts)] <- NA_integer_
    other_mod_counts <- pmax(0L, as.integer(other_mod_counts))
    dim(other_mod_counts) <- dim(beta)
    dimnames(other_mod_counts) <- dimnames(beta)

    if (is.null(positions)) {
        positions <- seq_len(n_sites) * 1000L
    }
    stopifnot(length(positions) == n_sites)

    if (is.null(sample_info)) {
        sample_info <- data.frame(
            sample_name = sample_names,
            condition = sample_names,
            replicate = seq_along(sample_names),
            stringsAsFactors = FALSE
        )
    }

    mod_type <- rep(mod_type, length.out = n_sites)
    motif <- rep(motif, length.out = n_sites)
    chrom <- rep(chrom, length.out = n_sites)
    strand <- rep(strand, length.out = n_sites)

    site_gr <- GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(start = positions, width = 1L),
        strand = strand,
        mod_type = factor(mod_type, levels = c("4mC", "5mC", "6mA")),
        motif = motif
    )
    if (!is.null(site_metadata)) {
        site_metadata <- S4Vectors::DataFrame(site_metadata)
        stopifnot(nrow(site_metadata) == n_sites)
        GenomicRanges::mcols(site_gr) <- cbind(
            GenomicRanges::mcols(site_gr),
            site_metadata
        )
    }
    GenomeInfoDb::seqinfo(site_gr) <- GenomeInfoDb::Seqinfo(
        seqnames = unique(chrom),
        seqlengths = stats::setNames(rep(seqlength, length(unique(chrom))), unique(chrom)),
        isCircular = FALSE
    )

    cd <- S4Vectors::DataFrame(sample_info)
    rownames(cd) <- cd$sample_name

    rse <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            methylation = beta,
            coverage = coverage,
            mod_counts = mod_counts,
            canonical_counts = canonical_counts,
            other_mod_counts = other_mod_counts
        ),
        rowRanges = site_gr,
        colData = cd
    )
    obj <- new("commaData", rse)
    S4Vectors::metadata(obj)$assay_defaults <- list(
        methylation = "methylation",
        coverage = "coverage",
        mod_counts = "mod_counts",
        canonical_counts = "canonical_counts",
        other_mod_counts = "other_mod_counts"
    )
    S4Vectors::metadata(obj)$assay_provenance <- list(
        methylation = commaKit:::.makeAssayLayerRecord(
            type = "filtered_beta",
            source = "test_fixture",
            role = "methylation",
            parent_assays = "coverage",
            method = "test_fixture",
            default_for = "methylation"
        ),
        coverage = commaKit:::.makeAssayLayerRecord(
            type = "observed_total_coverage",
            source = "test_fixture",
            role = "coverage",
            method = "test_fixture",
            default_for = "coverage"
        ),
        mod_counts = commaKit:::.makeAssayLayerRecord(
            type = "reconstructed_counts",
            source = "test_fixture",
            role = "mod_counts",
            parent_assays = c("methylation", "coverage"),
            method = "round_beta_times_coverage",
            default_for = "mod_counts"
        ),
        canonical_counts = commaKit:::.makeAssayLayerRecord(
            type = "reconstructed_counts",
            source = "test_fixture",
            role = "canonical_counts",
            parent_assays = c("coverage", "mod_counts"),
            method = "coverage_minus_mod_counts",
            default_for = "canonical_counts"
        ),
        other_mod_counts = commaKit:::.makeAssayLayerRecord(
            type = "observed_counts",
            source = "test_fixture",
            role = "other_mod_counts",
            parent_assays = "coverage",
            method = "test_fixture_zero_other_mod_counts",
            default_for = "other_mod_counts"
        )
    )
    obj
}

.make_two_modtype_fixture <- function(n_6ma = 8L, n_5mc = 4L,
                                      sample_names = c("samp1", "samp2"),
                                      conditions = NULL, replicate = NULL,
                                      seed = 2L) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    n_sites <- n_6ma + n_5mc
    n_samples <- length(sample_names)

    beta <- matrix(
        runif(n_sites * n_samples, 0.1, 0.9),
        nrow = n_sites,
        dimnames = list(NULL, sample_names)
    )
    coverage <- matrix(
        20L,
        nrow = n_sites,
        ncol = n_samples,
        dimnames = dimnames(beta)
    )
    if (is.null(conditions)) {
        conditions <- rep(c("control", "treatment"), length.out = n_samples)
    }
    if (is.null(replicate)) {
        replicate <- seq_len(n_samples)
    }
    sample_info <- data.frame(
        sample_name = sample_names,
        condition = conditions,
        replicate = replicate,
        stringsAsFactors = FALSE
    )

    .make_commaData_fixture(
        beta = beta,
        coverage = coverage,
        sample_info = sample_info,
        positions = c(seq_len(n_6ma) * 100L, seq_len(n_5mc) * 200L),
        mod_type = c(rep("6mA", n_6ma), rep("5mC", n_5mc)),
        motif = c(rep("GATC", n_6ma), rep("CCWGG", n_5mc)),
        strand = c(rep("+", n_6ma), rep("-", n_5mc))
    )
}

.make_diff_methyl_fixture <- function(n_sites = 20L, n_ctrl = 2L, n_treat = 1L,
                                      seed = 99L) {
    set.seed(seed)
    n_samples <- n_ctrl + n_treat
    n_diff <- n_sites %/% 2L
    beta_ctrl <- c(rep(0.9, n_diff), rep(0.5, n_sites - n_diff)) +
        matrix(rnorm(n_ctrl * n_sites, 0, 0.05), nrow = n_sites)
    beta_treat <- c(rep(0.2, n_diff), rep(0.5, n_sites - n_diff)) +
        rnorm(n_treat * n_sites, 0, 0.05)

    beta <- cbind(
        matrix(pmax(0, pmin(1, beta_ctrl)), nrow = n_sites, ncol = n_ctrl),
        matrix(pmax(0, pmin(1, beta_treat)), nrow = n_sites, ncol = n_treat)
    )
    colnames(beta) <- c(
        paste0("ctrl_", seq_len(n_ctrl)),
        paste0("treat_", seq_len(n_treat))
    )
    coverage <- matrix(30L, nrow = n_sites, ncol = n_samples,
                       dimnames = list(NULL, colnames(beta)))
    sample_info <- data.frame(
        sample_name = colnames(beta),
        condition = c(rep("control", n_ctrl), rep("treatment", n_treat)),
        replicate = seq_len(n_samples),
        stringsAsFactors = FALSE
    )

    .make_commaData_fixture(
        beta = beta,
        coverage = coverage,
        sample_info = sample_info,
        positions = seq_len(n_sites) * 100L,
        site_metadata = data.frame(
            is_diff = c(rep(TRUE, n_diff), rep(FALSE, n_sites - n_diff))
        )
    )
}
