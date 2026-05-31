# Integration tests for the core comma workflow.
#
# These tests intentionally exercise multiple exported functions together. Unit
# tests verify each function in isolation; this file verifies that their output
# contracts still compose into the workflow users actually run.

# ─── Core workflow: annotate -> test -> extract -> filter ─────────────────────

.write_pipeline_modkit_files <- function(object) {
    site_df <- as.data.frame(siteInfo(object))
    beta <- methylation(object)
    cov <- siteCoverage(object)

    tmp_dir <- tempfile("comma-integration-modkit-")
    dir.create(tmp_dir)

    mod_codes <- c("6mA" = "a", "5mC" = "m", "4mC" = "21839")
    files <- character(ncol(beta))
    names(files) <- colnames(beta)

    for (sample_name in colnames(beta)) {
        sample_beta <- beta[, sample_name]
        sample_cov <- cov[, sample_name]
        keep <- !is.na(sample_beta) & !is.na(sample_cov)

        n_valid <- as.integer(sample_cov[keep])
        n_mod <- as.integer(round(sample_beta[keep] * n_valid))
        n_canonical <- pmax(n_valid - n_mod, 0L)
        mod_code <- mod_codes[as.character(site_df$mod_type[keep])]
        motif <- site_df$motif[keep]
        mod_code <- ifelse(
            is.na(motif),
            mod_code,
            paste(mod_code, motif, "1", sep = ",")
        )

        bed_df <- data.frame(
            chrom = site_df$chrom[keep],
            start = as.integer(site_df$position[keep]) - 1L,
            end = as.integer(site_df$position[keep]),
            mod_code = mod_code,
            score = as.integer(round(sample_beta[keep] * 1000)),
            strand = site_df$strand[keep],
            thickStart = as.integer(site_df$position[keep]) - 1L,
            thickEnd = as.integer(site_df$position[keep]),
            itemRgb = "255,0,0",
            Nvalid_cov = n_valid,
            fraction_modified = round(sample_beta[keep] * 100, 2),
            Nmod = n_mod,
            Ncanonical = n_canonical,
            Nother_mod = 0L,
            Ndelete = 0L,
            Nfail = 0L,
            Ndiff = 0L,
            Nnocall = 0L,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        if (sample_name != colnames(beta)[1L]) {
            bed_df <- bed_df[rev(seq_len(nrow(bed_df))), , drop = FALSE]
        }

        file <- file.path(tmp_dir, paste0(sample_name, ".bed"))
        write.table(
            bed_df,
            file = file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
        files[[sample_name]] <- file
    }

    files
}

test_that("core workflow returns expected object and result contracts", {
    data(comma_example_data)

    annotated <- annotateSites(comma_example_data, keep = "overlap")
    dm <- diffMethyl(
        annotated,
        formula = ~ condition,
        mod_type = "6mA",
        method = "quasi_f"
    )
    res <- results(dm, mod_type = "6mA")
    sig <- filterResults(dm, padj = 0.05, delta_beta = 0.2, mod_type = "6mA")

    expect_s4_class(annotated, "commaData")
    expect_s4_class(dm, "commaData")
    expect_s3_class(res, "data.frame")
    expect_s3_class(sig, "data.frame")

    n_6ma <- sum(siteInfo(comma_example_data)$mod_type == "6mA")
    expect_equal(nrow(res), n_6ma)

    expected_result_cols <- c(
        "chrom", "position", "strand", "mod_type", "motif", "mod_context",
        "feature_types", "feature_names",
        "dm_pvalue", "dm_padj", "dm_delta_beta"
    )
    expect_true(all(expected_result_cols %in% colnames(res)))

    expect_true(all(sig$mod_type == "6mA"))
    expect_true(all(sig$dm_padj <= 0.05))
    expect_true(all(abs(sig$dm_delta_beta) >= 0.2))
    expect_false(any(is.na(sig$dm_padj)))
    expect_false(any(is.na(sig$dm_delta_beta)))

    # filterResults() should only return rows from results(). Use the stable
    # genomic key rather than row names, because results() preserves site keys
    # but callers should not rely on row-name behavior.
    res_keys <- paste(res$chrom, res$position, res$strand,
                      res$mod_type, res$motif, sep = ":")
    sig_keys <- paste(sig$chrom, sig$position, sig$strand,
                      sig$mod_type, sig$motif, sep = ":")
    expect_true(all(sig_keys %in% res_keys))
})

test_that("core workflow recovers ground-truth differential 6mA sites", {
    data(comma_example_data)

    annotated <- annotateSites(comma_example_data, keep = "overlap")
    dm <- diffMethyl(
        annotated,
        formula = ~ condition,
        mod_type = "6mA",
        method = "quasi_f"
    )
    res <- results(dm, mod_type = "6mA")
    sig <- filterResults(dm, padj = 0.05, delta_beta = 0.2, mod_type = "6mA")

    expect_true("is_diff" %in% colnames(res))
    expect_equal(sum(res$is_diff, na.rm = TRUE), 30L)

    delta_diff <- res$dm_delta_beta[res$is_diff & !is.na(res$dm_delta_beta)]
    delta_nondiff <- res$dm_delta_beta[!res$is_diff & !is.na(res$dm_delta_beta)]

    # The synthetic ground truth sets 6mA differential sites high in control
    # and low in treatment, so treatment - control should be strongly negative.
    expect_lt(median(delta_diff), -0.3)
    expect_lt(abs(median(delta_nondiff)), 0.1)

    n_true_diff <- sum(res$is_diff, na.rm = TRUE)
    n_detected <- sum(res$is_diff & !is.na(res$dm_pvalue) & res$dm_pvalue < 0.2,
                      na.rm = TRUE)
    expect_gte(n_detected, floor(n_true_diff * 0.5))

    # The pipeline should produce significant hits for this synthetic dataset.
    # If filterResults() regresses to return zero rows, that is a test failure,
    # not a skip — the whole point of this integration test is to verify the
    # workflow produces non-vacuous filtered results.
    expect_gt(nrow(sig), 0L)
    expect_gt(sum(sig$is_diff, na.rm = TRUE), 0L)

    precision <- mean(sig$is_diff, na.rm = TRUE)
    expect_gte(precision, 0.5)
})

test_that("full pipeline composes from modkit files through enrichment", {
    data(comma_example_data)

    files <- .write_pipeline_modkit_files(comma_example_data)
    on.exit(unlink(dirname(files[[1L]]), recursive = TRUE), add = TRUE)
    col_data <- as.data.frame(sampleInfo(comma_example_data))
    gff_file <- system.file("extdata", "example.gff3", package = "commaKit")
    expect_true(file.exists(gff_file))

    object <- commaData(
        files = files,
        colData = col_data,
        genome = c(chr_sim = 100000L),
        annotation = gff_file,
        caller = "modkit"
    )

    annotated <- annotateSites(object, keep = "overlap")
    dm <- diffMethyl(annotated, formula = ~ condition, mod_type = "6mA")
    res <- results(dm, mod_type = "6mA")
    sig <- filterResults(dm, padj = 0.05, delta_beta = 0.2, mod_type = "6mA")

    truth <- as.data.frame(siteInfo(comma_example_data))
    truth <- truth[truth$mod_type == "6mA", , drop = FALSE]
    truth_key <- paste(truth$chrom, truth$position, truth$strand,
                       truth$mod_type, truth$motif, sep = ":")
    res_key <- paste(res$chrom, res$position, res$strand,
                     res$mod_type, res$motif, sep = ":")
    truth_is_diff <- truth$is_diff[match(res_key, truth_key)]

    expect_s4_class(object, "commaData")
    expect_s4_class(annotated, "commaData")
    expect_s4_class(dm, "commaData")
    expect_s3_class(res, "data.frame")
    expect_s3_class(sig, "data.frame")
    expect_equal(nrow(res), 393L)
    expect_false(anyNA(truth_is_diff))
    expect_equal(sum(truth_is_diff), 30L)
    expect_gt(nrow(sig), 0L)

    sig_key <- paste(sig$chrom, sig$position, sig$strand,
                     sig$mod_type, sig$motif, sep = ":")
    sig_truth_is_diff <- truth$is_diff[match(sig_key, truth_key)]
    expect_true(all(sig_key %in% res_key))
    expect_false(anyNA(sig_truth_is_diff))
    expect_equal(sum(sig_truth_is_diff), 30L)
    expect_equal(sum(sig_truth_is_diff) / sum(truth$is_diff), 1)
    expect_gte(mean(sig_truth_is_diff), 0.95)

    term2gene <- data.frame(
        term = c("PATH:01", "PATH:01", "PATH:02", "PATH:02", "PATH:03"),
        gene = c("geneA", "geneB", "geneC", "geneD", "geneE"),
        stringsAsFactors = FALSE
    )
    enr <- enrichMethylation(
        dm,
        method = "ora",
        TERM2GENE = term2gene,
        minGSSize = 1L
    )

    expect_type(enr, "list")
    expect_true(all(c("go", "kegg") %in% names(enr)))
    expect_s4_class(enr$go, "enrichResult")
    expect_null(enr$kegg)
})
