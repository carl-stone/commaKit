# Integration tests for the core comma workflow.
#
# These tests intentionally exercise multiple exported functions together. Unit
# tests verify each function in isolation; this file verifies that their output
# contracts still compose into the workflow users actually run.

# ─── Core workflow: annotate -> test -> extract -> filter ─────────────────────

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

    if (nrow(sig) > 0L) {
        precision <- mean(sig$is_diff, na.rm = TRUE)
        expect_gte(precision, 0.5)
    } else {
        skip("No significant sites found; cannot evaluate precision.")
    }
})
