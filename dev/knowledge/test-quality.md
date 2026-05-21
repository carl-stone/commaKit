# Test Quality Audit

**Last updated:** 2026-05-21 (full audit, issue #124)
**Maintained by:** commaBot

## Audit Method

Each `test_that()` block classified as:
- **contract** — verifies behavioral contract (correct output for known input, correct error for bad input)
- **smoke** — only verifies the function runs without error (low value)
- **edge-case** — verifies behavior at boundaries (zero-length, single element, NA, circular)
- **integration** — verifies multiple components work together

Verdicts: **keep** (genuinely useful), **remove** (uninformative/redundant), **strengthen** (smoke that should be contract)

---

## Summary

| File | Tests | Contract | Smoke | Edge-case | Integration | Keep | Remove | Strengthen |
|------|-------|----------|-------|-----------|-------------|------|--------|------------|
| test-commaData.R | 30 | 25 | 0 | 3 | 2 | 30 | 0 | 0 |
| test-accessors.R | 57 | 49 | 1 | 7 | 2 | 56 | 0 | 1 |
| test-annotateSites.R | 35 | 31 | 0 | 4 | 0 | 35 | 0 | 0 |
| test-diffMethyl.R | 54 | 44 | 3 | 6 | 2 | 50 | 3 | 2 |
| test-results.R | 22 | 20 | 0 | 2 | 0 | 22 | 0 | 0 |
| test-mValues.R | 21 | 18 | 2 | 2 | 0 | 20 | 0 | 1 |
| test-writeBED.R | 20 | 17 | 0 | 2 | 0 | 20 | 0 | 0 |
| test-plot_coverage.R | 13 | 6 | 5 | 1 | 1 | 6 | 2 | 5 |
| test-plot_genome_track.R | 12 | 6 | 4 | 1 | 1 | 7 | 0 | 5 |
| test-plot_heatmap.R | 10 | 5 | 1 | 3 | 1 | 6 | 0 | 4 |
| test-plot_metagene.R | 10 | 5 | 3 | 2 | 1 | 6 | 0 | 4 |
| test-plot_methylation_distribution.R | 11 | 5 | 3 | 1 | 2 | 5 | 2 | 4 |
| test-plot_pca.R | 18 | 13 | 4 | 1 | 1 | 14 | 0 | 4 |
| test-plot_tss_profile.R | 25 | 12 | 5 | 7 | 2 | 18 | 2 | 5 |
| test-plot_volcano.R | 18 | 10 | 2 | 4 | 2 | 12 | 0 | 4 |
| test-slidingWindow.R | 18 | 16 | 1 | 0 | 0 | 17 | 0 | 1 |
| test-enrichMethylation.R | 62 | 45 | 6 | 8 | 0 | 61 | 1 | 5 |
| test-methylomeSummary.R | 18 | 16 | 0 | 1 | 0 | 18 | 0 | 0 |
| test-coverageDepth.R | 11 | 9 | 1 | 0 | 0 | 10 | 0 | 1 |
| test-varianceByDepth.R | 10 | 8 | 0 | 1 | 0 | 10 | 0 | 0 |
| test-integration.R | 2 | 0 | 0 | 0 | 2 | 2 | 0 | 0 |
| test-parse_dorado.R | 22 | 14 | 0 | 8 | 0 | 22 | 0 | 0 |
| test-parse_megalodon.R | 14 | 10 | 0 | 3 | 0 | 14 | 0 | 0 |
| test-parsers.R | 21 | 17 | 1 | 2 | 0 | 19 | 2 | 0 |
| test-findMotifSites.R | 19 | 17 | 0 | 1 | 0 | 18 | 1 | 0 |
| test-loadAnnotation.R | 18 | 16 | 0 | 2 | 0 | 18 | 0 | 0 |
| test-buildKEGGTermGene.R | 42 | 33 | 0 | 5 | 3 | 42 | 0 | 0 |
| test-genome_utils.R | 24 | 15 | 0 | 8 | 0 | 24 | 0 | 0 |
| **TOTAL** | **637** | **486** | **43** | **74** | **25** | **608** | **13** | **46** |

**Key findings:**
- 76% of tests are contract tests (486/637) — strong foundation
- 7% are smoke tests (43/637) — the main weakness, concentrated in plot functions
- 13 tests to remove (redundant)
- 46 tests to strengthen (smoke → contract)
- Plot tests are the weakest area: 27/117 (23%) are smoke, 35 need strengthening

---

## Tests to Remove (13)

| File | Test name | Reason |
|------|-----------|--------|
| test-diffMethyl.R | method='quasi_f' returns commaData with correct columns | Redundant with tests 1-4 that already cover quasi_f (the default) |
| test-diffMethyl.R | method='quasi_f' produces valid p-values in [0, 1] | Redundant with test 8 |
| test-diffMethyl.R | quasi_f and limma delta_beta are highly correlated | Duplicate of test 27 (same correlation, swapped order) |
| test-plot_coverage.R | per_sample = FALSE returns ggplot | Redundant with base test + unfaceted contract test |
| test-plot_methylation_distribution.R | per_sample = FALSE returns ggplot | Redundant with base test + single-mod contract test |
| test-plot_methylation_distribution.R | example data, mod_type = '6mA' | Redundant with example-data test; only adds mod_type filter |
| test-plot_tss_profile.R | works with comma_example_data, window = 1000L | Redundant with example-data test; window already covered by contract test |
| test-enrichMethylation.R | data.frame input runs without error when valid | Redundant with test above that exercises data.frame + checks warning |
| test-parsers.R | .parseDorado() stops with informative error for missing file | Duplicate of test-parse_dorado.R |
| test-parsers.R | .parseDorado() stops with informative error for non-character file | Duplicate of test-parse_dorado.R |
| test-findMotifSites.R | DNAStringSet multi-chromosome searches all sequences | Redundant with FASTA multi-chromosome test; only input type differs |

---

## Tests to Strengthen (46)

### Plot tests (35)

**Pattern:** Most plot tests only check `expect_s3_class(p, "ggplot")`. They should use `ggplot_build()` to verify data mappings.

| File | Test name | What to verify |
|------|-----------|----------------|
| test-plot_coverage.R | returns ggplot for valid input | Verify plot data contains coverage values |
| test-plot_coverage.R | mod_type filter accepted | Verify filtered data excludes non-matching mod_type |
| test-plot_coverage.R | works with a single-sample object | Verify single panel, correct sample name |
| test-plot_coverage.R | mod_type accepts character vector | Verify both mod types appear in plot data |
| test-plot_genome_track.R | returns ggplot for valid chromosome | Verify track contains expected point data |
| test-plot_genome_track.R | annotation = FALSE suppresses annotation track | Verify annotation layer absent (fewer layers) |
| test-plot_genome_track.R | mod_type filter returns ggplot without error | Verify filtered data excludes non-matching mod types |
| test-plot_genome_track.R | start only (no end) accepted | Verify data starts at specified position |
| test-plot_heatmap.R | returns ggplot for valid input | Verify tile count matches n_sites |
| test-plot_heatmap.R | n_sites larger than available sites clamps silently | Verify all available sites shown |
| test-plot_heatmap.R | NA beta values handled without error | Verify NA cells render as grey/missing |
| test-plot_metagene.R | returns ggplot for valid feature type | Verify metagene data points present |
| test-plot_metagene.R | mod_type filter accepted | Verify filtering reduces plotted data |
| test-plot_metagene.R | n_bins parameter accepted | Verify bin count affects x-axis breaks |
| test-plot_methylation_distribution.R | returns ggplot for valid input | Verify density/histogram layer exists |
| test-plot_methylation_distribution.R | mod_type filter returns ggplot | Verify filtered data excludes other mod type |
| test-plot_methylation_distribution.R | NAs in beta values are silently excluded | Verify fewer points in ggplot_build data |
| test-plot_pca.R | returns ggplot for valid input | Verify scatter layer has n sample points |
| test-plot_pca.R | color_by argument accepted | Verify colour aesthetic mapping present |
| test-plot_pca.R | shape_by = NULL accepted | Verify no shape aesthetic mapped |
| test-plot_pca.R | mod_type filter reduces sites used | Verify PCA computed on fewer sites |
| test-plot_tss_profile.R | returns ggplot for valid input | Verify TSS-relative positions plotted |
| test-plot_tss_profile.R | returns ggplot when mod_type is specified | Verify filtering applied |
| test-plot_tss_profile.R | color_by = 'mod_type' returns ggplot | Verify colour aesthetic maps to mod_type |
| test-plot_tss_profile.R | color_by = 'regulatory_element' with valid types | Verify regulatory coloring applied |
| test-plot_tss_profile.R | color_by = 'none' + facet_by = 'mod_type' + show_smooth | Verify interaction (no colour, faceted, smooth) |
| test-plot_volcano.R | returns ggplot for valid input | Verify point layer has expected row count |
| test-plot_volcano.R | custom thresholds accepted | Verify threshold lines reflect custom values |
| test-plot_volcano.R | rows with NA padj are excluded without error | Verify fewer points plotted |
| test-plot_volcano.R | rows with NA delta_beta handled without error | Verify fewer points plotted |

### Non-plot tests (11)

| File | Test name | What to verify |
|------|-----------|----------------|
| test-accessors.R | coverage() returns a matrix | Add `is.numeric()` or `is.integer()` check |
| test-diffMethyl.R | p_adjust_method = 'bonferroni' is accepted | Verify bonferroni padj >= pvalue |
| test-diffMethyl.R | works with comma_example_data | Verify return type or column presence |
| test-mValues.R | alpha = 1 runs without error | Verify output matrix is valid (finite, correct dims) |
| test-slidingWindow.R | circular=FALSE works without error | Verify boundary behavior or compare vs circular=TRUE |
| test-enrichMethylation.R | ORA with TERM2NAME works without error | Verify TERM2NAME appears in enrichResult |
| test-enrichMethylation.R | mod_type filter passes through to results() | Verify filtering changes gene set |
| test-enrichMethylation.R | feature_type = 'gene' runs without error | Verify gene-filtered results differ |
| test-enrichMethylation.R | feature_type = NULL includes all features | Verify NULL equals default |
| test-enrichMethylation.R | gene_role='target' uses target genes | Verify target genes actually used |
| test-enrichMethylation.R | gene_role='regulator' uses regulator genes | Verify regulator genes actually used |
| test-coverageDepth.R | method='median' works without error | Verify median values differ from mean |

---

## Per-File Audit

### test-commaData.R (30 tests: 25 contract, 3 edge-case, 2 integration)

All 30 keep. Strong contract coverage for S4 validity, constructor errors, mod_type/mod_context filtering, and on-demand computation. No smoke tests.

### test-accessors.R (57 tests: 49 contract, 1 smoke, 7 edge-case, 2 integration)

56 keep, 1 strengthen. Nearly all contract tests. One smoke test: `coverage() returns a matrix` should also verify storage mode.

### test-annotateSites.R (35 tests: 31 contract, 4 edge-case)

All 35 keep. Excellent contract coverage for all four keep modes (all/overlap/proximity/metagene), strand-aware positions, list-column types, and error contracts.

### test-diffMethyl.R (54 tests: 44 contract, 3 smoke, 6 edge-case, 2 integration)

50 keep, 3 remove, 2 strengthen. Strong ground-truth recovery tests. Three redundant quasi_f tests to remove. Two smoke tests to strengthen (bonferroni, comma_example_data).

### test-results.R (22 tests: 20 contract, 2 edge-case)

All 22 keep. Clean contract coverage for return structure, filtering, thresholds, and error contracts.

### test-mValues.R (21 tests: 18 contract, 2 smoke, 2 edge-case)

20 keep, 1 strengthen. Good formula verification with hand-computed values. One smoke test (alpha=1) to strengthen.

### test-writeBED.R (20 tests: 17 contract, 2 edge-case)

All 20 keep. Excellent — verifies exact BED format (0-based coords, score formula, RGB values), NA exclusion, and error contracts.

### test-plot_coverage.R (13 tests: 6 contract, 5 smoke, 1 edge-case, 1 integration)

6 keep, 2 remove, 5 strengthen. Faceting and axis-label contracts are good. Five smoke tests need data verification via ggplot_build.

### test-plot_genome_track.R (12 tests: 6 contract, 4 smoke, 1 edge-case, 1 integration)

7 keep, 5 strengthen. Start/end filtering contract is good. Four smoke tests need data verification.

### test-plot_heatmap.R (10 tests: 5 contract, 1 smoke, 3 edge-case, 1 integration)

6 keep, 4 strengthen. n_sites contract is good. Three edge-case tests (clamping, NA) only check class — should verify data.

### test-plot_metagene.R (10 tests: 5 contract, 3 smoke, 2 edge-case, 1 integration)

6 keep, 4 strengthen. x-axis range contract is good. Three smoke tests need data verification.

### test-plot_methylation_distribution.R (11 tests: 5 contract, 3 smoke, 1 edge-case, 2 integration)

5 keep, 2 remove, 4 strengthen. Faceting contracts are good. Two redundant parameter-variation tests to remove.

### test-plot_pca.R (18 tests: 13 contract, 4 smoke, 1 edge-case, 1 integration)

14 keep, 4 strengthen. return_data path is well-tested (5 contract tests). Four smoke tests need data verification.

### test-plot_tss_profile.R (25 tests: 12 contract, 5 smoke, 7 edge-case, 2 integration)

18 keep, 2 remove, 5 strengthen. Good contract coverage for x-axis range, vline, faceting, smooth layers. Five smoke tests need data verification.

### test-plot_volcano.R (18 tests: 10 contract, 2 smoke, 4 edge-case, 2 integration)

12 keep, 4 strengthen. Layer-count contracts and faceting contracts are good. Two edge-case NA tests need data verification.

### test-slidingWindow.R (18 tests: 16 contract, 1 smoke)

17 keep, 1 strengthen. Strong contract coverage including known-value test. One smoke test (circular=FALSE) needs boundary verification.

### test-enrichMethylation.R (62 tests: 45 contract, 6 smoke, 8 edge-case)

61 keep, 1 remove, 5 strengthen. Excellent internal function coverage (.parseTargetGenes, .extractGeneRoles). Five smoke tests need verification.

### test-methylomeSummary.R (18 tests: 16 contract, 1 edge-case)

All 18 keep. Strong contract coverage with exact value verification.

### test-coverageDepth.R (11 tests: 9 contract, 1 smoke)

10 keep, 1 strengthen. Good contract coverage for log2 transform. One smoke test (method='median') needs verification.

### test-varianceByDepth.R (10 tests: 8 contract, 1 edge-case)

All 10 keep. Good contract coverage.

### test-integration.R (2 tests: 2 integration)

Both keep. End-to-end pipeline tests.

### test-parse_dorado.R (22 tests: 14 contract, 8 edge-case)

All 22 keep. Strong contract and edge-case coverage for Dorado parser.

### test-parse_megalodon.R (14 tests: 10 contract, 3 edge-case)

All 14 keep. Good contract coverage for Megalodon parser.

### test-parsers.R (21 tests: 17 contract, 1 smoke, 2 edge-case)

19 keep, 2 remove. Two duplicate error tests already covered in test-parse_dorado.R.

### test-findMotifSites.R (19 tests: 17 contract, 1 edge-case)

18 keep, 1 remove. One redundant multi-chromosome test to remove.

### test-loadAnnotation.R (18 tests: 16 contract, 2 edge-case)

All 18 keep. Strong contract coverage for GFF3 and BED parsing.

### test-buildKEGGTermGene.R (42 tests: 33 contract, 5 edge-case, 3 integration)

All 42 keep. Excellent coverage of KEGG API, caching, ID mapping, and integration with enrichMethylation.

### test-genome_utils.R (24 tests: 15 contract, 8 edge-case)

All 24 keep. Strong coverage of .validateGenomeInfo, .circularIndex, and .makeSeqinfo.

---

## What We Haven't Verified

1. **Plot data mappings.** Most plot tests don't inspect `ggplot_build(p)$data[[1]]` to verify that the right columns map to the right aesthetics. A silent column drop or wrong mapping would pass all current tests.

2. **slidingWindow circular correctness.** The test checks that `circular=TRUE` differs from `circular=FALSE`, but doesn't verify the actual smoothed values at chromosome boundaries. For a site at position 99500 on a 100kb chromosome with window=1000, is the smoothed value mathematically correct?

3. **Integration correctness.** Only 2 integration tests exist. No test runs the full pipeline: `commaData() -> annotateSites() -> diffMethyl() -> results() -> filterResults() -> enrichMethylation()`. A breaking change in one function's output format could silently break downstream consumers.

4. **Parser edge cases.** Tests use `example_modkit.bed` and `comma_example_data`. We haven't verified that parsers handle real production data edge cases (malformed lines, unexpected chromosomes, missing fields).

5. **Small-sample behavior.** `diffMethyl()` with exactly 2 samples per condition has 1 residual df for quasi-F. Is this still valid? What happens with exactly 2 samples total? Not tested.

6. **Performance limits.** How many sites can `diffMethyl()` handle? What's the memory footprint for 50K sites x 6 samples? Unknown.

---

## How to Improve Confidence

1. **Strengthen plot tests** (35 tests): Use `ggplot_build()` to verify data mappings, point counts, aesthetic mappings, and layer structure.
2. **Remove redundant tests** (13 tests): Clean up duplicates that add no coverage.
3. **Add integration test**: Full pipeline on `comma_example_data`, verify the 30 ground-truth diff sites are recovered end-to-end.
4. **Add slidingWindow circular boundary test**: Known-value test at chromosome edges.
5. **Strengthen enrichment smoke tests** (5 tests): Verify filtering actually changes gene sets.

See GitHub Issues for prioritized work items.
