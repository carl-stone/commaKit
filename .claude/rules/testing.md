---
paths:
  - "tests/testthat/**/*.R"
---

# Testing

## Framework

`testthat` edition 3. Tests live in `tests/testthat/`. Run with:

```r
devtools::test()
# or for a specific file:
testthat::test_file("tests/testthat/test-annotateSites.R")
```

There is **no** `tests/testthat/helper-fixtures.R` — fixtures are defined inline within each test file, or `comma_example_data` is loaded via `data(comma_example_data)`.

## Test Fixture: `comma_example_data`

Synthetic `commaData` created by `data-raw/create_example_data.R` (`set.seed(42)`):

- **300 sites**: 200 × 6mA, 100 × 5mC
- **3 samples**: ctrl_1, ctrl_2, treat_1
- **2 conditions**: control (n=2), treatment (n=1)
- **Genome**: chr_sim, 100 kb
- **Ground truth**: ~30 of 200 6mA sites are differentially methylated (control ~0.9, treatment ~0.25); marked in `rowData$is_diff`
- **Annotation**: 5 simulated genes (GRanges)

## Test File Coverage

> **Keep this table current.** Add a row whenever a new test file is created.

| File | Coverage |
|---|---|
| `test-commaData.R` | S4 class validity, constructor, bad inputs, show() — ~20 tests |
| `test-parsers.R` | Modkit column mapping, mod codes, coverage filter — ~15 tests |
| `test-accessors.R` | Matrix shape, value ranges, multi-mod-type, subsetting — ~20 tests |
| `test-genome_utils.R` | `.validateGenomeInfo()`, `.circularIndex()`, `.makeSeqinfo()` |
| `test-load_annotation.R` | GFF3/BED parsing, feature_type filtering, feature_subtype preservation from GFF3 attribute |
| `test-find_motif_sites.R` | Motif search, both strands, palindromic motifs |
| `test-parse_megalodon.R` | `.parseMegalodon()` aggregation, mod_type requirement |
| `test-annotateSites.R` | unified keep="all/overlap/proximity/metagene", strand-aware rel_position, frac_position NA for outside, metadata_cols parallel column — ~20 tests |
| `test-slidingWindow.R` | stat modes, circular wrap, genome-size inference — ~15 tests |
| `test-methylomeSummary.R` | per-sample stats, mod_type filtering, all-NA sample — ~11 tests |
| `test-coverageAnalysis.R` | `coverageDepth()` windowing, `varianceByDepth()` bins — ~8 tests |
| `test-writeBED.R` | file creation, track header, 9-col BED, RGB, NA exclusion, errors — ~20 tests |
| `test-diffMethyl.R` | `diffMethyl()` basic, statistical correctness, edge cases, ground-truth recovery — ~30 tests |
| `test-results.R` | `results()` and `filterResults()`: shape, filtering, thresholds, errors — ~23 tests |
| `test-parse_dorado.R` | `.cigarToRefPos()`, `.parseMmTag()` (ML boundary values), `.parseDorado()` errors — ~21 tests |
| `test-plot_distribution.R` | `plot_methylation_distribution()`: returns ggplot, mod_type filter, per_sample |
| `test-plot_genome_track.R` | `plot_genome_track()`: returns ggplot/patchwork, windowing, annotation |
| `test-plot_metagene.R` | `plot_metagene()`: returns ggplot, feature normalization |
| `test-plot_volcano.R` | `plot_volcano()`: returns ggplot, thresholds, coloring |
| `test-plot_heatmap.R` | `plot_heatmap()`: returns ggplot, top-N sites, sample annotation |
| `test-m_values.R` | `mValues()`: formula, NA/zero-coverage propagation, alpha, mod_type — ~24 tests |
| `test-plot_pca.R` | `plot_pca()`: returns ggplot, color_by/shape_by, return_data — ~22 tests |
| `test-plot_coverage.R` | `plot_coverage()`: returns ggplot, per_sample mode |
| `test-plot_tss_profile.R` | TSS window, color_by modes, regulatory_element fallback, loess, faceting — ~16 tests |
| `test-enrichment.R` | `.siteToGeneMap()`, `.computeGeneScores()`, `.parseTargetGenes()`, `.parseRegulatorGenes()`, `.extractGeneRoles()`, ORA+GSEA, gene_role, multi-feature_type, data.frame input — ~60 tests |

## Required Coverage

Every exported function needs tests for:
- Valid input → correct output (shape, type, values)
- Invalid input → informative error message (not a cryptic R error)
- Edge cases: empty data, NA handling, single sample, single site
