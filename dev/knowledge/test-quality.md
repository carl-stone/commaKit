# Test Quality — What We Know

**Last updated:** 2026-05-14
**Maintained by:** commaBot

---

## Strong Tests — High Confidence in Correctness

These functions have tests that verify **output is correct**, not just that the function runs.

| Function | Evidence | Confidence |
|----------|----------|------------|
| `diffMethyl()` | Ground-truth recovery (30 known diff sites out of 393 6mA), delta_beta sign checks, cross-backend correlation (limma vs quasi_f r>0.95), edge cases (zero coverage, perfect separation, single-condition sites), methylKit crash regression test | **High** |
| `annotateSites()` | Strand-aware rel_position signs verified, frac_position in [0,1] verified, multi-feature overlaps, list-column types checked (CharacterList/IntegerList), window boundary behavior | **High** |
| `mValues()` | Formula verified with known inputs: beta=0.8, cov=10, alpha=0.5 → M=log2(8.5/2.5), NA propagation, zero coverage → NA, alpha parameter validation | **High** |
| `writeBED()` | 0-based coordinates verified (position-1), score=round(beta*1000) verified, exact RGB values per band verified, NA exclusion, BED9 format (8 tabs) | **High** |

These tests use **known inputs with expected outputs**. A bug would fail the test.

---

## Weak Tests — Confidence in "Doesn't Crash"

These functions have smoke tests only. They verify the function returns the right **type**, but not that the **data is correct**.

| Function | Test Coverage | Risk |
|----------|---------------|------|
| `plot_coverage()` | Checks `ggplot` class, layer count | Silent column drop or wrong mapping would pass |
| `plot_methylation_distribution()` | Checks `ggplot` class, layer count | Silent column drop or wrong mapping would pass |
| `plot_pca()` | Checks `ggplot` class, layer count | Silent column drop or wrong mapping would pass |
| `plot_genome_track()` | Checks `ggplot` class, layer count | Silent column drop or wrong mapping would pass |
| `plot_metagene()` | Checks `ggplot` class, layer count | Silent column drop or wrong mapping would pass |
| `plot_tss_profile()` | Checks `ggplot` class, layer count | Silent column drop or wrong mapping would pass |
| `plot_volcano()` | Checks `ggplot` class, layer count | Silent column drop or wrong mapping would pass |
| `plot_heatmap()` | Checks `Heatmap` class | Silent column drop or wrong mapping would pass |
| `slidingWindow()` | Checks output structure, value ranges | Circular wrap correctness not verified |
| `enrichMethylation()` | Uses fake TERM2GENE (smart for CI) | Real clusterProfiler integration untested |

**What this means:** If `plot_volcano()` silently dropped all points with `padj < 0.05`, the test would still pass. The plot would just look empty.

---

## Documentation Quality — In Sync with Code

Verified 2026-05-14 by comparing `R/*.R` source to `man/*.Rd`:

- All 34 exported functions have complete `@param`, `@return`, `@examples`
- No stubs, no TODOs, no "document me" placeholders
- Examples actually run (verified by execution)
- Return values match what the docs claim

**Assessment:** Better than most R packages. Claude Code did competent work here.

**Remaining issue:** Parameter type inconsistency is documented correctly but still a UX friction point. See `known-issues.md`.

---

## What We Haven't Verified

1. **Plot data mappings.** We haven't inspected `ggplot_build(p)$data[[1]]` to verify that `plot_volcano()` actually maps `dm_padj` to the y-axis and `dm_delta_beta` to the x-axis. This would require adding data verification tests.

2. **slidingWindow circular correctness.** The test checks that `circular=TRUE` differs from `circular=FALSE`, but doesn't verify the actual smoothed values at chromosome boundaries. For a site at position 99500 on a 100kb chromosome with window=1000, is the smoothed value mathematically correct?

3. **Integration correctness.** No test runs the full pipeline: `commaData() → annotateSites() → diffMethyl() → results() → filterResults() → enrichMethylation()`. A breaking change in one function's output format could silently break downstream consumers.

4. **Parser edge cases.** Tests use `example_modkit.bed` and `comma_example_data`. We haven't verified that parsers handle real production data edge cases (malformed lines, unexpected chromosomes, missing fields).

5. **Small-sample behavior.** `diffMethyl()` with exactly 2 samples per condition has 1 residual df for quasi-F. Is this still valid? What happens with exactly 2 samples total? Not tested.

6. **Performance limits.** How many sites can `diffMethyl()` handle? What's the memory footprint for 50K sites × 6 samples? Unknown.

---

## Assumptions Made During Audits

1. **Test suite is representative.** We read ~8 of ~20 test files. We assumed the unread plot tests follow the same smoke-test pattern as `plot_volcano.R`. Spot-checking suggests this is true.

2. **Ground truth is correct.** Tests assert that 30 of 393 6mA sites are differentially methylated. We assumed `set.seed(1312)` in `data-raw/create_example_data.R` produces what it claims. We didn't audit that script.

3. **No test mocking hides bugs.** We assumed tests aren't mocking or patching functions in ways that hide real bugs. Spot-checking suggests this is true.

4. **Examples in docs run.** Some examples are marked `\donttest` or `\dontrun`, which means even R CMD check doesn't verify them. We tested a sample manually; they worked.

---

## How to Improve Confidence

To move a function from "weak tests" to "strong tests":

1. **Add data verification.** For plots, inspect `ggplot_build(p)$data[[1]]` and check that key values match input.
2. **Add integration test.** Run full pipeline on `comma_example_data`, verify the 30 ground-truth diff sites are recovered.
3. **Add boundary tests.** For `slidingWindow`, test specific positions near chromosome ends with known expected values.

See `BACKLOG.md` for prioritized work items.
