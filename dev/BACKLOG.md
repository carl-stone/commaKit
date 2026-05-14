# comma Backlog — All Work Items

**Last updated:** 2026-05-14
**Maintained by:** commaBot

This is the single source of truth for everything that needs to be done. Carl reads this to see what's pending. commaBot reads this to know what to work on next.

---

## Priority Order

Items are sorted by priority. Top = most important.

| ID | Title | Priority | Status | Size | Source |
|----|-------|----------|--------|------|--------|
| **T-001** | Add integration test across full pipeline | 1 | proposed | M | Test audit |
| **P-001** | Add troubleshooting guide for data import | 2 | proposed | M | Package health audit |
| **P-002** | Unify `mod_type` parameter type across functions | 3 | proposed | M | Package health audit |
| **P-003** | Add method selection guidance for diffMethyl | 4 | proposed | S | Package health audit |
| **T-002** | Add data verification to plot tests | 5 | proposed | L | Test audit |
| **P-004** | Audit source for `diag()` scalar trap | 6 | proposed | S | Package health audit |
| **P-005** | Audit source for `rename()` masking | 7 | proposed | S | Package health audit |
| **T-003** | Verify slidingWindow circular correctness | 8 | proposed | M | Test audit |
| **T-004** | Test enrichMethylation with real clusterProfiler | 9 | proposed | M | Test audit |
| **P-006** | Add real-world example dataset | 10 | proposed | M | Package health audit |
| **P-007** | Add performance expectations documentation | 11 | proposed | S | Package health audit |
| **P-008** | Add "Getting Help" section to README | 12 | proposed | S | Package health audit |
| **T-005** | Test parser edge cases with real data | 13 | proposed | L | Test audit |
| **B-001** | Register Zenodo DOI | 14 | proposed | S | Pre-submission |
| **B-002** | Version bump to 0.99.0 | 15 | proposed | S | Pre-submission |
| **B-003** | Add NEWS.md 0.99.0 entry | 16 | proposed | S | Pre-submission |
| **B-004** | Submit to contributions.bioconductor.org | 17 | proposed | S | Pre-submission |
| **P-009** | Clean up `.codex/` directory | 18 | proposed | S | ROADMAP N5 |

---

## Work Item Details

### T-001: Add integration test across full pipeline

**Problem:** No test runs the full workflow end-to-end. A breaking change in one function's output format could silently break downstream consumers.

**Proposed solution:** Write one test that:
1. Takes `comma_example_data`
2. Runs: `annotateSites() → diffMethyl() → results() → filterResults()`
3. Verifies the 30 ground-truth diff sites are recovered

**Size:** M (medium) — moderate effort, clear scope

**Risk if skipped:** Unknown integration bugs could ship.

---

### P-001: Add troubleshooting guide for data import

**Problem:** New users will get stuck at data import. Chromosome name mismatches, malformed BED lines, missing MM/ML tags — all common issues. Error messages exist but don't guide users to solutions.

**Proposed solution:** Either:
- Add a `vignette("troubleshooting")` article, or
- Add a "Common Issues" section to README

Should cover:
- Chromosome name mismatches (`chr1` vs `1`)
- Malformed modkit BED lines
- Missing MM/ML tags in Dorado BAM
- Coverage threshold too aggressive → all sites NA

**Size:** M (medium) — requires writing, not coding

**Risk if skipped:** Users will email Carl instead of finding answers.

---

### P-002: Unify `mod_type` parameter type across functions

**Problem:** `mod_type` accepts a vector in some functions and a single string in others. Users pass `c("6mA", "5mC")` to `plot_coverage()` and get a cryptic error.

**Proposed solution:** Two options:
1. Make all functions accept vectors (larger change)
2. Document the distinction clearly and improve error messages (smaller change)

**Size:** M (medium) — requires API decision and implementation

**Risk if skipped:** Ongoing UX friction.

---

### P-003: Add method selection guidance for diffMethyl

**Problem:** `diffMethyl()` offers four backends (methylKit, limma, quasi_f, beta-binomial) but doesn't explain when to use which. Users default to methylKit without knowing if it's appropriate.

**Proposed solution:** Add a "Choosing a Statistical Backend" section to vignette or a separate article. Should explain:
- Small n (2 samples per group): which method?
- Zero-inflated data: which method?
- Complex designs: which method?

**Size:** S (small) — documentation only

**Risk if skipped:** Users pick wrong method.

---

### T-002: Add data verification to plot tests

**Problem:** Plot tests only verify that a ggplot object is returned. They don't verify that the plot data maps to the correct columns. A silent drop or wrong mapping would pass.

**Proposed solution:** For each plot function, add one test that:
1. Builds the plot: `p <- plot_volcano(res)`
2. Inspects the data: `d <- ggplot_build(p)$data[[1]]`
3. Verifies key values match input

**Size:** L (large) — 8 plot functions to update

**Risk if skipped:** Silent plot bugs could ship.

---

### P-004: Audit source for `diag()` scalar trap

**Problem:** If `diag(x)` appears with scalar `x`, it creates an x×x identity matrix instead of a 1×1 matrix. Silent bug. ROADMAP claims "no known usage" but unverified.

**Proposed solution:** Run `grep -r "diag(" R/*.R` and verify no scalar usage exists.

**Size:** S (small) — grep and verify

**Risk if skipped:** Potential silent bug.

---

### P-005: Audit source for `rename()` masking

**Problem:** If code uses `rename()` without `dplyr::` prefix, it may call `S4Vectors::rename()` which has different semantics. ROADMAP claims "no known usage" but unverified.

**Proposed solution:** Run `grep -r "rename(" R/*.R` and verify all calls are `dplyr::rename()`.

**Size:** S (small) — grep and verify

**Risk if skipped:** Potential silent bug.

---

### T-003: Verify slidingWindow circular correctness

**Problem:** Test verifies that `circular=TRUE` differs from `circular=FALSE` but doesn't verify the actual smoothed values at chromosome boundaries are mathematically correct.

**Proposed solution:** Add test with:
- Known chromosome size (e.g., 100kb)
- Known site positions near boundaries (e.g., position 500, position 99500)
- Known window size (e.g., 1000bp)
- Expected smoothed values computed by hand

**Size:** M (medium) — requires math verification

**Risk if skipped:** Circular logic may be wrong.

---

### T-004: Test enrichMethylation with real clusterProfiler

**Problem:** Tests use fake TERM2GENE to avoid network dependencies. Real clusterProfiler integration is untested beyond "returns an S4 object".

**Proposed solution:** Add test that:
- Runs ORA with real clusterProfiler on `comma_example_data`
- Skips on CI if no network or clusterProfiler unavailable
- Verifies output structure and at least one GO term returned

**Size:** M (medium) — requires test infrastructure

**Risk if skipped:** Enrichment may fail in production.

---

### P-006: Add real-world example dataset

**Problem:** `comma_example_data` is synthetic. Users can't verify their data "looks right" by comparison. No demonstration of realistic chromosome names, gene IDs, pathway names.

**Proposed solution:** Add a small real dataset (e.g., subset of E. coli K-12):
- ~1000 sites on actual E. coli genome
- Real gene symbols and GO terms
- Provenance documented

**Size:** M (medium) — requires data curation and licensing check

**Risk if skipped:** Lower user confidence.

---

### P-007: Add performance expectations documentation

**Problem:** Unknown limits. How many sites? What memory? What runtime?

**Proposed solution:** Add benchmark section to vignette or README:
- "Tested with 10K sites × 6 samples: ~30 seconds, ~500MB memory"
- "Not tested beyond 100K sites"
- Guidance on when to subsample

**Size:** S (small) — documentation only

**Risk if skipped:** Users may attempt unsupported scales.

---

### P-008: Add "Getting Help" section to README

**Problem:** README doesn't tell users where to go for help.

**Proposed solution:** Add section with:
- GitHub Issues link
- Bioconductor support site (when submitted)
- Email contact (optional)

**Size:** S (small) — documentation only

**Risk if skipped:** Users have no support channel.

---

### T-005: Test parser edge cases with real data

**Problem:** Tests use `example_modkit.bed` only. Real production data may have malformed lines, unexpected chromosomes, missing fields.

**Proposed solution:** Add test file with edge cases:
- Malformed lines (should warn and skip)
- Unknown chromosomes (should error or warn)
- Missing fields (should error)
- Mixed modification types (should handle)

**Size:** L (large) — requires test data creation and parser updates

**Risk if skipped:** Parser may crash on real data.

---

### B-001: Register Zenodo DOI

**Problem:** Bioconductor requires a DOI for the package.

**Proposed solution:** Register on Zenodo before submission.

**Size:** S (small) — administrative

**Dependencies:** None (can do anytime)

---

### B-002: Version bump to 0.99.0

**Problem:** Current version is 0.8.0.9000 (dev). Bioconductor requires 0.99.0 for new package submission.

**Proposed solution:** Update DESCRIPTION when ready to submit.

**Size:** S (small) — one-line change

**Dependencies:** All blockers resolved (already done)

---

### B-003: Add NEWS.md 0.99.0 entry

**Problem:** NEWS.md needs an entry for the Bioconductor submission version.

**Proposed solution:** Add 0.99.0 section summarizing changes since 0.8.0.

**Size:** S (small) — documentation

**Dependencies:** None

---

### B-004: Submit to contributions.bioconductor.org

**Problem:** Package is not yet submitted to Bioconductor.

**Proposed solution:** Submit after Carl registers on support site.

**Size:** S (small) — administrative

**Dependencies:** Carl must register on support.bioconductor.org

---

### P-009: Clean up `.codex/` directory

**Problem:** `.codex/` directory exists but is in .Rbuildignore. Not urgent but messy.

**Proposed solution:** Remove `.codex/` or document why it exists.

**Size:** S (small) — cleanup

**Dependencies:** None

---

## Status Definitions

- **proposed**: Not yet accepted by Carl
- **accepted**: Carl agrees it should be done, not yet started
- **in-progress**: commaBot is working on it
- **done**: Complete
- **wontfix**: Decided not to do

---

## How to Update This File

When adding new work items:
1. Add to the table with a new ID (increment the number)
2. Fill in all columns
3. Add a detail section below the table
4. Update STATUS.md to reflect the change

When changing status:
1. Update the table's Status column
2. Update STATUS.md to move the item between columns
3. If done, add to STATUS.md "Recently Completed" section
