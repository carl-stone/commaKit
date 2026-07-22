---
type: Issue Knowledge Base
title: Known Issues
description: Durable bugs, gotchas, edge cases, and resolved issue notes for commaKit.
resource: dev/knowledge/known-issues.md
tags: [bugs, gotchas, edge-cases, tests]
timestamp: 2026-06-15T00:00:00Z
status: current
owner: Carl Stone
---

# Known Issues — Bugs, Gotchas, Edge Cases

**Last updated:** 2026-06-23
**Maintained by:** commaBot

**Currentness note:** This file preserves useful historical gotchas. For live
API shape and exported functions, prefer [Architecture](architecture.md) and
the current source.

---

## Confirmed Bugs

None confirmed at this time. All known issues are potential, documented behavior, or edge cases.

---

## Potential Bugs (Verified)

### P-001: `diag()` scalar trap — VERIFIED SAFE

**Risk:** If `diag(x)` appears in the source with scalar `x`, it creates an x×x identity matrix instead of a 1×1 matrix. Silent bug.

**Status:** Verified 2026-05-30. `grep -r "diag(" R/*.R` returns no matches. No `diag()` usage in the codebase.

**Action:** None needed. Documented as a convention gotcha for future contributors.

---

### P-002: `S4Vectors::rename()` vs `dplyr::rename()` masking — VERIFIED SAFE

**Risk:** If the code uses `rename()` without `dplyr::` prefix, it may call `S4Vectors::rename()` which has different semantics.

**Status:** Verified 2026-05-30. All `rename()` calls in the codebase use the `dplyr::rename()` prefix explicitly.

**Action:** None needed. Documented as a convention gotcha for future contributors.

---

## Known Gotchas (Documented Behavior That Trips Users)

### G-001: methylKit crashes on zero-coverage sites

**What happens:** `methylKit::calculateDiffMeth()` crashes when a site has zero coverage in all samples after filtering.

**Fix:** commaKit wraps this and assigns `p = 1` (consistent with null hypothesis). Regression test exists in `test-diffMethyl.R`.

**User impact:** Shouldn't crash anymore, but may be surprising that such sites get `padj = 1` instead of `NA`.

---

### G-002: CI pinned to R 4.5

**Why:** S4Vectors C API breaks in R 4.6.0. Bioconductor hasn't patched yet.

**User impact:** None for users. Developers must use R 4.5 for CI.

---

### G-003: `org.EcK12.eg.db` requires `::` syntax in examples

**What:** Can't use `library(org.EcK12.eg.db)` in examples. Must use `org.EcK12.eg.db::org.EcK12.eg.db`.

**Why:** R CMD check doesn't attach suggested packages during example evaluation.

**User impact:** Examples work, but users may be confused by the unusual syntax.

---

### G-004: Non-ASCII characters cause R CMD check notes

**What:** Characters like `×` in documentation cause encoding notes.

**Fix:** Replaced with ASCII equivalents (e.g., `x`).

**User impact:** None. Cosmetic only.

---

### G-005: `mod_type` parameter type varied by function — RESOLVED/STALE

**What:** Older docs recorded that some functions accepted a character vector
while others accepted only a single string.

**Current state:** The main exported filtering surfaces now document
`mod_type` as a character vector or `NULL`, routed through shared site-filter
helpers.

**Action:** Do not copy this older inconsistency into new docs. Check the
specific roxygen block before changing an API.

---

### G-006: `\donttest{}` examples that reference external files — FIXED

**What:** Four roxygen examples used `\donttest{}` but referenced files not included in the package. CI runs `--run-donttest`, so these examples executed and crashed during R CMD check.

**Affected functions:** `commaData()`, `loadAnnotation()`, `findMotifSites()`

**Fix:** Changed `\donttest{}` to `\dontrun{}` for all four. `\dontrun{}` means "show in docs but never execute," which is correct for examples that need user-provided files.

**Lesson:** `\donttest{}` is not safe for examples that need external files. Use `\dontrun{}` instead.

---

### G-007: `gene_role` default uses match-arg pattern

**What:** The default is `gene_role = c("target", "regulator", "both")` with `match.arg()`.

**Why:** R idiom for restricting to specific values.

**User impact:** Slightly different from typical "default = 'target'" pattern. No actual issue, just unusual.

---

## Edge Cases (Tested But Worth Knowing)

### E-001: Zero-variance sites in `diffMethyl()`

**What happens:** Sites with identical methylation across all samples get `padj = 1` or `NA` depending on backend.

**Tested:** Yes, edge case tests exist.

---

### E-002: Perfect separation in `diffMethyl()`

**What happens:** Sites where one group is all 0 and the other is all 1.

**Tested:** Yes, handled by methylKit wrapper.

---

### E-003: Single-condition sites

**What happens:** Sites that only appear in one condition (all samples in the other condition have `NA`).

**Tested:** Yes, edge case tests exist.

---

### E-004: Circular chromosome wrap in `slidingWindow()` — VERIFIED

**What happens:** Positions near chromosome boundaries wrap around when `circular = TRUE`.

**Tested:** Yes. Boundary behavior now has known-value tests for circular wrap-around.

**Confidence:** High for the tested bacterial-genome use cases.

---

## Missing Tests (Unknown if Bug)

### M-001: Full-pipeline integration — COVERED

**What:** Older audit notes recorded that no test ran `commaData() -> annotateSites() -> diffMethyl() -> results() -> filterResults() -> enrichMethylation()`.

**Current state:** `tests/testthat/test-integration.R` now exercises the end-to-end workflow. Keep strengthening this file when result-layer or enrichment contracts change.

**Confidence:** Medium. The pipeline shape is covered, but real biological data remain out of scope for the bundled synthetic fixture.

---

### M-002: Real-world parser edge cases — PARTIALLY COVERED

**What:** Real production data may have malformed lines, unexpected chromosomes,
missing fields, and non-standard delimiter shapes.

**Current state:** The modkit parser now uses an explicit tab delimiter
(`sep = "\t"`) so blank required fields are detected rather than silently
shifting later bedMethyl values. Beta is computed from authoritative count
fields (`Nmod / Nvalid_cov`) instead of the derived `fraction_modified`
percentage, and zero-valid-coverage rows with undefined fractions are dropped
by the coverage filter. Tests cover blank-field rejection, space-separated file
rejection, count-based beta, zero-coverage drops, and malformed partial rows.

**Remaining gaps:**
- No tests with real production modkit files beyond the bundled synthetic
  fixture.
- No coverage for unexpected chromosome names that do not match a user-supplied
  genome (the constructor catches mismatches, but parser-level behaviour is not
  tested with exotic names).
- No tests for very large bedMethyl files (performance and memory).

**Confidence:** Medium. The pipeline shape and specific edge cases above are
covered, but real biological data and performance at scale remain out of scope.

---

### M-003: Performance at scale

**What:** No tests for memory/runtime with 50K+ sites.

**Risk:** Package may not scale to real bacterial genomes.

**Confidence:** Unknown. Not tested.

---

## How to Use This Document

- **If you hit a crash:** Check "Confirmed Bugs" first, then "Known Gotchas"
- **If you see unexpected behavior:** May be in "Edge Cases"
- **If adding new tests:** Focus on "Missing Tests" section
- **If auditing code:** Start with "Potential Bugs" section
