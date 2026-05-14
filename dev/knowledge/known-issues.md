# Known Issues — Bugs, Gotchas, Edge Cases

**Last updated:** 2026-05-14
**Maintained by:** commaBot

---

## Confirmed Bugs

None confirmed at this time. All known issues are potential or unverified.

---

## Potential Bugs (Unverified)

### P-001: `diag()` scalar trap — status unknown

**Risk:** If `diag(x)` appears in the source with scalar `x`, it creates an x×x identity matrix instead of a 1×1 matrix. Silent bug.

**Status:** ROADMAP claims "no known usage" but commaBot hasn't verified this by grepping the source.

**Action needed:** `grep -r "diag(" R/*.R` and verify no scalar usage exists.

---

### P-002: `S4Vectors::rename()` vs `dplyr::rename()` masking — status unknown

**Risk:** If the code uses `rename()` without `dplyr::` prefix, it may call `S4Vectors::rename()` which has different semantics.

**Status:** ROADMAP claims "no known usage" but commaBot hasn't verified this.

**Action needed:** `grep -r "rename(" R/*.R` and verify all calls are `dplyr::rename()`.

---

## Known Gotchas (Documented Behavior That Trips Users)

### G-001: methylKit crashes on zero-coverage sites

**What happens:** `methylKit::calculateDiffMeth()` crashes when a site has zero coverage in all samples after filtering.

**Fix:** comma wraps this and assigns `p = 1` (consistent with null hypothesis). Regression test exists at `test-diffMethyl.R` lines 484-509.

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

### G-005: `mod_type` parameter type varies by function

**What:** Some functions accept a character vector, others accept only a single string.

| Function | `mod_type` type |
|----------|-----------------|
| `diffMethyl()` | Character vector |
| `subset()` | Character vector |
| `plot_coverage()` | Single string |
| `mValues()` | Single string |
| `methylomeSummary()` | Single string |

**Why:** Functions that split by type internally need to iterate over types. Functions that just filter don't.

**User impact:** User passes `c("6mA", "5mC")` to `plot_coverage()`, gets cryptic error.

**Mitigation:** Documented correctly in each function's `@param`. But still a UX friction point.

---

### G-006: `\donttest{}` examples that reference external files — FIXED

**What:** Four roxygen examples used `\donttest{}` but referenced files not included in the package. CI runs `--run-donttest`, so these examples executed and crashed during R CMD check.

**Affected functions:**
- `commaData()` — referenced `ctrl_1_modkit.bed`, `treat_1_modkit.bed`
- `loadAnnotation()` — referenced `my_genome.gff3`
- `findMotifSites()` — referenced `MG1655.fa` and `BSgenome.Ecoli.NCBI.20080805`

**Fix:** Changed `\donttest{}` to `\dontrun{}` for all four. `\dontrun{}` means "show in docs but never execute," which is correct for examples that need user-provided files.

**Fixed in:** PR #60, commit dd27a75.

**Lesson:** `\donttest{}` is not safe for examples that need external files. Use `\dontrun{}` instead. This was already documented in conventions gotcha #9 but was not in the backlog or known-issues until CI exposed it.

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

### E-004: Circular chromosome wrap in `slidingWindow()`

**What happens:** Positions near chromosome boundaries wrap around when `circular = TRUE`.

**Tested:** Partially. Test verifies difference from `circular = FALSE` but doesn't verify exact values.

**Confidence:** Medium. Logic looks correct but unverified at boundary positions.

---

## Missing Tests (Unknown if Bug)

### M-001: Integration across full pipeline

**What:** No test runs `commaData() → annotateSites() → diffMethyl() → results() → filterResults() → enrichMethylation()`.

**Risk:** Breaking change in one function's output format could silently break downstream consumers.

**Confidence:** Unknown. Individual functions tested, but integration is not.

---

### M-002: Real-world parser edge cases

**What:** Tests use `example_modkit.bed` only. Real production data may have malformed lines, unexpected chromosomes, missing fields.

**Risk:** Parser crashes on real data.

**Confidence:** Low for production use.

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
