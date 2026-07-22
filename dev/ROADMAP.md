# commaKit Roadmap — Strategic Direction

**Last updated:** 2026-06-16
**Current version:** 0.2.0
**Public name:** commaKit (Comparative Microbial Methylomics Analysis Kit)
**R package namespace:** commaKit
**GitHub repo:** carl-stone/commaKit

This file is the strategic roadmap: where commaKit is going and why. GitHub Issues are the work log for bugs, blockers, and deferred readiness items.

---

## Strategic Goals

These are the reasons behind the work. They determine priority order.

1. **Correctness** — commaKit must produce results you can trust. If the stats are wrong, nothing else matters. This means real discriminating tests, not just smoke tests; known-quantity verification; and audits of known R gotchas.

2. **Usability** — Claire should be able to use commaKit independently. She needs docs, clear error messages, method selection guidance, and a package that doesn't surprise her with silent failures.

3. **Robustness** — commaKit should handle real data, not just a 588-site toy example. Edge cases, large genomes, malformed bedMethyl, production-scale site counts.

4. **Publishability** — Bioconductor-ready when the time comes. This is a low priority. We'd rather have a great, stable package installable from GitHub than a rushed Bioconductor submission.

---

## Milestone Sequence

Milestones are coherent groups of work that advance one or more strategic goals. They are ordered by priority. Each milestone has a GitHub Milestone for issue tracking.

### 1. Schema v2 — COMPLETE (correctness + robustness)

Restructured the `commaData` class for a stronger foundation. All 14 issues closed, merged to main as v0.2.0 on 2026-05-21.

**Issues:** #92–#97, #99, #105–#111
**GitHub Milestone:** commaData Schema v2
**Version:** 0.2.0

Key changes:
- SE → RangedSummarizedExperiment
- genomeInfo → Seqinfo
- annotation/motifSites → metadata()
- mod_context derived on demand, not stored
- mod_type as factor, min_coverage stored, site key convention
- No-rownames alignment (findOverlaps-based)

### 2. Test Quality — in progress (correctness)

The test suite has too many smoke tests and not enough discriminating tests. We don't really believe the tests yet — they verify "doesn't crash" but not "produces correct results."

**Open issues:** #75, #153, #233
**Completed issues:** #73, #74, #124, #126-#129
**Notable PRs:** #132, #133, #190, #204

What this milestone looks like when done:
- Remaining plot/fallback tests verify data mappings and optional-dependency behavior
- `enrichMethylation()` is exercised with real `clusterProfiler` plus real biological identifiers, not only synthetic TERM2GENE tables
- Integration test across full pipeline remains reliable and meaningful

Already completed:
- Full test audit and smoke-test cleanup
- Full-pipeline integration test
- slidingWindow circular correctness verified with known boundary values

See `dev/knowledge/test-quality.md` for the full audit.

### 3. Code Quality Audits — cleanup follow-through (correctness)

Thermonuclear review completed 2026-05-25. Findings filed as issues #135–#163, index #164.

**Merged PRs:**
- #165 — diffMethyl multi-level formula support (#135–#138)
- #166 — import/enrichment hardening (#139–#143, #146)
- #174, #177–#187, #192, #194, #197–#201, #203 — follow-up cleanup, plotting, docs, fixtures, and helper extraction

All behavior-oriented code-quality issues from the thermonuclear review have either been resolved or moved into explicit readiness buckets. The remaining repo-hygiene follow-through is #162 and should close when the buildignore/dependency cleanup lands.

### 4. Circle Ops — COMPLETE (correctness + robustness)

Audited and documented behavior at circular genome boundaries. `slidingWindow()` now reads per-chromosome circularity from `Seqinfo`, defaults missing circularity metadata to circular for bacterial genomes, and has known-value boundary tests.

**Issues:** #112, #122, #129
**Key functions:** slidingWindow(), plot_metagene(), plot_tss_profile(), annotateSites() with proximity method

### 5. Layered Assays and Result Layers — COMPLETE (robustness)

Layered assay provenance and named differential-methylation result layers are now part of the v0.2.0 codebase. Running `diffMethyl()` with explicit `result_name` values preserves multiple result layers; unnamed runs continue to update the compatibility/default layer.

**Issues:** #118, #167, #212

### 6. Technical Rename to commaKit — COMPLETE (publishability)

Renamed package namespace, docs, vignettes, pkgdown configuration, and public references from `comma`/`CoMMA` to `commaKit`. The exported `commaData` class and `comma_example_data` dataset keep their names for API continuity.

**Issues:** #168–#173
**PR:** #207

### 7. Usability (usability)

Make commaKit usable by someone other than Carl. Documentation, guidance, error messages.

**Open issue:** #68
**Completed issues:** #62, #64, #134, #161

What this milestone looks like when done:
- Performance expectations documented
- Claire can work through the getting-started vignette and the troubleshooting guide without Carl's help

Already completed:
- Troubleshooting guide for data import
- Method selection guidance for diffMethyl backends
- commaData object vignette
- Import troubleshooting guide
- Onboarding guide function names synchronized with the API

### 8. Real-World Readiness (robustness + publishability)

Make commaKit handle real data and be ready for broader distribution.

**Issues:** #67, #68, #76, #206, #232

This is the lowest priority milestone. Bioconductor submission is way down the list.

---

## Versioning Policy

- **Dev versions:** `x.y.z.9000` (Bioconductor convention for development)
- **Releases:** `x.y.z` (no `.9000` suffix)
- **Minor bumps** (0.2.0, 0.3.0): coherent feature sets, API changes, milestones
- **Patch bumps** (0.2.1, 0.2.2): accumulated small fixes, stable snapshots
- **Version bumps are deliberate** — Carl decides when to cut a release
- **0.99.0** is reserved for actual Bioconductor submission
- **Don't claim stability before it's earned** — no premature 1.0

---

## Post-v1.0 Feature Roadmap

These are aspirational future features. They are not current commitments.

### v1.1 — Effect Size Shrinkage

Goal: Bring DESeq2-style shrinkage thinking to methylation effect sizes.

- `lfcShrink()` equivalent for `delta_beta`
- Empirical Bayes prior on delta_beta
- `plot_effect_size()` visualization

Why: Stabilizes noisy site-level effect sizes in small-sample experiments.

### v1.2 — DMR Calling

Goal: Add region-level differential methylation, not just site-level testing.

- `callDMR()` function
- Interface to bsseq/DSS DMR methods or custom sliding-window method
- `plot_manhattan()` genome-wide DM landscape

Why: Biological interpretation often happens at regions/genes, not individual bases.

### v1.3 — Batch Effects & Complex Designs

Goal: Support more realistic experimental designs.

- `~ batch + condition` in all backends
- Multi-factor formula support
- Contrast specification in `results()`

Why: Real experiments have batches, strains, timepoints, and interactions.

### v1.4 — QC Report

Goal: Give users an automatic, citeable QC summary.

- `commaQC()` — runs all QC checks, stores results in metadata
- `qcReport()` — printable HTML/PDF summary

Why: Users need to know whether their data are usable before testing.

### v1.5 — VST & IHW

Goal: Better transformations and multiple-testing power.

- Variance-stabilizing transformation (VST/rlog equivalent)
- Independent hypothesis weighting (IHW package)

Why: Improves exploratory analysis and potentially increases detection power.

---

## References

- `dev/knowledge/test-quality.md` — what tests are strong, weak, or missing
- `dev/knowledge/known-issues.md` — bugs, gotchas, edge cases
- `dev/knowledge/design-decisions.md` — why the package is designed this way
- `dev/knowledge/git-discipline.md` — branching and versioning conventions
- `dev/knowledge/branching-releases.md` — release strategy
- [GitHub Issues](https://github.com/carl-stone/commaKit/issues) — tactical work items
