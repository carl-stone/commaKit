# comma Roadmap — Strategic Direction

**Last updated:** 2026-05-18
**Current version:** 0.1.0.9000 (dev)
**Next release:** 0.2.0 (schema v2 milestone)

This file is the strategic roadmap: where comma is going and why. Tactical work items are tracked in [GitHub Issues](https://github.com/carl-stone/CoMMA/issues).

---

## Strategic Goals

These are the reasons behind the work. They determine priority order.

1. **Correctness** — comma must produce results you can trust. If the stats are wrong, nothing else matters. This means real discriminating tests, not just smoke tests; known-quantity verification; and audits of known R gotchas.

2. **Usability** — Claire should be able to use comma independently. She needs docs, clear error messages, method selection guidance, and a package that doesn't surprise her with silent failures.

3. **Robustness** — comma should handle real data, not just the 588-site toy example. Edge cases, large genomes, weird callers, production-scale site counts.

4. **Publishability** — Bioconductor-ready when the time comes. This is a low priority. We'd rather have a great, stable package installable from GitHub than a rushed Bioconductor submission.

---

## Milestone Sequence

Milestones are coherent groups of work that advance one or more strategic goals. They are ordered by priority. Each milestone has a GitHub Milestone for issue tracking.

### 1. Schema v2 — *in progress* (correctness + robustness)

Restructure the `commaData` class for a stronger foundation. The current SE-based class with custom slots works but creates friction for future improvements.

**Issues:** #92 (done), #93, #94, #95, #96, #97, #98, #99
**GitHub Milestone:** commaData Schema v2
**Target version:** 0.2.0

Key changes:
- SE → RangedSummarizedExperiment (done, PR #100)
- genomeInfo → Seqinfo (#93)
- annotation/motifSites → metadata() (#94)
- mod_context derived on demand, not stored (#95)
- mod_type as factor, caller/min_coverage stored, site key convention (#96–#98)

### 2. Test Quality — *next* (correctness)

The test suite has too many smoke tests and not enough discriminating tests. We don't really believe the tests yet — they verify "doesn't crash" but not "produces correct results."

**Issues:** #73, #74, #75
**GitHub Milestone:** Test Quality (to be created)

What this milestone looks like when done:
- Plot tests verify data mappings (not just ggplot class)
- slidingWindow circular correctness verified with known boundary values
- enrichMethylation tested with real clusterProfiler (not fake TERM2GENE)
- Integration test across full pipeline is reliable and meaningful

See `dev/knowledge/test-quality.md` for the full audit.

### 3. Code Quality Audits (correctness)

Known R gotchas that could cause silent bugs. These are small, scoped, and high-confidence.

**Issues:** #65, #66
**GitHub Milestone:** Code Quality Audits (to be created)

What this milestone looks like when done:
- `diag()` scalar trap audited and either safe or guarded
- `rename()` masking audited and either safe or guarded
- Both findings documented in `dev/knowledge/known-issues.md`

### 4. Usability (usability)

Make comma usable by someone other than Carl. Documentation, guidance, error messages.

**Issues:** #62, #64, #68
**GitHub Milestone:** Usability (to be created)

What this milestone looks like when done:
- Troubleshooting guide for data import
- Method selection guidance for diffMethyl backends
- Performance expectations documented
- Claire can work through the getting-started vignette and the troubleshooting guide without Carl's help

### 5. Real-World Readiness (robustness + publishability)

Make comma handle real data and be ready for broader distribution.

**Issues:** #67, #70, #76, #77
**GitHub Milestone:** Real-World Readiness (to be created)

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

- `dev/PRD.md` — v1.0 product requirements and scope
- `dev/VISION.md` — long-term dream package
- `dev/README.md` — how the dev directory is organized
- [GitHub Issues](https://github.com/carl-stone/CoMMA/issues) — tactical work items
- `dev/knowledge/test-quality.md` — what tests are strong, weak, or missing
- `dev/knowledge/known-issues.md` — bugs, gotchas, edge cases
- `dev/knowledge/design-decisions.md` — why the package is designed this way
