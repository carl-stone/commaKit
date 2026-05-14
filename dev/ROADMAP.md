# comma Roadmap — Strategic Direction

**Last updated:** 2026-05-14  
**Current version:** 0.1.0.9000 (fresh-start dev baseline)  
**Submission target:** 0.99.0 when ready for Bioconductor

This file is the strategic roadmap: where comma is going. Tactical work items live in `dev/BACKLOG.md`. Current work status lives in `dev/STATUS.md`.

---

## Product Positioning

comma aims to be the **DESeq2 of bacterial methylation analysis**: an opinionated Bioconductor-style package that takes users from Oxford Nanopore methylation calls to QC, annotation, differential methylation, enrichment, visualization, and export.

Core principles:
- Modification-type agnostic: 6mA, 5mC, 4mC, and future modifications
- Bacterial-first: circular chromosomes, compact genomes, overlapping features
- Bioconductor-native: S4, `SummarizedExperiment`, `GenomicRanges`
- End-to-end workflow: one coherent object and pipeline

---

## v1.0 Status

**Engineering status:** Feature-complete for current v1.0 scope.  
**Submission status:** Not submitting yet. Package reset to `0.1.0.9000` on 2026-05-14 to mark the first disciplined package-development phase. `0.99.0` is reserved for eventual Bioconductor submission.

| Area | Status | Notes |
|------|--------|-------|
| Core feature set | ✅ Complete | 32/32 exports implemented |
| Tests | 🟡 Mixed | 938 passing; core tests strong; plot/integration gaps remain |
| Documentation | ✅ In sync | Function docs match code; examples sampled and pass |
| R CMD check | ✅ Clean | 0 errors, 0 warnings, 0 notes |
| BiocCheck | 🟡 Mostly ready | Administrative tasks remain |
| User guidance | 🟡 Needs improvement | Troubleshooting and method-selection guidance missing |

For detailed quality findings, see:
- `dev/knowledge/test-quality.md`
- `dev/knowledge/known-issues.md`
- `dev/BACKLOG.md`

---

## Bioconductor Submission Path

When Carl decides comma is ready for Bioconductor submission:

| Task | Status | Owner | Notes |
|------|--------|-------|-------|
| R CMD check clean | ✅ Done | commaBot | 0 errors/warnings/notes |
| BiocCheck clean enough | ✅ Mostly done | commaBot | One admin requirement remains |
| Bundled data < 5 MB | ✅ Done | commaBot | 40 KB total (`data/` + `inst/extdata/`) |
| biocViews declared | ✅ Done | commaBot | Appropriate current set |
| Two vignettes | ✅ Done | commaBot | Getting started + multiple modification types |
| NEWS.md | ✅ Done | commaBot | Needs 0.99.0 entry at submission time |
| Carl registers on Bioc support site | ⏳ Pending | Carl | Required by Bioconductor |
| Zenodo DOI | ⏳ Pending | Carl/commaBot | Before submission |
| Version bump to 0.99.0 | ⏳ Pending | commaBot | Only when actually submitting to Bioconductor |
| Submit at contributions.bioconductor.org | ⏳ Pending | Carl/commaBot | Final step |

Tactical submission tasks are tracked in `dev/BACKLOG.md` as B-series items.

---

## Post-v1.0 Feature Roadmap

These are aspirational future versions. They are not current commitments.

### v1.1 — Effect Size Shrinkage

Goal: Bring DESeq2-style shrinkage thinking to methylation effect sizes.

Potential features:
- `lfcShrink()` equivalent for `delta_beta`
- Empirical Bayes prior on delta_beta
- `plot_effect_size()` visualization

Why it matters:
- Stabilizes noisy site-level effect sizes in small-sample experiments
- Helps users avoid over-interpreting extreme estimates from low-depth sites

---

### v1.2 — DMR Calling

Goal: Add region-level differential methylation, not just site-level testing.

Potential features:
- `callDMR()` function
- Interface to bsseq/DSS DMR methods or custom sliding-window method
- `plot_manhattan()` genome-wide DM landscape

Why it matters:
- Biological interpretation often happens at regions/genes, not individual bases
- Bacterial methylation may cluster around motifs, promoters, or operons

---

### v1.3 — Batch Effects & Complex Designs

Goal: Support more realistic experimental designs.

Potential features:
- `~ batch + condition` in all backends
- Multi-factor formula support
- Contrast specification in `results()`

Why it matters:
- Real experiments have batches, strains, timepoints, and interactions
- Current simple condition-vs-condition model is not enough long-term

---

### v1.4 — QC Report

Goal: Give users an automatic, citeable QC summary.

Potential features:
- `commaQC()` — runs all QC checks, stores results in metadata
- `qcReport()` — printable HTML/PDF summary

Why it matters:
- Users need to know whether their data are usable before testing
- Claire and future lab members need a standard QC workflow

---

### v1.5 — VST & IHW

Goal: Better transformations and multiple-testing power.

Potential features:
- Variance-stabilizing transformation (VST/rlog equivalent)
- Independent hypothesis weighting (IHW package)

Why it matters:
- Improves exploratory analysis and potentially increases detection power

---

## References

- `dev/PRD.md` — v1.0 product requirements and scope
- `dev/VISION.md` — long-term dream package
- `dev/BACKLOG.md` — tactical prioritized work items
- `dev/STATUS.md` — current work status
- `dev/knowledge/design-decisions.md` — why the package is designed this way
