# comma Product Requirements Document (PRD)

**Package:** comma (Comparative Methylomics for Microbial Analysis)
**Author:** Carl Stone, Vanderbilt University
**Current version:** 0.1.0.9000
**Target version:** 0.1.0 internal checkpoint; 0.99.0 for eventual Bioconductor submission
**Date:** 2026-05-14

---

## 1. Product Vision

comma is the **DESeq2 of bacterial methylation analysis** — a single, opinionated R package that takes you from raw Oxford Nanopore methylation calls to publication-quality results, with a unified data container, multiple statistical backends, and rich visualization. It is modification-type agnostic (6mA, 5mC, 4mC), works with any monoploid bacterial genome, and follows the DESeq2-style workflow: construct → QC → annotate → test → visualize → export.

### Why comma exists

Existing methylation tools (methylKit, DSS, nanomethyR) were built for bisulfite sequencing of eukaryotic CpG methylation. They don't handle:
- **Multiple modification types** in one dataset (6mA + 5mC + 4mC from the same Nanopore run)
- **Modification context** as a first-class variable (6mA_GATC vs 6mA_ACCACC are biologically distinct)
- **Bacterial genome structure** (circular chromosomes, operons, small genomes, no CpG island paradigm)
- **Oxford Nanopore-native formats** (modkit pileup, Dorado BAM, MM/ML tags)
- **End-to-end workflow** from raw calls to enrichment analysis in one package

comma fills this gap.

---

## 2. Target Users

| User | Needs | Pain points |
|------|-------|-------------|
| Microbial epigenomic researchers | Full workflow from Nanopore calls to differential methylation | Cobbling together methylKit + custom scripts + ggplot2 |
| Bioinformaticians supporting microbiology labs | Reproducible, scriptable pipeline | No unified container; every lab reinvents the wheel |
| Bioconductor users | Familiar S4/SummarizedExperiment interface | Existing tools don't extend SummarizedExperiment |
| Carl's lab & collaborators | Ship papers | Manual analysis is slow and error-prone |

---

## 3. Core Workflow (v1.0)

```
commaData()          →  Construct unified data object
  ↓
methylomeSummary()   →  Per-sample QC
coverageDepth()      →  Coverage diagnostics
varianceByDepth()    →  Variance vs depth check
plot_coverage()      →  Visualize coverage
plot_methylation_distribution()
plot_pca()           →  Sample-level QC
  ↓
annotateSites()      →  Map sites to genomic features
plot_metagene()      →  Feature-level profiles
plot_tss_profile()   →  TSS-centered profiles
plot_genome_track()  →  Genome browser view
  ↓
diffMethyl()         →  Differential methylation testing
results()            →  Extract results table
filterResults()      →  Filter to significant sites
plot_volcano()       →  Volcano plot
plot_heatmap()       →  Heatmap of top DM sites
  ↓
enrichMethylation()  →  GO/KEGG enrichment
buildKEGGTermGene()  →  Offline KEGG path
buildKEGGGeneIDMap() →  Gene ID mapping
  ↓
writeBED()           →  Export to BED
```

---

## 4. Feature Inventory — Current State

### 4.1 Data Infrastructure

| Feature | Status | Notes |
|---------|--------|-------|
| `commaData` S4 class (extends SummarizedExperiment) | ✅ Done | Validity method, show() |
| modkit pileup parser | ✅ Done | Primary format |
| Dorado BAM parser (MM/ML tags) | ✅ Done | CIGAR decoding |
| Megalodon parser (legacy) | ✅ Done | Per-read aggregation |
| `commaData()` constructor | ✅ Done | Multi-sample, multi-mod-type |
| `expected_mod_contexts` filtering | ✅ Done | Drops unexpected mod_type × motif combos |
| `mod_context` rowData column | ✅ Done | Required in current data model |
| All accessors (methylation, coverage, sampleInfo, siteInfo, modTypes, modContexts, motifs, genome, annotation, motifSites) | ✅ Done | S4 methods |
| Subsetting (`[`, `subset()`) | ✅ Done | By mod_type, mod_context, condition, chrom |
| `comma_example_data` | ✅ Done | 588 sites, 6 samples, 2 mod types |

### 4.2 Analysis

| Feature | Status | Notes |
|---------|--------|-------|
| `annotateSites()` (list-column design) | ✅ Done | 4 keep modes: all/overlap/proximity/metagene |
| `loadAnnotation()` (GFF3/BED) | ✅ Done | feature_subtype preservation |
| `findMotifSites()` | ✅ Done | Both strands, palindromic |
| `slidingWindow()` | ✅ Done | Circular wrap, genome-size inference |
| `methylomeSummary()` | ✅ Done | Per-sample stats |
| `coverageDepth()` | ✅ Done | Windowed depth |
| `varianceByDepth()` | ✅ Done | Binned variance |
| `mValues()` | ✅ Done | M-value transformation |
| `diffMethyl()` (beta-binomial) | ✅ Done | Default backend |
| `diffMethyl()` (methylKit) | ✅ Done | Fisher's exact / logistic regression |
| `diffMethyl()` (limma) | ✅ Done | eBayes |
| `diffMethyl()` (quasi_f) | ✅ Done | Quasi-binomial F-test |
| `results()` / `filterResults()` | ✅ Done | DESeq2-style extraction |
| `enrichMethylation()` (GO/KEGG ORA + GSEA) | ✅ Done | gene_role = target/regulator/both |
| `buildKEGGTermGene()` | ✅ Done | Offline KEGG (2 API calls + RDS cache) |
| `buildKEGGGeneIDMap()` | ✅ Done | Symbol ↔ KEGG ID mapping |
| `writeBED()` | ✅ Done | 9-col BED with RGB |

### 4.3 Visualization

| Feature | Status | Notes |
|---------|--------|-------|
| `plot_coverage()` | ✅ Done | Per-sample or combined |
| `plot_methylation_distribution()` | ✅ Done | Beta value density |
| `plot_pca()` | ✅ Done | M-value PCA, return_data |
| `plot_genome_track()` | ✅ Done | Annotation track below |
| `plot_metagene()` | ✅ Done | Feature-normalized profile |
| `plot_tss_profile()` | ✅ Done | Regulatory element coloring |
| `plot_volcano()` | ✅ Done | DM volcano |
| `plot_heatmap()` | ✅ Done | ComplexHeatmap backend |

### 4.4 Documentation & Testing

| Feature | Status | Notes |
|---------|--------|-------|
| Roxygen2 docs for all exports | ✅ Done | All 32 exports documented |
| `?comma` package-level docs | ✅ Done | Bioconductor requirement |
| Vignette: getting-started | ✅ Done | End-to-end workflow |
| Vignette: multiple-modification-types | ✅ Done | Joint 6mA + 5mC |
| Test suite | ✅ Done | 938 tests, 0 failures |
| NEWS.md | ✅ Done | Full version history |
| README.Rmd | ✅ Done | Workflow showcase |
| pkgdown site | ✅ Done | GitHub Pages deployment |
| CI: R CMD check | ✅ Done | GitHub Actions |
| CI: pkgdown build | ✅ Done | GitHub Actions |
| CI: Rmd render | ✅ Done | GitHub Actions |

---

## 5. Gap Analysis — v1.0 Blockers

### 5.1 Bioconductor Submission Requirements

| Requirement | Status | Action needed |
|-------------|--------|---------------|
| `R CMD check --as-cran` zero errors/warnings | ⏳ Verify | Run and fix any issues |
| `BiocCheck::BiocCheck()` zero errors | ⏳ Run | Run and fix |
| Bundled data < 5 MB | ⏳ Verify | Check `data/` + `inst/extdata/` |
| Zenodo DOI | ⏳ Register | Before submission |
| Version policy documented | ✅ Done | Reset to 0.1.0.9000; 0.99.0 reserved for Bioconductor submission |
| `biocViews` complete | ⏳ Review | Current: Sequencing, Epigenetics, Coverage, DifferentialMethylation — may need more |
| Submit at contributions.bioconductor.org | ⏳ Pending | Final step |

### 5.2 Test Coverage Gaps

18 exported functions lack dedicated test files:

| Function | Issue |
|----------|-------|
| `buildKEGGGeneIDMap` | Tested within `test-buildKEGGTermGene.R` but no dedicated file |
| `coverageDepth` | Tested within `test-coverageAnalysis.R` (different name) |
| `enrichMethylation` | Tested within `test-enrichment.R` (different name) |
| `filterResults` | Tested within `test-results.R` (different name) |
| `findMotifSites` | Tested within `test-find_motif_sites.R` (different name) |
| `loadAnnotation` | Tested within `test-load_annotation.R` (different name) |
| `mValues` | Tested within `test-m_values.R` (different name) |
| `methylation` | Tested within `test-accessors.R` (accessor) |
| `modContexts` | Tested within `test-accessors.R` (accessor) |
| `modTypes` | Tested within `test-accessors.R` (accessor) |
| `motifSites` | Tested within `test-accessors.R` (accessor) |
| `motifs` | Tested within `test-accessors.R` (accessor) |
| `plot_methylation_distribution` | Tested within `test-plot_distribution.R` (different name) |
| `sampleInfo` | Tested within `test-accessors.R` (accessor) |
| `siteInfo` | Tested within `test-accessors.R` (accessor) |
| `subset` | Tested within `test-accessors.R` (accessor) |
| `varianceByDepth` | Tested within `test-coverageAnalysis.R` (different name) |

**Note:** Most of these are tested but under different file names (snake_case vs camelCase). The naming convention in `.claude/rules/testing.md` says `test-functionName.R`. The actual test files use snake_case for some functions. This is a naming inconsistency, not a coverage gap — but it should be cleaned up for Bioconductor review.

### 5.3 Documentation Gaps

| Issue | Action |
|-------|--------|
| README still says "300 sites" for example data | Update to 588 sites |
| README roadmap is outdated (stops at 0.5.0) | Update to reflect current state |
| Historical version numbers were informal and overstated maturity | Reset to 0.1.0.9000 fresh-start baseline |
| Vignettes don't cover enrichment workflow | Add enrichment section or vignette |
| No vignette for KEGG offline path | Add or extend existing vignette |

### 5.4 Code Quality

| Issue | Action |
|-------|--------|
| CI pinned to R 4.5 (S4Vectors C API break in R 4.6) | Monitor Bioconductor patch; update when ready |
| methylKit crashes on zero-variance/zero-coverage sites | Fixed (fc0fa10, 7a28e86) but add regression tests |
| `diag()` scalar trap risk | Audit for any usage |
| `S4Vectors::rename()` vs `dplyr::rename()` conflict | Audit for any usage |

---

## 6. v1.0 Release Criteria

The package ships when ALL of the following are true:

1. **R CMD check passes** with zero errors and zero warnings
2. **BiocCheck passes** with zero errors
3. **All 32 exported functions** have test files matching the naming convention
4. **All vignettes knit** without errors
5. **README is current** (reflects actual features, example data size, roadmap)
6. **AGENTS.md version matches DESCRIPTION**
7. **Bundled data < 5 MB**
8. **Zenodo DOI registered**
9. **Version policy documented** and DESCRIPTION/NEWS.md consistent

---

## 7. Out of Scope for v1.0

These are explicitly deferred. Do not implement without discussion:

- Multi-species comparative methylomics
- Integration with transcriptomics (RNA-seq correlation)
- Motif discovery (de novo)
- Phage/plasmid methylation analysis
- Shiny interactive browser
- Python or command-line interface
- Genome browser track export beyond BED (bigWig, etc.)
- DMR calling (region-level differential methylation)
- Effect size shrinkage (lfcShrink-style)
- Variance-stabilizing transformation (VST/rlog-style)
- Batch effect correction
- Multi-factor experimental designs beyond ~condition
- Per-read methylation analysis
- Haplotype-resolved methylation
- Cell-type deconvolution
- Integration with single-cell methylation data
