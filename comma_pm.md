# comma — Project Management Document

**Comprehensive Methylation Analysis for Bacterial Nanopore Data**

> This document is the authoritative reference for the design, scope, and development roadmap of the `comma` R package. It is intended to be read by human developers and AI coding agents alike. All architectural decisions, rationale, and implementation guidance are recorded here.

-----

## 1. Background and Current State

### 1.1 History

`comma` was originally developed as `CoMMA` (Comparison of Microbial Methylated Adenines) by Carl Stone (Vanderbilt University) to support the analysis presented in a 2022 preprint on adenine methylation in *E. coli* K-12 MG1655. The package was built to reproduce a specific set of figures and analyses from that manuscript — it was never designed as a general-purpose tool.

The name `comma` replaces `CoMMA` for two reasons: (1) it is far easier to type, and (2) the original acronym implied exclusive focus on adenine methylation, which no longer reflects the package’s scope (see Section 1.3). The new name stands for **Comparative Methylomics for Microbial Analysis**, preserving the spirit of comparative analysis while being modification-type agnostic.

### 1.2 What the Current Codebase Contains

The repository (`carl-stone/CoMMA` on GitHub, to be renamed `carl-stone/comma`) as of early 2023 contains the following:

**Formally exported functions (NAMESPACE):**

- `annotateMethylSites()` — joins methylation site positions to a genomic feature table by overlap
- `annotateTSS()` — assigns relative position to transcription start sites within a user-defined window
- `methylRollingMedian()` — genome-wide sliding window median methylation, with `"exact"` (zoo-based) and `"fast"` (loop-based) modes

**Implemented but not exported (`R/` directory):**

- `annotateTTS()` — transcription termination site analog of `annotateTSS`; supports Rho-independent, Rho-dependent, and TU-inferred terminators
- `methylRollingMean()` — mean-based sliding window; faster matrix implementation
- `calculateMethylSiteDepth()` — windowed average sequencing depth
- `varByCoverage()` — variance in methylation change as a function of coverage depth
- `writeBED()` — exports methylation data as BED file; **hardcoded to absolute paths on the developer’s machine; non-functional for any other user**

**Root-level scripts (not part of the package):**

- `functions.R` — earlier, undocumented versions of the above functions; developmental scratchpad
- `methylKitGATC.R` — 513-line analysis script from the original manuscript; contains the most sophisticated logic in the codebase (multi-sample differential methylation via methylKit, Euler diagrams, sigma factor binding site analysis) but is not packaged, uses hardcoded Box Drive paths, and is not usable by others

**Tests:**

- One test file exists (`tests/testthat/test-annotateMethylSites.R`) containing only the default placeholder: `expect_equal(2 * 2, 4)`. No actual package functionality is tested.

**Bundled data:**

- Nine `.rda` datasets (~9.5 MB) in `data/`, all specific to *E. coli* K-12 MG1655. These are used in README examples.

**CI/CD:**

- GitHub Actions workflows for `rcmdcheck` (R 3.6.3 and 4.1.1, macOS) and auto-rendering of `.Rmd` files

**Package metadata:**

- Version: 0.1.0
- Hard imports: `methylKit`, `forcats`, `ggplot2`, `tidyr`, `dplyr`
- License: MIT

### 1.3 The Critical Scope Change: Beyond Adenine

When `CoMMA` was originally conceived, it was focused exclusively on N6-methyladenine (6mA), the modification made by Dam methyltransferase at GATC motifs in *E. coli*. Since then it has become clear that other methylation types — most importantly N4-methylcytosine (4mC) and 5-methylcytosine (5mC) — have significant and understudied biological roles even in organisms like *E. coli*. Nanopore sequencing with tools like Dorado can now detect all three modification types simultaneously from a single sequencing run.

**`comma` must be modification-type agnostic from the ground up.** Every data structure, function signature, and analysis module must treat methylation type as a first-class parameter, not an assumption. The package should work equally well for:

- 6mA (Dam, e.g., GATC in *E. coli*)
- 4mC (various bacterial MTases)
- 5mC (dcm and analogs, e.g., CCWGG in *E. coli*)
- Any other modification type that nanopore callers produce

This is not merely a software generalization — it is a scientific one. The package’s value proposition is comparative methylomics across all modification types, not just adenine.

### 1.4 Known Weaknesses to Carry Forward

These are specific issues in the current code that must be addressed and not replicated in new development:

1. **No data object.** Functions pass raw dataframes with no enforced schema. Column names are assumed but not validated.
1. **Organism specificity.** Genome size is hardcoded to 4,641,652 (MG1655). Annotations are bundled as MG1655-specific `.rda` files. There is no mechanism to use any other organism.
1. **Input format lock-in.** Only pre-formatted Megalodon output tables are supported. Megalodon is deprecated; current ONT standard is Dorado + modkit.
1. **Performance.** `annotateMethylSites()` and `annotateTSS()` use nested for-loops over every genomic position — this is O(n×m) and unacceptably slow on any reasonably sized genome.
1. **Circular genome handling.** The sliding window wrap-around logic is not generalized and assumes MG1655 size.
1. **No differential methylation function.** The most important analysis capability exists only in an unpackaged script with hardcoded paths.
1. **No tests.** Zero functional tests exist.
1. **`writeBED()` is broken for any user other than the original developer.**

-----

## 2. Design Philosophy

### 2.1 The DESeq2 Model

The guiding analogy for `comma`‘s architecture is DESeq2 for RNA-seq. DESeq2’s power comes not only from its statistical methods but from its data infrastructure:

- A single, well-defined input object (`DESeqDataSet`) that carries all data and metadata
- A clear, linear pipeline: construct → normalize → fit → test → results
- Clean separation between data, transformation, and inference
- Every function knows what it receives; the user always knows what they’re working with
- Output objects that are themselves usable by downstream functions

`comma` must follow this model. The central object (`commaData`, see Section 3.1) is the most important design decision in the entire package. Everything else — parsers, analysis functions, visualization — is built around it.

### 2.2 Principles

**Modification-type agnostic.** Every function that touches methylation data must accept a `mod_type` argument or infer it from the `commaData` object. No function should assume 6mA, GATC, or any other specific modification or motif.

**Genome-agnostic.** No organism-specific data should be bundled in the package except small synthetic datasets for testing. Genome annotations are loaded from standard formats (GFF3, BED) at runtime. Genome size is always a user-supplied parameter or inferred from input files.

**Standard input formats.** The package should speak the formats that nanopore methylation tools actually produce. The primary targets are Dorado (`.bam` with MM/ML tags) and modkit (`.bed` output), with Megalodon supported for backward compatibility. Input parsing should be isolated in a dedicated module so new caller formats can be added without touching analysis code.

**Tidy output.** All analysis functions should return tidy dataframes (or `commaData` objects that contain them) suitable for direct use with ggplot2 and tidyverse tools. Users should never need to reshape output before plotting.

**ggplot2-native visualization.** All built-in plot functions return `ggplot` objects, not rendered images. This means users can add their own layers, themes, and annotations on top of any `comma` plot.

**Performance.** Any function that iterates over genomic positions must be implemented with vectorized operations or Rcpp, not R for-loops. The current loop-based annotation functions are unacceptable for production use.

**Bioconductor-ready.** Given that a core dependency (`methylKit`) is on Bioconductor, and given the target audience of bioinformaticians, `comma` should be designed for Bioconductor submission from the start. This means following Bioconductor coding standards, using `SummarizedExperiment` as a base class or parallel, writing proper vignettes, and passing BiocCheck.

-----

## 3. Architecture

### 3.1 The `commaData` Object (Do This First)

The `commaData` S4 class is the foundation of the entire package. It is analogous to `DESeqDataSet` or `SummarizedExperiment`. All analysis functions accept a `commaData` object and most return one (possibly with added slots).

**Required slots:**

```
commaData
├── methylation        # SummarizedExperiment-like: sites × samples matrix of beta values
├── coverage           # sites × samples matrix of read depth at each site
├── rowData            # per-site metadata: chromosome, position, strand, motif, mod_type
├── colData            # per-sample metadata: sample_name, condition, replicate, caller, file_path
├── genomeInfo         # chromosome names and sizes (equivalent of Seqinfo)
├── annotation         # GRanges object of genomic features, loaded from GFF3/BED
├── motifSites         # GRanges of all motif instances in the genome (e.g., all GATC positions)
└── metadata           # list: package version, creation date, any user-defined fields
```

**Key design decisions:**

- `rowData` must include `mod_type` (e.g., `"6mA"`, `"5mC"`, `"4mC"`) as a first-class column. A single `commaData` object can contain multiple modification types simultaneously, filtered by `mod_type` as needed.
- `methylation` stores beta values (proportion methylated reads, 0–1) not raw counts. Raw counts (modified reads, total reads) are stored separately if needed.
- Sites with coverage below a minimum threshold (set at object creation, default: 5×) are stored as `NA` not removed, so the site universe is consistent across samples.
- The object should have a `show()` method that prints a clean summary (n sites, n samples, modification types present, organism, etc.)

**Constructor function:**

```r
commaData(
  files,           # named character vector or data frame: sample_name → file_path
  colData,         # data frame with at minimum: sample_name, condition, replicate
  genome,          # BSgenome object, FASTA path, or named integer vector of chr sizes
  annotation,      # GFF3 file path or GRanges object (optional but recommended)
  mod_type,        # character: "6mA", "5mC", "4mC", or NULL to auto-detect from files
  motif,           # character: regex motif (e.g., "GATC") or NULL to use all sites
  min_coverage,    # integer: minimum read depth to include a site (default: 5)
  caller           # character: "dorado", "modkit", "megalodon" (auto-detected if possible)
)
```

### 3.2 Module Structure

The package is organized into five functional modules. These correspond to the major analysis phases a user would move through.

```
R/
├── data_input/
│   ├── parse_dorado.R          # Parse Dorado BAM with MM/ML tags
│   ├── parse_modkit.R          # Parse modkit BED output
│   ├── parse_megalodon.R       # Parse Megalodon output (backward compat)
│   └── commaData_constructor.R # commaData S4 class definition and constructor
│
├── genome/
│   ├── load_annotation.R       # GFF3/BED → GRanges; standardize feature types
│   ├── find_motif_sites.R      # FASTA + motif regex → GRanges of all motif positions
│   └── genome_utils.R          # Circular genome arithmetic, seqinfo helpers
│
├── characterization/
│   ├── methylome_summary.R     # Distribution stats, saturation, per-sample QC
│   ├── sliding_window.R        # methylRollingMedian/Mean, generalized
│   ├── feature_annotation.R    # annotateSites() — replaces annotateMethylSites/TSS/TTS
│   └── coverage_analysis.R     # calculateMethylSiteDepth, varByCoverage
│
├── differential/
│   ├── diff_methylation.R      # Main differential methylation interface
│   ├── beta_binomial.R         # Beta-binomial model implementation
│   ├── methylkit_wrapper.R     # methylKit-based testing (alternative method)
│   └── multiple_testing.R      # Genome-wide multiple comparison correction
│
└── visualization/
    ├── plot_distribution.R     # Beta value density, ECDF plots
    ├── plot_genome.R           # Genome browser-style methylation track
    ├── plot_volcano.R          # Differential methylation volcano plot
    ├── plot_heatmap.R          # Heatmap of diff methylated sites across samples
    ├── plot_metagene.R         # TSS/TTS/feature proximity metagene plots
    └── plot_pca.R              # PCA on methylation profiles
```

-----

## 4. Implementation Roadmap

Development is organized into five phases. Each phase must be substantially complete before beginning the next, because later phases depend on the data infrastructure built in earlier ones.

-----

### Phase 1 — Data Infrastructure (Foundation)

**Goal:** Define the `commaData` object and build input parsers. This is the most important phase. Nothing else is meaningful without a stable, well-defined data object.

**Tasks:**

1. **Define the `commaData` S4 class** with all slots described in Section 3.1. Write `validity()` method to enforce schema. Write `show()` method.
1. **Write the modkit parser** (`parse_modkit.R`). modkit BED output is the current ONT standard and should be the primary supported format. The modkit `pileup` command produces a BED file with columns: chrom, start, end, mod_code, score, strand, coverage, mod_frequency, n_mod, n_canonical, n_other_mod, n_delete, n_fail, n_diff, n_no_call. Parse this into the internal format, handling multiple `mod_code` values (a = 6mA, m = 5mC, 21839 = 4mC).
1. **Write the Dorado BAM parser** (`parse_dorado.R`). Reads `.bam` files with MM/ML tags directly using `Rsamtools`. This is lower priority than modkit (since modkit is the recommended intermediate step) but important for users who want to start from BAM.
1. **Write the Megalodon parser** (`parse_megalodon.R`). Adapts the existing input handling from the current codebase. Lower priority but needed for backward compatibility with existing data.
1. **Write `find_motif_sites()`** — takes a genome (as `BSgenome` or FASTA path) and a motif regex, returns a `GRanges` of all motif instances on both strands. This is the bacterial equivalent of CpG island finding. Example: `find_motif_sites(genome = "MG1655.fa", motif = "GATC")`.
1. **Write `load_annotation()`** — reads a GFF3 or BED file, returns a `GRanges` with standardized feature type vocabulary. Should handle both NCBI and Ensembl GFF3 conventions.
1. **Write the `commaData()` constructor** that calls the appropriate parser based on `caller` argument, runs `find_motif_sites()` if `motif` is specified, loads annotation if provided, applies `min_coverage` filtering (setting low-coverage sites to NA), and returns a valid `commaData` object.
1. **Write accessor functions:** `methylation()`, `coverage()`, `sampleInfo()`, `siteInfo()`, `modTypes()`, `genome()`, `annotation()`. These should follow Bioconductor accessor conventions.
1. **Write subsetting methods:** `[` for sites, `[[` for samples, `subset()` by condition/mod_type/chromosome.

**What “done” looks like:** A user can go from a modkit BED file and a GFF3 to a `commaData` object in three lines of R. The object prints cleanly and can be subsetted.

-----

### Phase 2 — Genome Generalization

**Goal:** Remove all organism-specific assumptions from the codebase. After this phase, the package should work on any sequenced bacterial genome with a GFF3.

**Tasks:**

1. **Audit every function** for hardcoded genome size (4641652), hardcoded chromosome name (“U00096.3”, “NC_000913.3”), hardcoded feature type names (“Transcription-Units”, “DNA-Binding-Sites”), and MG1655-specific logic. Replace all of these with parameters drawn from the `commaData` object.
1. **Remove organism-specific bundled data.** Delete or move the nine MG1655 `.rda` files in `data/`. Replace with a small synthetic dataset (`comma_example_data`) — a simulated 100kb genome with 3 samples, 2 conditions, 6mA and 5mC sites — used in vignettes and tests. This should be <500KB.
1. **Generalize circular genome arithmetic.** The sliding window wrap-around in `methylRollingMedian()` assumes MG1655 size. Refactor `genome_utils.R` to provide general circular indexing given any genome size from `commaData@genomeInfo`.
1. **Generalize the motif system.** Currently GATC is implicitly assumed everywhere. The `commaData` object carries `motifSites` for the user’s motif of interest. All annotation and windowing functions should use these rather than any hardcoded motif.

**What “done” looks like:** The package installs and runs its core functions on a GFF3 + modkit output from *Caulobacter crescentus*, *Bacillus subtilis*, or any other organism without modification.

-----

### Phase 3 — Refactor Existing Functions

**Goal:** Take the existing analysis functions (annotation, sliding window) and rebuild them properly: accepting `commaData`, using vectorized operations, returning tidy output.

**Tasks:**

1. **Rewrite `annotateSites()`** (replaces `annotateMethylSites`, `annotateTSS`, `annotateTTS`). Uses `GenomicRanges::findOverlaps()` instead of nested for-loops — this is the critical performance fix. A vectorized overlap join that takes 30 seconds in a loop should take milliseconds with `findOverlaps`. Accepts a `commaData` object and returns it with an updated `rowData` containing feature annotations. The function should support:
- Overlap annotation (is this site within a gene, promoter, etc.)
- Proximity annotation (distance to nearest TSS, TTS, or any feature)
- Metagene position (relative position within feature as fraction 0–1)
1. **Rewrite `methylRollingMedian()` and `methylRollingMean()`** to accept `commaData`, use genome size from `commaData@genomeInfo`, operate per-sample (returning all samples), and use `zoo::rollapply` consistently (the fast matrix loop implementation can be an internal helper).
1. **Rewrite `calculateMethylSiteDepth()`** — generalize window size and return a tidy dataframe.
1. **Rewrite `varByCoverage()`** — currently uses non-standard column names; generalize to accept column name arguments.
1. **Generalize `writeBED()`** — remove hardcoded paths entirely. Accept a `commaData` object, an output path, and an optional sample name. This should be a proper exported function.
1. **Move `functions.R` and `methylKitGATC.R`** from the repo root. `functions.R` can be deleted (superseded). `methylKitGATC.R` should move to `inst/scripts/` as a historical reference, with a comment noting that its logic has been incorporated into the package properly.

**What “done” looks like:** The three currently-exported functions have been replaced with better implementations. The root-level scripts are gone. All `R/` functions accept `commaData` and return tidy output.

-----

### Phase 4 — Differential Methylation

**Goal:** Package the differential methylation analysis that currently lives in `methylKitGATC.R` into a proper, exported, documented module. This is the scientifically most important phase.

**Tasks:**

1. **Design the `diffMethyl()` interface.** This is the main user-facing function, analogous to `DESeq2::DESeq()`:

```r
diffMethyl(
  object,          # commaData object
  design,          # formula: ~ condition, or ~ condition + batch
  mod_type,        # which modification type to test: "6mA", "5mC", etc.
  method,          # "beta_binomial" (default) or "methylkit"
  min_coverage,    # minimum per-site coverage to include in testing (default: 10)
  min_samples,     # minimum samples with coverage at a site to include (default: all)
  ...              # passed to underlying method
)
```

Returns an updated `commaData` object with differential methylation results stored in a `results` slot.

1. **Implement the beta-binomial model** (`beta_binomial.R`). methylKit uses logistic regression, which treats the outcome as binary (each read either is or isn’t methylated). A beta-binomial model is more appropriate because it accounts for overdispersion in methylation rates across reads at the same site — the exact same problem that negative binomial models solve for RNA-seq counts. The model: for each site, the number of modified reads follows Binomial(n, p) where p is Beta-distributed across samples within a group. Test H0: p_group1 = p_group2 using a likelihood ratio test or Wald test. Effect size should be reported as difference in mean beta values (Δβ) — biologically interpretable and on a [−1, 1] scale.
1. **Wrap methylKit** (`methylkit_wrapper.R`) as an alternative method. Some users will prefer it for familiarity. The wrapper should accept a `commaData` object, extract the data methylKit needs, run the methylKit pipeline, and return results in the same tidy format as the beta-binomial method.
1. **Implement genome-wide multiple testing correction** (`multiple_testing.R`). Benjamini-Hochberg across all tested sites is the default. Also offer q-value (Storey) for users who want FDR control. Be explicit in documentation about what the denominator is (all sites tested, not just significant ones).
1. **Write a `results()` method** that extracts the differential methylation table from a `commaData` object as a tidy dataframe with columns: chrom, position, strand, mod_type, mean_beta_condition1, mean_beta_condition2, delta_beta, stat, pvalue, padj, any feature annotations present in rowData.
1. **Write `filterResults()`** — convenience function to filter results by |Δβ| and padj thresholds.

**What “done” looks like:** A user with a `commaData` object can run `dm <- diffMethyl(object, ~ condition)` and get a results table with per-site statistics, effect sizes, and adjusted p-values. The analysis handles multiple modification types simultaneously if present.

-----

### Phase 5 — Visualization, Tests, and Release Preparation

**Goal:** Complete the visualization layer, write real tests, write vignettes, and prepare for Bioconductor submission.

**Tasks:**

**Visualization:**

All plot functions follow the naming convention `plot_*()` and return `ggplot` objects.

- `plot_methylation_distribution(object, mod_type, per_sample)` — beta value density or ECDF per sample; the current README example, generalized
- `plot_genome_track(object, chromosome, start, end, mod_type, annotation)` — genome browser-style plot: position on x-axis, beta on y-axis, feature annotation as colored rectangles below
- `plot_metagene(object, feature, mod_type, window)` — average methylation relative to a feature midpoint (TSS, TTS, or any feature in annotation); replaces the `annotateTSS` + manual plotting workflow
- `plot_volcano(results, delta_beta_threshold, padj_threshold)` — standard volcano plot for differential methylation results
- `plot_heatmap(results, object, n_sites, annotation_cols)` — heatmap of top differentially methylated sites across samples, with optional sample/site annotation bars
- `plot_pca(object, mod_type, color_by, shape_by)` — PCA on per-sample methylation profiles; essential for QC and exploratory analysis
- `plot_coverage(object, per_sample)` — coverage distribution across sites; QC plot

**Tests:**

Every exported function needs tests. Use `testthat` with a small synthetic `commaData` object (built from the `comma_example_data` bundled dataset) as the test fixture. Minimum test coverage targets:

- Constructor: valid input succeeds, invalid input throws informative error, all parsers produce identical structure
- `annotateSites()`: correct feature assignment, correct relative position calculation, handles sites with no overlapping features
- `methylRollingMedian()`: output length equals genome size, handles circular wrap-around, NAs handled correctly
- `diffMethyl()`: correct identification of simulated differentially methylated sites, multiple testing correction applied, handles missing data
- All `plot_*()` functions: returns a `ggplot` object without error

**Vignettes:**

Write at minimum two vignettes:

1. **“Getting Started with comma”** — end-to-end workflow using `comma_example_data`: create `commaData` → characterize methylome → identify differentially methylated sites → visualize results. This should be completable in under 5 minutes on a laptop.
1. **“Working with Multiple Modification Types”** — demonstrates the 6mA + 5mC joint analysis workflow; shows how to subset `commaData` by `mod_type`, compare modification patterns, and visualize both types simultaneously.

**Documentation:**

- Every exported function: full roxygen2 documentation with `@param` for every argument, `@return` describing the output, `@examples` that run without error
- `annotateTSS()` in the current codebase has stub documentation (“A dataframe.”, “A string.”) — this pattern must not appear in any new code
- Package-level documentation (`?comma`) explaining the overall workflow

**Bioconductor submission preparation:**

- Run `BiocCheck::BiocCheck()` and address all errors and warnings
- Ensure package passes `R CMD check --as-cran` with no errors or warnings
- Reduce installed package size: bundled data should be <5MB total
- Add `NEWS.md` with version history
- Register the package DOI via Zenodo before submission

-----

## 5. What to Cut or Defer

**Cut entirely:**

- `writeBED()` in its current form (hardcoded paths). Replace with a generalized version in Phase 3 or defer to a later release.
- The bundled MG1655 `.rda` datasets (all nine). Replace with `comma_example_data`.
- `functions.R` from the repo root.

**Defer to post-1.0:**

- Dorado BAM parsing (modkit covers the immediate need)
- Multi-organism joint analysis (comparing methylomes across species)
- Integration with genome browsers via track export formats beyond BED (bigWig, etc.)
- Shiny-based interactive genome browser

-----

## 6. File and Directory Structure (Target State)

```
comma/
├── DESCRIPTION
├── NAMESPACE
├── NEWS.md
├── README.md
├── README.Rmd
├── LICENSE.md
├── .github/
│   └── workflows/
│       ├── r-cmd-check.yml       # rcmdcheck on push/PR, multiple R versions
│       └── render-rmarkdown.yaml
├── R/
│   ├── commaData_class.R         # S4 class definition, show(), validity()
│   ├── commaData_constructor.R   # commaData() constructor
│   ├── accessors.R               # methylation(), coverage(), sampleInfo(), etc.
│   ├── parse_modkit.R
│   ├── parse_dorado.R
│   ├── parse_megalodon.R
│   ├── load_annotation.R
│   ├── find_motif_sites.R
│   ├── genome_utils.R
│   ├── annotateSites.R
│   ├── sliding_window.R
│   ├── methylome_summary.R
│   ├── coverage_analysis.R
│   ├── diffMethyl.R
│   ├── beta_binomial.R
│   ├── methylkit_wrapper.R
│   ├── multiple_testing.R
│   ├── results_methods.R
│   ├── plot_distribution.R
│   ├── plot_genome_track.R
│   ├── plot_metagene.R
│   ├── plot_volcano.R
│   ├── plot_heatmap.R
│   ├── plot_pca.R
│   └── plot_coverage.R
├── data/
│   └── comma_example_data.rda    # Small synthetic dataset for vignettes and tests
├── inst/
│   ├── extdata/
│   │   ├── example_modkit.bed    # Small modkit output used in vignette
│   │   ├── example.gff3          # GFF3 annotation for example genome
│   │   └── example_genome.fa     # Small synthetic FASTA (~100kb)
│   └── scripts/
│       └── methylKitGATC_historical.R   # Original analysis script, for reference only
├── man/
│   └── [roxygen-generated .Rd files]
├── tests/
│   └── testthat/
│       ├── test-commaData.R
│       ├── test-parsers.R
│       ├── test-annotateSites.R
│       ├── test-sliding_window.R
│       ├── test-diffMethyl.R
│       └── test-plots.R
└── vignettes/
    ├── getting-started.Rmd
    └── multiple-modification-types.Rmd
```

-----

## 7. Dependencies

**Hard imports (Imports in DESCRIPTION):**

|Package               |Purpose                                                       |
|----------------------|--------------------------------------------------------------|
|`GenomicRanges`       |Core genomic interval arithmetic; replaces for-loop annotation|
|`SummarizedExperiment`|Base class infrastructure for `commaData`                     |
|`IRanges`             |Range operations (via GenomicRanges)                          |
|`Rsamtools`           |BAM file parsing for Dorado input                             |
|`methylKit`           |Differential methylation (alternative method)                 |
|`zoo`                 |Rolling window operations in `methylRollingMedian`            |
|`ggplot2`             |All visualization                                             |
|`dplyr`               |Data manipulation                                             |
|`tidyr`               |Data reshaping                                                |
|`BiocGenerics`        |Bioconductor generic methods                                  |

**Soft dependencies (Suggests in DESCRIPTION):**

|Package              |Purpose                                        |
|---------------------|-----------------------------------------------|
|`BSgenome`           |Genome sequence access for `find_motif_sites`  |
|`rtracklayer`        |GFF3 import via `import()`                     |
|`ComplexHeatmap`     |Heatmap visualization (preferred over pheatmap)|
|`ggrepel`            |Labeled points in volcano plots                |
|`testthat`           |Testing                                        |
|`knitr` / `rmarkdown`|Vignettes                                      |

**Dependencies to avoid:** Do not add `tidyverse` as a dependency — import individual packages (`dplyr`, `tidyr`, etc.) instead. This is a Bioconductor requirement.

-----

## 8. Naming Conventions

- **Class:** `commaData` (lowercase c, camelCase)
- **Constructor:** `commaData()` (same as class name, Bioconductor convention)
- **Analysis functions:** `verbNoun()` camelCase — `annotateSites()`, `diffMethyl()`, `findMotifSites()`
- **Plot functions:** `plot_noun()` snake_case with `plot_` prefix — `plot_volcano()`, `plot_metagene()`
- **Internal functions:** prefix with `.` — `.parseBetaValues()`, `.circularIndex()`
- **Arguments:** snake_case throughout — `mod_type`, `min_coverage`, `position_col`

-----

## 9. Version Plan

|Version|Content                                                                           |
|-------|----------------------------------------------------------------------------------|
|0.2.0  |Phase 1 + 2: `commaData` object, modkit parser, genome generalization             |
|0.3.0  |Phase 3: refactored annotation and sliding window functions                       |
|0.4.0  |Phase 4: differential methylation (`diffMethyl`, beta-binomial, methylKit wrapper)|
|0.5.0  |Phase 5: full visualization suite, tests, vignettes                               |
|1.0.0  |Bioconductor submission; all checks passing; full documentation                   |

-----

## 10. Out of Scope (for 1.0)

The following are scientifically interesting extensions but are explicitly deferred beyond 1.0 to keep the project tractable:

- Multi-species comparative methylomics (comparing methylomes across different organisms)
- Integration with transcriptomics data (correlating methylation changes with RNA-seq)
- Motif discovery (finding enriched sequence contexts around differentially methylated sites)
- Phage and plasmid methylation (non-chromosomal elements)
- Shiny interactive browser
- Python or command-line interface

-----

*Document version: 1.0 — written March 2026*
*Author: Carl Stone, with project management by Claude (Anthropic)*