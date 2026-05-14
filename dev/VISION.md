# comma Vision — The Dream Package

**What comma would look like if it were as good as DESeq2 for bacterial methylation.**

This document describes the aspirational target. Not everything here ships in v1.0 — see PRD.md for the shippable milestone. But the vision informs the architecture: build nothing that conflicts with these goals.

---

## 1. What makes DESeq2 the gold standard

DESeq2 succeeded because it solved five problems at once:

1. **One data container** — `DESeqDataSet` holds raw counts, design, size factors, dispersions, and results. No loose data frames.
2. **One workflow** — `DESeqDataSetFromMatrix()` → `DESeq()` → `results()` → `lfcShrink()`. Three function calls, one object.
3. **Shrinkage** — Empirical Bayes priors stabilize dispersion and fold-change estimates. Small samples don't produce garbage.
4. **Rich diagnostics** — `plotPCA`, `plotMA`, `plotCounts`, dispersion plots, Cook's distance. You can see what the model is doing.
5. **Extensibility** — Complex designs, custom contrasts, interaction terms, batch correction. The object grows with your analysis.

comma should match this standard for methylation analysis.

---

## 2. The Dream Workflow

```r
# One object, one pipeline
cd <- commaData(files, colData, genome = "ecoli.fa",
                annotation = "ecoli.gff3",
                expected_mod_contexts = list("6mA" = c("GATC", "ACCACC"),
                                            "5mC" = "CCWGG"))

# QC — automatic, returns a summary you can cite
cd <- commaQC(cd)                    # runs all QC checks, stores results
qcReport(cd)                         # printable QC report

# Annotation — already best-in-class
cd <- annotateSites(cd, keep = "metagene")

# Differential methylation — shrinkage, like DESeq2
cd <- diffMethyl(cd, ~ condition)     # default: beta-binomial with shrinkage
cd <- diffMethyl(cd, ~ condition + batch, method = "limma")  # complex designs
res <- results(cd, contrast = c("condition", "treat", "ctrl"))
res <- lfcShrink(cd, res)            # shrink delta_beta estimates

# DMR calling — region-level, not just site-level
dmr <- callDMR(cd, method = "bsseq")  # or "dss" or "comb-p"

# Enrichment — gene-level, with proper universe
enr <- enrichMethylation(cd, ont = "BP", gene_role = "both")

# Visualization — every plot is publication-ready
plot_volcano(res, label_top = 10)
plot_metagene(cd, feature = "gene", color_by = "condition")
plot_tss_profile(cd, regulatory = TRUE)
plot_genome_track(cd, "chr", 1e6, 2e6, tracks = c("methylation", "annotation", "dmr"))
plot_manhattan(res)                  # genome-wide DM landscape
plot_effect_size(res)                # shrinkage visualization

# Export — everything is reproducible
exportResults(cd, dir = "results/")  # BED, CSV, and RDS
writeBED(cd, "dm_sites.bed")
```

---

## 3. Feature Roadmap — Dream comma

### Tier 1: Essential for best-in-class (v1.1–v1.2)

| Feature | Why it matters | Implementation notes |
|---------|---------------|----------------------|
| **Effect size shrinkage** (`lfcShrink()` equivalent) | Small-sample delta_beta estimates are noisy. Shrinkage stabilizes ranking and visualization. | Empirical Bayes prior on delta_beta, like DESeq2's apeglm/ashr. Could use beta-binomial posterior. |
| **DMR calling** (region-level differential methylation) | Site-level testing misses coordinated methylation changes across regions. Bacterial genomes have small genes — DMRs often span entire operons. | Interface to bsseq/DSS DMR methods, or custom sliding-window with cluster-based p-value combination (comb-p/Stouffer). |
| **Batch effect correction** | Real experiments have batch structure. Currently `diffMethyl()` only handles simple `~ condition` designs. | Support `~ batch + condition` in all backends. limma backend already handles this; extend to beta-binomial. |
| **Comprehensive QC report** | Users need a single artifact they can cite in methods sections. methylKit has `getMethylationStats()` but nothing integrated. | `commaQC()` runs all QC checks, stores results in `metadata(cd)$qc`. `qcReport()` generates a printable summary. |
| **Vignette: enrichment workflow** | The enrichment pipeline is the most complex part of comma and has no vignette coverage. | Extend getting-started vignette or add a third vignette. |

### Tier 2: Differentiating features (v1.3–v2.0)

| Feature | Why it matters | Implementation notes |
|---------|---------------|----------------------|
| **Manhattan plot** (`plot_manhattan()`) | Genome-wide view of DM landscape. Standard in GWAS; missing in methylation tools. | Chromosome on x, -log10(padj) on y, color by mod_context. |
| **Effect size visualization** (`plot_effect_size()`) | Show shrinkage effect: raw vs shrunk delta_beta. Builds trust in the model. | Like DESeq2's MA plot with shrinkage overlay. |
| **Multi-factor designs** | Time-course, strain × condition, dose-response. | Extend formula interface; ensure all backends handle >2 conditions. |
| **Variance-stabilizing transformation** | M-values are heteroscedastic at extremes. VST/rlog-style transform improves PCA and clustering. | Arcsinh-square-root or beta-binomial VST. |
| **Independent hypothesis weighting** (IHW) | More power than BH adjustment when covariates (coverage, methylation level) predict power. | Interface to IHW package, like DESeq2. |
| **Per-read methylation analysis** | Nanopore gives per-read calls. Aggregating to per-site loses information about cell-state heterogeneity. | Store per-read data in a separate slot; provide aggregation methods. |
| **Operon-aware annotation** | Bacterial genes are organized in operons. DM at one site in an operon affects the whole transcription unit. | `annotateSites()` gains `group_by = "operon"` mode; enrichment at operon level. |
| **Circular genome plots** | Bacterial chromosomes are circular. Linear genome tracks miss the origin/terminus junction. | `plot_circos()` or `plot_genome_track(circular = TRUE)`. |
| **Reproducible pipeline** (targets integration) | Users want `tar_make()` to run the full workflow. | Provide a `targets.R` template or `comma_pipeline()` helper. |

### Tier 3: Long-term vision (v2.0+)

| Feature | Why it matters | Implementation notes |
|---------|---------------|----------------------|
| **Multi-species comparative methylomics** | Compare methylation patterns across related strains/species. | Synteny-based alignment of methylation maps. |
| **Integration with transcriptomics** | Correlate methylation changes with gene expression. | Interface to DESeq2 results; joint visualization. |
| **Haplotype-resolved methylation** | Nanopore reads span haplotypes. Methylation can differ between alleles. | Requires variant calling + per-haplotype methylation aggregation. |
| **Cell-type deconvolution** | Mixed bacterial populations have different methylation patterns per cell type. | Like EpiDISH for eukaryotes; per-read methylation as cell-type signal. |
| **Shiny browser** | Interactive exploration for non-programmers. | `commaExplorer()` Shiny app. |
| **Python interface** | Python users in microbiology. | reticulate wrapper or separate Python package. |
| **Single-cell methylation** | Nanopore can resolve single-cell methylation (scNanoCOOLSeq, etc.). | New data container extending commaData. |

---

## 4. Competitive Landscape

| Feature | comma | methylKit | DSS | nanomethyR | DiffMethylTools |
|---------|-------|-----------|-----|-----------|-----------------|
| **Data container** | commaData (S4/SE) | methylRaw/methylBase (S4) | BSseq (S4/SE) | ModBamResult (S4/SE) | None (data.frames) |
| **Input format** | modkit, Dorado, Megalodon | Bismark BAM | Bismark | modkit/Dorado BAM | Various |
| **Modification types** | 6mA, 5mC, 4mC | 5mC only | 5mC only | 5mC only | 5mC only |
| **mod_context** | First-class | No | No | No | No |
| **Bacterial genomes** | Native | No (eukaryotic) | No | No | No |
| **DM backends** | 4 (beta-binomial, methylKit, limma, quasi_f) | 2 (Fisher, logistic) | 1 (beta-binomial Wald) | 1 (beta-binomial) | Multiple |
| **Effect size shrinkage** | Planned | No | No | No | No |
| **DMR calling** | Planned | Tile-based | Yes (DML→DMR) | No | Yes |
| **Enrichment** | Built-in (GO/KEGG) | External (genomation) | No | No | Built-in |
| **Visualization** | 8 plot functions | 3 (histogram, clustering, PCA) | 1 (heatmap) | 1 (genome track) | Multiple |
| **QC report** | Planned | Basic stats | No | No | No |
| **Vignettes** | 2 | 1 | 1 | 1 | 1 |
| **Bioconductor** | In progress | Yes (since 2016) | Yes (since 2013) | Yes (since 2023) | No |

**comma's unique advantages:**
1. Only package that handles 6mA, 5mC, and 4mC simultaneously
2. Only package with mod_context as a first-class variable
3. Only package designed for bacterial (monoploid, circular) genomes
4. Only package with native Nanopore format support (modkit, Dorado)
5. Only package with built-in enrichment (no external genomation dependency)
6. Most visualization functions of any methylation package
7. Most DM backends (4 vs 1–2)

**comma's gaps vs competition:**
1. No DMR calling (DSS has it, methylKit has tile-based)
2. No effect size shrinkage (no package has it — opportunity)
3. No batch effect handling (DESeq2 has it)
4. Not yet on Bioconductor (credibility gap)
5. No per-read analysis (nanomethyR has some support)

---

## 5. Architecture Principles

These principles should guide all development, from v1.0 through the dream:

1. **commaData is the single source of truth.** Every function takes it and returns it (or a tidy data frame derived from it). No loose data frames floating around.
2. **mod_context is never optional.** Every analysis must respect that 6mA_GATC and 6mA_ACCACC are biologically distinct. Mixing them is a scientific error.
3. **Effect sizes on beta scale.** Delta_beta (0–1) is interpretable. M-values are for internal computation only.
4. **Multiple testing is genome-wide.** Correction across all mod_contexts, not per-context.
5. **Shrinkage before reporting.** Raw estimates are noisy; shrinkage stabilizes ranking and visualization.
6. **Bacterial genomes are circular.** Sliding windows, annotation, and genome tracks must handle origin/terminus.
7. **Every plot returns a ggplot.** Users can customize further. No side-effect plotting.
8. **Bioconductor idioms.** S4 classes, explicit imports, roxygen2, testthat ed3. No tidyverse in Imports.

---

## 6. Success Metrics

How we know comma has arrived:

1. **Bioconductor acceptance** — on the release branch
2. **Citation count** — cited in ≥5 published papers within 2 years of v1.0
3. **Community adoption** — ≥100 GitHub stars, questions on Bioc Support site
4. **Lab productivity** — Carl's lab and collaborators can go from data to results in <1 day
5. **Feature parity** — matches or exceeds methylKit on every dimension that matters for bacterial methylation
