# Package index

## All functions

- [`annotateSites()`](https://carl-stone.github.io/commaKit/reference/annotateSites.md)
  : Annotate methylation sites relative to genomic features
- [`annotation(`*`<commaData>`*`)`](https://carl-stone.github.io/commaKit/reference/annotation-commaData-method.md)
  : Accessor for genomic feature annotation
- [`assayLayers()`](https://carl-stone.github.io/commaKit/reference/assayLayers.md)
  : List assay layers and defaults
- [`assayProvenance()`](https://carl-stone.github.io/commaKit/reference/assayProvenance.md)
  : Accessor for assay provenance metadata
- [`buildKEGGGeneIDMap()`](https://carl-stone.github.io/commaKit/reference/buildKEGGGeneIDMap.md)
  : Build a KEGG gene ID map for symbol translation
- [`buildKEGGTermGene()`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md)
  : Build a KEGG term-to-gene mapping for use with enrichMethylation()
- [`canonicalCounts()`](https://carl-stone.github.io/commaKit/reference/canonicalCounts.md)
  : Accessor for observed canonical-read counts
- [`commaData-class`](https://carl-stone.github.io/commaKit/reference/commaData-class.md)
  : commaData: the central data object for the commaKit package
- [`commaData()`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  : Create a commaData object from methylation calling output files
- [`commaKit-package`](https://carl-stone.github.io/commaKit/reference/commaKit-package.md)
  [`commaKit`](https://carl-stone.github.io/commaKit/reference/commaKit-package.md)
  [`comma`](https://carl-stone.github.io/commaKit/reference/commaKit-package.md)
  : commaKit: Comparative Microbial Methylomics Analysis Kit
- [`comma_example_data`](https://carl-stone.github.io/commaKit/reference/comma_example_data.md)
  : Synthetic example methylation dataset for the commaKit package
- [`coverage(`*`<commaData>`*`)`](https://carl-stone.github.io/commaKit/reference/coverage-commaData-method.md)
  : Deprecated coverage accessor for commaData objects
- [`coverageDepth()`](https://carl-stone.github.io/commaKit/reference/coverageDepth.md)
  : Windowed sequencing depth across the genome
- [`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
  : Identify differentially methylated sites between conditions
- [`enrichMethylation()`](https://carl-stone.github.io/commaKit/reference/enrichMethylation.md)
  : Gene set enrichment analysis of differential methylation results
- [`filterResults()`](https://carl-stone.github.io/commaKit/reference/filterResults.md)
  : Filter differential methylation results by significance thresholds
- [`filterSites()`](https://carl-stone.github.io/commaKit/reference/filterSites.md)
  : Filter a commaData object by condition, modification type, or
  chromosome
- [`findMotifSites()`](https://carl-stone.github.io/commaKit/reference/findMotifSites.md)
  : Find all instances of a sequence motif in a genome
- [`genome(`*`<commaData>`*`)`](https://carl-stone.github.io/commaKit/reference/genome-commaData-method.md)
  : GenomeInfoDb genome accessor compatibility method for commaData
- [`genomeSizes()`](https://carl-stone.github.io/commaKit/reference/genomeSizes.md)
  : Accessor for chromosome size information
- [`loadAnnotation()`](https://carl-stone.github.io/commaKit/reference/loadAnnotation.md)
  : Load genomic feature annotations from a GFF3 or BED file
- [`mValues()`](https://carl-stone.github.io/commaKit/reference/mValues.md)
  : Compute M-values from a commaData object
- [`methylation()`](https://carl-stone.github.io/commaKit/reference/methylation.md)
  : Accessor for the methylation (beta value) matrix
- [`methylomeSummary()`](https://carl-stone.github.io/commaKit/reference/methylomeSummary.md)
  : Summarize per-sample methylation and coverage distributions
- [`minCoverage()`](https://carl-stone.github.io/commaKit/reference/minCoverage.md)
  : Accessor for the minimum coverage threshold
- [`modContexts()`](https://carl-stone.github.io/commaKit/reference/modContexts.md)
  : Return the modification contexts present in a commaData object
- [`modCounts()`](https://carl-stone.github.io/commaKit/reference/modCounts.md)
  : Accessor for observed modified-read counts
- [`modTypes()`](https://carl-stone.github.io/commaKit/reference/modTypes.md)
  : Return the modification types present in a commaData object
- [`motifSites()`](https://carl-stone.github.io/commaKit/reference/motifSites.md)
  : Accessor for motif site positions
- [`motifs()`](https://carl-stone.github.io/commaKit/reference/motifs.md)
  : Accessor for sequence context motifs present in a commaData object
- [`otherModCounts()`](https://carl-stone.github.io/commaKit/reference/otherModCounts.md)
  : Accessor for observed non-target modified-read counts
- [`plot_coverage()`](https://carl-stone.github.io/commaKit/reference/plot_coverage.md)
  : Plot coverage depth distribution
- [`plot_genome_track()`](https://carl-stone.github.io/commaKit/reference/plot_genome_track.md)
  : Genome browser-style methylation track plot
- [`plot_heatmap()`](https://carl-stone.github.io/commaKit/reference/plot_heatmap.md)
  : Heatmap of top differentially methylated sites
- [`plot_metagene()`](https://carl-stone.github.io/commaKit/reference/plot_metagene.md)
  : Metagene plot of methylation across genomic features
- [`plot_methylation_distribution()`](https://carl-stone.github.io/commaKit/reference/plot_methylation_distribution.md)
  : Plot methylation beta value distributions
- [`plot_pca()`](https://carl-stone.github.io/commaKit/reference/plot_pca.md)
  : PCA of methylation profiles
- [`plot_tss_profile()`](https://carl-stone.github.io/commaKit/reference/plot_tss_profile.md)
  : TSS-centered methylation profile
- [`plot_volcano()`](https://carl-stone.github.io/commaKit/reference/plot_volcano.md)
  : Volcano plot for differential methylation results
- [`resultLayers()`](https://carl-stone.github.io/commaKit/reference/resultLayers.md)
  : List differential methylation result layers
- [`results()`](https://carl-stone.github.io/commaKit/reference/results.md)
  : Extract differential methylation results as a tidy data frame
- [`sampleInfo()`](https://carl-stone.github.io/commaKit/reference/sampleInfo.md)
  : Accessor for per-sample metadata
- [`siteCoverage()`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md)
  : Accessor for the sequencing coverage (read depth) matrix
- [`siteInfo()`](https://carl-stone.github.io/commaKit/reference/siteInfo.md)
  : Accessor for per-site metadata
- [`slidingWindow()`](https://carl-stone.github.io/commaKit/reference/slidingWindow.md)
  : Sliding window methylation summary along the genome
- [`` `[`( ``*`<commaData>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://carl-stone.github.io/commaKit/reference/sub-commaData-ANY-ANY-ANY-method.md)
  : Subset a commaData object by sites and/or samples
- [`subset(`*`<commaData>`*`)`](https://carl-stone.github.io/commaKit/reference/subset.commaData.md)
  : Deprecated subset method for commaData objects
- [`summarizeRegions()`](https://carl-stone.github.io/commaKit/reference/summarizeRegions.md)
  : Summarize methylation counts over genomic regions
- [`varianceByDepth()`](https://carl-stone.github.io/commaKit/reference/varianceByDepth.md)
  : Methylation variance as a function of sequencing depth
- [`writeBED()`](https://carl-stone.github.io/commaKit/reference/writeBED.md)
  : Export methylation data as a BED file
