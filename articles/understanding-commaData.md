# Understanding the commaData object

## Why `commaData` matters

Most commaKit workflow functions start from a `commaData` object, and
many return an updated `commaData` object with added analysis layers.
Other functions expose those layers as tables, plots, or enrichment
result lists. The object is the stable substrate for the workflow: rows
are genomic methylation sites, columns are samples, assays store
measurements, and later analyses add named layers to the same biological
coordinates.

``` r

library(commaKit)
data(comma_example_data)

comma_example_data
#> class: commaData
#> sites: 588 | samples: 6
#> mod types: 5mC, 6mA
#> motifs: CCWGG, GATC
#> mod contexts: 5mC_CCWGG, 6mA_GATC
#> conditions: control, treatment
#> genome: 1 chromosome (100,000 bp total)
#> annotation: 5 features
#> motif sites: none
#> min_coverage: 5
```

The built-in example data are synthetic, but they exercise the same
object model used for real data. They contain multiple modification
types and multiple sequence contexts, so the examples below avoid
assuming that a data set is specific to one methylation chemistry.

## Object anatomy

`commaData` extends Bioconductor’s `RangedSummarizedExperiment`. The
core dimensions therefore have the same interpretation throughout the
package:

- rows are one-base methylation sites;
- columns are biological samples;
- assays are site-by-sample matrices;
- `rowRanges()` stores genomic coordinates;
- `colData()` stores sample metadata;
- `metadata()` stores package-level context such as annotations, motif
  sites, caller information, and analysis parameters.

``` r

dim(comma_example_data)
#> [1] 588   6
SummarizedExperiment::assayNames(comma_example_data)
#> [1] "methylation"      "coverage"         "mod_counts"       "canonical_counts"
```

The main package accessors expose the two required assay matrices. Beta
values are proportions in `[0, 1]`; coverage values are read depths.

``` r

methylation(comma_example_data)[1:5, 1:3]
#>         ctrl_1    ctrl_2    ctrl_3
#> [1,] 0.7575443 0.8107216 0.8870622
#> [2,] 0.9738114 0.9121264 0.8963508
#> [3,] 0.8678683 0.8581721 0.9193607
#> [4,] 0.9569443 0.9232934 0.8870824
#> [5,] 0.9165347 0.8104115 0.9441796
siteCoverage(comma_example_data)[1:5, 1:3]
#>      ctrl_1 ctrl_2 ctrl_3
#> [1,]     76     93     51
#> [2,]    124     83    126
#> [3,]     88    149     55
#> [4,]     31     37     15
#> [5,]     73    130    110
```

Sample metadata live in `colData()` and are also available through
[`sampleInfo()`](https://carl-stone.github.io/commaKit/reference/sampleInfo.md).

``` r

sampleInfo(comma_example_data)
#>         sample_name condition replicate
#> ctrl_1       ctrl_1   control         1
#> ctrl_2       ctrl_2   control         2
#> ctrl_3       ctrl_3   control         3
#> treat_1     treat_1 treatment         1
#> treat_2     treat_2 treatment         2
#> treat_3     treat_3 treatment         3
```

## Genomic row identity

Schema v2 stores genomic identity in `rowRanges()`, not in string row
names. Each row is a one-base `GRanges` range with methylation-site
metadata in the range metadata columns.

``` r

rr <- SummarizedExperiment::rowRanges(comma_example_data)
rr[1:4]
#> GRanges object with 4 ranges and 3 metadata columns:
#>       seqnames    ranges strand | mod_type       motif   is_diff
#>          <Rle> <IRanges>  <Rle> | <factor> <character> <logical>
#>   [1]  chr_sim       443      + |      6mA        GATC     FALSE
#>   [2]  chr_sim       512      + |      6mA        GATC     FALSE
#>   [3]  chr_sim      1024      + |      6mA        GATC     FALSE
#>   [4]  chr_sim      1073      + |      6mA        GATC     FALSE
#>   -------
#>   seqinfo: 1 sequence (1 circular) from an unspecified genome
GenomeInfoDb::seqinfo(rr)
#> Seqinfo object with 1 sequence (1 circular) from an unspecified genome:
#>   seqnames seqlengths isCircular genome
#>   chr_sim      100000       TRUE   <NA>
```

For human-readable inspection,
[`siteInfo()`](https://carl-stone.github.io/commaKit/reference/siteInfo.md)
reconstructs a flat table from `rowRanges()` and its metadata columns.

``` r

si <- siteInfo(comma_example_data)
si[1:5, c(
  "chrom", "position", "strand", "mod_type", "motif",
  "mod_context", "site_key"
)]
#> DataFrame with 5 rows and 7 columns
#>         chrom  position      strand mod_type       motif mod_context
#>   <character> <integer> <character> <factor> <character> <character>
#> 1     chr_sim       443           +      6mA        GATC    6mA_GATC
#> 2     chr_sim       512           +      6mA        GATC    6mA_GATC
#> 3     chr_sim      1024           +      6mA        GATC    6mA_GATC
#> 4     chr_sim      1073           +      6mA        GATC    6mA_GATC
#> 5     chr_sim      1536           -      6mA        GATC    6mA_GATC
#>                 site_key
#>              <character>
#> 1 chr_sim:443:+:6mA:GATC
#> 2 chr_sim:512:+:6mA:GATC
#> 3 chr_sim:1024:+:6mA:G..
#> 4 chr_sim:1073:+:6mA:G..
#> 5 chr_sim:1536:-:6mA:G..
```

The `site_key` column is a display label. Alignment between objects and
samples is based on genomic ranges and site metadata, not on matching
string keys.

## Modification type, context, and `mod_context`

The `mod_context` value should be read as **modification + context**. It
is not just a synonym for motif.

- `mod_type` answers which chemical or base modification is measured,
  such as `6mA`, `5mC`, or `4mC`.
- `motif` answers which sequence or biological context the site belongs
  to, such as `GATC`, `CCWGG`, `CpG`, or a future context definition.
- `mod_context` answers which biological class of sites is being
  analyzed, such as `6mA_GATC` or `5mC_CCWGG`.

``` r

modTypes(comma_example_data)
#> [1] "5mC" "6mA"
motifs(comma_example_data)
#> [1] "CCWGG" "GATC"
modContexts(comma_example_data)
#> [1] "5mC_CCWGG" "6mA_GATC"
```

The built-in data illustrate that a `mod_context` is the combination of
both parts.

``` r

context_table <- unique(as.data.frame(si[, c("mod_type", "motif", "mod_context")]))
context_table[order(context_table$mod_context), ]
#>     mod_type motif mod_context
#> 394      5mC CCWGG   5mC_CCWGG
#> 1        6mA  GATC    6mA_GATC
```

When a caller does not provide motif information, `motif` may be `NA`
and the computed `mod_context` falls back to the modification type
alone. Code that uses
[`modContexts()`](https://carl-stone.github.io/commaKit/reference/modContexts.md)
rather than hard-coded labels remains compatible with single-context,
multi-context, and future-context data sets.

This distinction is important for differential methylation. Analyses
should not silently pool distinct `mod_context` groups, because
methylation at different modification-plus-context units can reflect
different biological systems.

## Metadata and package-level context

Object-level context is stored in `metadata()`. In this example,
annotation and motif-site containers are available without reading any
external files.

``` r

md <- S4Vectors::metadata(comma_example_data)
names(md)
#> [1] "annotation"       "motifSites"       "min_coverage"     "assay_defaults"  
#> [5] "assay_provenance"

c(
  annotation_features = length(annotation(comma_example_data)),
  motif_sites = length(motifSites(comma_example_data))
)
#> annotation_features         motif_sites 
#>                   5                   0

vapply(md, function(x) paste(class(x), collapse = "/"), character(1))
#>       annotation       motifSites     min_coverage   assay_defaults 
#>        "GRanges"        "GRanges"        "integer"           "list" 
#> assay_provenance 
#>           "list"
```

## Analysis layers

A useful way to learn commaKit is to think in layers:

- the measurement layer contains methylation and coverage assays;
- the coordinate layer is `rowRanges()`;
- the sample layer is `colData()`;
- the annotation layer adds feature relationships to sites;
- the statistical layer adds differential methylation results;
- plotting, filtering, and enrichment functions interpret those layers.

The object dimensions stay anchored to the same sites and samples while
layers are added or recomputed. The v1 policy is deliberately small and
strong:

- **Raw evidence assays are canonical.** `methylation`, `coverage`,
  `mod_counts`, `canonical_counts`, and `other_mod_counts` store caller
  or reconstructed evidence. commaKit functions do not mutate these raw
  assays in place. Users can still edit assays manually with standard
  `SummarizedExperiment` tools, but package-level transformations should
  not silently rewrite the raw evidence layer.
- **Transformations are named layers.** Normalized beta values,
  M-values, or future Seurat-like transformations should be added as
  explicit assay layers with provenance and parent assays, rather than
  replacing the raw assays.
- **Filtering returns subset objects.**
  [`filterSites()`](https://carl-stone.github.io/commaKit/reference/filterSites.md)
  and `[` follow Bioconductor subsetting semantics: they return a
  `commaData` object with assay layers subset to the requested rows or
  samples. Site-level result layers are trimmed for row filters; sample
  filters keep those site-level results attached to the remaining sites.
  commaKit does not create hidden lazy views.

``` r

base_cols <- colnames(siteInfo(comma_example_data))

annotated <- annotateSites(comma_example_data, keep = "overlap")
annotation_cols <- setdiff(colnames(siteInfo(annotated)), base_cols)
annotation_cols
#> [1] "feature_types" "feature_names"

identical(dim(comma_example_data), dim(annotated))
#> [1] TRUE
```

Layering is intended to be explicit: rerunning an unnamed compatibility
layer updates the active/default view, while explicitly named layers can
coexist for comparison. In either case, commaKit never creates a second
copy of the biological rows.

``` r

annotated_again <- annotateSites(annotated, keep = "overlap")
identical(colnames(siteInfo(annotated)), colnames(siteInfo(annotated_again)))
#> [1] TRUE
identical(dim(annotated), dim(annotated_again))
#> [1] TRUE
```

Differential methylation adds a statistical result layer. By default,
[`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
works over the modification contexts present in the object, which keeps
distinct modification-plus-context units separate.

``` r

dm <- diffMethyl(annotated, formula = ~condition, method = "quasi_f")
dm_cols <- setdiff(colnames(siteInfo(dm)), colnames(siteInfo(annotated)))
dm_cols
#> [1] "dm_pvalue"              "dm_padj"                "dm_delta_beta"         
#> [4] "dm_mean_beta_control"   "dm_mean_beta_treatment"

S4Vectors::metadata(dm)$diffMethyl_params[c("method", "p_adjust_method")]
#> $method
#> [1] "quasi_f"
#> 
#> $p_adjust_method
#> [1] "BH"
resultLayers(dm)
#> DataFrame with 1 row and 18 columns
#>          name        role                   type      source is_default
#>   <character> <character>            <character> <character>  <logical>
#> 1  diffMethyl  diffMethyl differential_methyla..  diffMethyl       TRUE
#>        method     formula   reference   treatment     mod_context
#>   <character> <character> <character> <character> <CharacterList>
#> 1     quasi_f  ~condition     control   treatment                
#>          mod_type           motif p_adjust_method min_coverage     alpha
#>   <CharacterList> <CharacterList>     <character>    <integer> <numeric>
#> 1                                              BH            5       0.5
#>                           result_cols              timestamp package_version
#>                       <CharacterList>            <character>     <character>
#> 1 dm_pvalue,dm_padj,dm_delta_beta,... 2026-08-26 19:06:09 ..           0.2.0
```

[`results()`](https://carl-stone.github.io/commaKit/reference/results.md)
exposes the active result layer as a tidy data frame while preserving
the site metadata needed to interpret each row. When you keep multiple
named runs, pass `result` or `result_name` to retrieve an older layer.

``` r

res <- results(dm)
res <- res[order(res$dm_padj), ]
head(res[, c(
  "chrom", "position", "strand", "mod_type", "motif",
  "mod_context", "dm_padj", "dm_delta_beta"
)])
#>       chrom position strand mod_type motif mod_context      dm_padj
#> 64  chr_sim    16504      +      6mA  GATC    6mA_GATC 6.309184e-07
#> 196 chr_sim    50176      -      6mA  GATC    6mA_GATC 6.309184e-07
#> 287 chr_sim    70003      -      6mA  GATC    6mA_GATC 6.309184e-07
#> 347 chr_sim    86016      +      6mA  GATC    6mA_GATC 6.309184e-07
#> 249 chr_sim    61440      +      6mA  GATC    6mA_GATC 8.942641e-07
#> 63  chr_sim    16384      -      6mA  GATC    6mA_GATC 1.248983e-06
#>     dm_delta_beta
#> 64     -0.6142911
#> 196    -0.7336497
#> 287    -0.7050844
#> 347    -0.6743832
#> 249    -0.7090099
#> 63     -0.7178796
```

The important pattern is that downstream functions do not replace the
central object model. They add, expose, or interpret layers attached to
the same genomic methylation-site substrate.

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Etc/UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] commaKit_0.2.0   BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] SummarizedExperiment_1.42.0 gtable_0.3.6               
#>  [3] xfun_0.60                   bslib_0.12.0               
#>  [5] ggplot2_4.0.3               htmlwidgets_1.6.4          
#>  [7] Biobase_2.72.0              lattice_0.22-9             
#>  [9] vctrs_0.7.3                 tools_4.6.1                
#> [11] generics_0.1.4              stats4_4.6.1               
#> [13] tibble_3.3.1                pkgconfig_2.0.3            
#> [15] Matrix_1.7-6                RColorBrewer_1.1-3         
#> [17] S7_0.2.2                    desc_1.4.3                 
#> [19] S4Vectors_0.50.2            lifecycle_1.0.5            
#> [21] compiler_4.6.1              farver_2.1.2               
#> [23] textshaping_1.0.5           statmod_1.5.2              
#> [25] Seqinfo_1.2.0               GenomeInfoDb_1.48.0        
#> [27] htmltools_0.5.9             sass_0.4.10                
#> [29] yaml_2.3.12                 pkgdown_2.2.1              
#> [31] pillar_1.11.1               jquerylib_0.1.4            
#> [33] DelayedArray_0.38.2         cachem_1.1.0               
#> [35] limma_3.68.5                abind_1.4-8                
#> [37] tidyselect_1.2.1            digest_0.6.39              
#> [39] dplyr_1.2.1                 bookdown_0.47              
#> [41] fastmap_1.2.0               grid_4.6.1                 
#> [43] cli_3.6.6                   SparseArray_1.12.2         
#> [45] magrittr_2.0.5              S4Arrays_1.12.0            
#> [47] UCSC.utils_1.8.0            scales_1.4.0               
#> [49] rmarkdown_2.31              XVector_0.52.0             
#> [51] httr_1.4.8                  matrixStats_1.5.0          
#> [53] otel_0.2.0                  ragg_1.5.2                 
#> [55] zoo_1.9-0                   evaluate_1.0.5             
#> [57] knitr_1.51                  GenomicRanges_1.64.0       
#> [59] IRanges_2.46.0              rlang_1.3.0                
#> [61] glue_1.8.1                  BiocManager_1.30.27        
#> [63] BiocGenerics_0.58.1         jsonlite_2.0.0             
#> [65] R6_2.6.1                    MatrixGenerics_1.24.0      
#> [67] systemfonts_1.3.2           fs_2.1.0
```
