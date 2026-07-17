# Working with Multiple Modification Types

## Introduction

Bacterial genomes frequently carry multiple simultaneous DNA methylation
marks:

- **6mA** (N6-methyladenine) — deposited by Dam methyltransferase at
  GATC motifs; functions in mismatch repair and cell-cycle regulation.
- **5mC** (5-methylcytosine) — deposited by Dcm methyltransferase at
  CCWGG motifs; roles in gene regulation are less well characterized in
  bacteria.
- **4mC** (N4-methylcytosine) — deposited by restriction-modification
  systems at various sequence motifs; protects against foreign DNA.

commaKit stores all modification types in a single `commaData` object,
enabling joint analysis, comparison, and visualization across
modification types.

This vignette demonstrates multi-modification workflows using the
built-in `comma_example_data` dataset, which contains both 6mA and 5mC
sites.

``` r

library(commaKit)
data(comma_example_data)
```

## Exploring Multiple Modification Types

### What modifications are present?

``` r

modTypes(comma_example_data)
#> [1] "5mC" "6mA"
```

### Count sites per modification type

``` r

table(siteInfo(comma_example_data)$mod_type)
#> 
#> 4mC 5mC 6mA 
#>   0 195 393
```

### Per-modification methylation summary

[`methylomeSummary()`](https://carl-stone.github.io/commaKit/reference/methylomeSummary.md)
reports stats per sample. Combine it with `mod_type` filtering to
compare distributions:

``` r

methylomeSummary(comma_example_data, mod_type = "6mA")[
  ,
  c("sample_name", "condition", "mean_beta", "median_beta", "n_covered")
]
#>   sample_name condition mean_beta median_beta n_covered
#> 1      ctrl_1   control 0.9036813   0.9173248       393
#> 2      ctrl_2   control 0.9029957   0.9139209       393
#> 3      ctrl_3   control 0.9014299   0.9133574       393
#> 4     treat_1 treatment 0.8557059   0.9107859       393
#> 5     treat_2 treatment 0.8491166   0.9039562       393
#> 6     treat_3 treatment 0.8527682   0.9074161       393
```

``` r

methylomeSummary(comma_example_data, mod_type = "5mC")[
  ,
  c("sample_name", "condition", "mean_beta", "median_beta", "n_covered")
]
#>   sample_name condition mean_beta median_beta n_covered
#> 1      ctrl_1   control 0.7885015   0.8101719       195
#> 2      ctrl_2   control 0.8052173   0.8270815       195
#> 3      ctrl_3   control 0.7879713   0.8214410       195
#> 4     treat_1 treatment 0.7956813   0.8264400       195
#> 5     treat_2 treatment 0.8122951   0.8483001       195
#> 6     treat_3 treatment 0.8107686   0.8296301       195
```

## Distribution Comparison Across Modification Types

[`plot_methylation_distribution()`](https://carl-stone.github.io/commaKit/reference/plot_methylation_distribution.md)
automatically facets by modification type when multiple types are
present:

``` r

plot_methylation_distribution(comma_example_data, per_sample = TRUE)
```

![Beta value density per sample, faceted by modification
type.](multiple-modification-types_files/figure-html/plot-dist-1.png)

Beta value density per sample, faceted by modification type.

6mA sites in bacteria commonly show a bimodal distribution (near 0 or
near 1) because restriction-modification systems methylate recognition
sequences with high efficiency. 5mC patterns can be more variable.

## Subsetting by Modification Type

Use
[`filterSites()`](https://carl-stone.github.io/commaKit/reference/filterSites.md)
to extract a single modification type:

``` r

only_6ma <- filterSites(comma_example_data, mod_type = "6mA")
only_5mc <- filterSites(comma_example_data, mod_type = "5mC")

cat("6mA sites:", nrow(methylation(only_6ma)), "\n")
#> 6mA sites: 393
cat("5mC sites:", nrow(methylation(only_5mc)), "\n")
#> 5mC sites: 195
```

## Differential Methylation: Testing Each Modification Type

### Test all types in one call

By default,
[`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
tests all modification types present in the object. Use
`results(..., mod_type = "6mA")` to extract type-specific results:

``` r

dm_all <- diffMethyl(comma_example_data, formula = ~condition)
res_6ma <- results(dm_all, mod_type = "6mA")
res_5mc <- results(dm_all, mod_type = "5mC")

cat(
  "6mA: significant sites (padj < 0.05, |Δβ| ≥ 0.2):",
  sum(!is.na(res_6ma$dm_padj) & res_6ma$dm_padj < 0.05 &
    abs(res_6ma$dm_delta_beta) >= 0.2, na.rm = TRUE), "\n"
)
#> 6mA: significant sites (padj < 0.05, |Δβ| ≥ 0.2): 31
cat(
  "5mC: significant sites (padj < 0.05, |Δβ| ≥ 0.2):",
  sum(!is.na(res_5mc$dm_padj) & res_5mc$dm_padj < 0.05 &
    abs(res_5mc$dm_delta_beta) >= 0.2, na.rm = TRUE), "\n"
)
#> 5mC: significant sites (padj < 0.05, |Δβ| ≥ 0.2): 11
```

### Volcano plots per modification type

``` r

p_6ma <- plot_volcano(res_6ma) +
  ggplot2::ggtitle("6mA differential methylation")
p_5mc <- plot_volcano(res_5mc) +
  ggplot2::ggtitle("5mC differential methylation")

# Arrange side by side if patchwork is available
if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(p_6ma, p_5mc, ncol = 1L)
} else {
  print(p_6ma)
  print(p_5mc)
}
```

![6mA differential methylation (top) and 5mC
(bottom).](multiple-modification-types_files/figure-html/volcanos-1.png)

6mA differential methylation (top) and 5mC (bottom).

## PCA Across Modification Types

[`plot_pca()`](https://carl-stone.github.io/commaKit/reference/plot_pca.md)
converts beta values to M-values via
[`mValues()`](https://carl-stone.github.io/commaKit/reference/mValues.md)
before PCA, which stabilizes variance for sites near 0 or 1. By default
all sites are used; pass `mod_type` to restrict to one modification
type:

``` r

plot_pca(comma_example_data, color_by = "condition")
```

![PCA using all modification
types.](multiple-modification-types_files/figure-html/pca-all-1.png)

PCA using all modification types.

``` r

plot_pca(comma_example_data, mod_type = "6mA", color_by = "condition")
```

![PCA using 6mA sites
only.](multiple-modification-types_files/figure-html/pca-6ma-1.png)

PCA using 6mA sites only.

Use `return_data = TRUE` to retrieve the PC scores for custom plots or
downstream analysis:

``` r

pca_df <- plot_pca(comma_example_data, mod_type = "6mA", return_data = TRUE)
# data.frame with PC1, PC2, sample_name, condition, replicate, ...
head(pca_df[, c("sample_name", "condition", "PC1", "PC2")])
#>   sample_name condition       PC1         PC2
#> 1      ctrl_1   control -15.40587  11.1961851
#> 2      ctrl_2   control -15.28274   5.0977438
#> 3      ctrl_3   control -16.11112 -14.9438670
#> 4     treat_1 treatment  13.82789  -8.7254793
#> 5     treat_2 treatment  15.15304   0.4254472
#> 6     treat_3 treatment  17.81880   6.9499703
attr(pca_df, "percentVar")
#>  PC1  PC2 
#> 44.4 15.1
```

## Genome Track with Both Modification Types

When `mod_type = NULL`,
[`plot_genome_track()`](https://carl-stone.github.io/commaKit/reference/plot_genome_track.md)
displays sites from all modification types, using different colors to
distinguish them:

``` r

plot_genome_track(comma_example_data,
  chromosome = "chr_sim",
  start = 1L, end = 30000L, annotation = FALSE
)
```

![Genome track showing 6mA (red) and 5mC (blue)
sites.](multiple-modification-types_files/figure-html/genome-track-1.png)

Genome track showing 6mA (red) and 5mC (blue) sites.

## Heatmaps per Modification Type

Filter the
[`results()`](https://carl-stone.github.io/commaKit/reference/results.md)
output to a single modification type before passing it to
[`plot_heatmap()`](https://carl-stone.github.io/commaKit/reference/plot_heatmap.md):

``` r

plot_heatmap(res_6ma, dm_all, n_sites = 20L)
```

![Heatmap of top 20 differentially methylated 6mA
sites.](multiple-modification-types_files/figure-html/heatmap-6ma-1.png)

Heatmap of top 20 differentially methylated 6mA sites.

## Session Information

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
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
#> [1] commaKit_0.2.0   BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] bitops_1.0-9                rlang_1.3.0                
#>   [3] magrittr_2.0.5              otel_0.2.0                 
#>   [5] matrixStats_1.5.0           compiler_4.5.2             
#>   [7] mgcv_1.9-4                  systemfonts_1.3.2          
#>   [9] vctrs_0.7.3                 reshape2_1.4.5             
#>  [11] stringr_1.6.0               pkgconfig_2.0.3            
#>  [13] crayon_1.5.3                fastmap_1.2.0              
#>  [15] XVector_0.50.0              labeling_0.4.3             
#>  [17] Rsamtools_2.26.0            rmarkdown_2.31             
#>  [19] UCSC.utils_1.6.1            ragg_1.5.2                 
#>  [21] xfun_0.60                   cachem_1.1.0               
#>  [23] cigarillo_1.0.0             GenomeInfoDb_1.46.2        
#>  [25] jsonlite_2.0.0              DelayedArray_0.36.1        
#>  [27] BiocParallel_1.44.0         parallel_4.5.2             
#>  [29] R6_2.6.1                    bslib_0.11.0               
#>  [31] stringi_1.8.7               RColorBrewer_1.1-3         
#>  [33] limma_3.66.0                rtracklayer_1.70.1         
#>  [35] GenomicRanges_1.62.1        jquerylib_0.1.4            
#>  [37] numDeriv_2016.8-1.1         Rcpp_1.1.2                 
#>  [39] Seqinfo_1.0.0               bookdown_0.47              
#>  [41] SummarizedExperiment_1.40.0 knitr_1.51                 
#>  [43] zoo_1.8-15                  R.utils_2.13.0             
#>  [45] IRanges_2.44.0              Matrix_1.7-4               
#>  [47] splines_4.5.2               tidyselect_1.2.1           
#>  [49] qvalue_2.42.0               abind_1.4-8                
#>  [51] yaml_2.3.12                 codetools_0.2-20           
#>  [53] curl_7.1.0                  lattice_0.22-7             
#>  [55] tibble_3.3.1                plyr_1.8.9                 
#>  [57] Biobase_2.70.0              withr_3.0.3                
#>  [59] S7_0.2.2                    coda_0.19-4.1              
#>  [61] evaluate_1.0.5              desc_1.4.3                 
#>  [63] mclust_6.1.3                Biostrings_2.78.0          
#>  [65] pillar_1.11.1               BiocManager_1.30.27        
#>  [67] MatrixGenerics_1.22.0       stats4_4.5.2               
#>  [69] generics_0.1.4              RCurl_1.98-1.19            
#>  [71] emdbook_1.3.14              S4Vectors_0.48.1           
#>  [73] ggplot2_4.0.3               scales_1.4.0               
#>  [75] gtools_3.9.5                glue_1.8.1                 
#>  [77] tools_4.5.2                 BiocIO_1.20.0              
#>  [79] data.table_1.18.4           GenomicAlignments_1.46.0   
#>  [81] fs_2.1.0                    mvtnorm_1.4-2              
#>  [83] XML_3.99-0.23               grid_4.5.2                 
#>  [85] bbmle_1.0.25.1              bdsmatrix_1.3-7            
#>  [87] patchwork_1.3.2             nlme_3.1-168               
#>  [89] restfulr_0.0.17             cli_3.6.6                  
#>  [91] textshaping_1.0.5           fastseg_1.56.0             
#>  [93] S4Arrays_1.10.1             methylKit_1.36.0           
#>  [95] dplyr_1.2.1                 gtable_0.3.6               
#>  [97] R.methodsS3_1.8.2           sass_0.4.10                
#>  [99] digest_0.6.39               BiocGenerics_0.56.0        
#> [101] SparseArray_1.10.10         rjson_0.2.23               
#> [103] htmlwidgets_1.6.4           farver_2.1.2               
#> [105] htmltools_0.5.9             pkgdown_2.2.1              
#> [107] R.oo_1.27.1                 lifecycle_1.0.5            
#> [109] httr_1.4.8                  statmod_1.5.2              
#> [111] MASS_7.3-65
```
