# Troubleshooting Data Import

## Overview

Most import problems come from one of four places:

1.  The file is not modkit pileup bedMethyl output.
2.  Sample names in `files` and `colData$sample_name` do not match
    exactly.
3.  Genome information is missing, unnamed, or uses chromosome names
    that do not match the methylation files.
4.  Annotation files are missing optional dependencies or use an
    unexpected file extension.

Start by checking each input file and confirming that the sample names
are identical across `files` and `colData`.

``` r

library(commaKit)

files <- c(
  ctrl_1 = "ctrl_1.modkit.bed",
  treat_1 = "treat_1.modkit.bed"
)

col_data <- data.frame(
  sample_name = c("ctrl_1", "treat_1"),
  condition = c("control", "treatment"),
  replicate = c(1L, 1L)
)

obj <- commaData(
  files = files,
  colData = col_data,
  genome = c(chr1 = 4641652L)
)
```

## Modkit BED Format Problems

commaKit expects modkit pileup bedMethyl output with at least 18 columns
and no header row. The parser uses the first 18 fields:

- columns 1 to 3: chromosome, 0-based start, and end
- column 4: modification code, often `code,motif,position`
- column 6: strand
- column 10: valid coverage
- column 11: percent modified (non-authoritative; beta is computed from
  count fields)

The parser converts 0-based BED starts to 1-based genomic positions and
computes beta from authoritative count fields (`Nmod / Nvalid_cov`)
rather than the `fraction_modified` percent column. Zero-valid-coverage
rows with missing or undefined `fraction_modified` are dropped by the
coverage filter rather than failing import.

Common symptoms and fixes:

| Symptom | Likely cause | Fix |
|----|----|----|
| `expected at least 18` | The file is not modkit pileup bedMethyl output, or it was truncated. | Re-run `modkit pileup` and pass the BED output, not the BAM. |
| Unknown `mod_code` warning | Column 4 contains codes outside the supported modkit map. | Confirm the file is a methylation pileup and decide whether the code should map to `6mA`, `5mC`, or `4mC`. |
| Missing required fields | One or more rows are partial or malformed after parsing the first 18 bedMethyl columns. | Inspect the reported row, then regenerate or filter the pileup before import. |
| No sites remain | Coverage or `mod_type` filtering removed every row. | Lower `min_coverage`, check `mod_type`, and inspect the input file for nonzero coverage. |
| Duplicate methylation site rows | The file has repeated rows for the same chromosome, position, strand, modification type, and motif. | Aggregate or regenerate the per-sample pileup before importing. |

If the modkit `mod_code` field includes motif context, such as
`a,GATC,1`, commaKit stores `mod_type = "6mA"` and `motif = "GATC"`.
Older or derived files without motif context can still import, but
`motif` will be `NA`.

``` r

first_rows <- read.table("ctrl_1.modkit.bed", nrows = 5, fill = TRUE)
ncol(first_rows)
first_rows[, c(1, 2, 3, 4, 6, 10, 11)]
```

## Sample Metadata Problems

`files` must be a named character vector, and its names must match
`colData$sample_name`. The constructor requires `sample_name` and
`replicate` columns in `colData`; `condition` is optional sample
metadata used by grouping and design workflows when present.

``` r

setdiff(names(files), col_data$sample_name)
setdiff(col_data$sample_name, names(files))
```

If either result is non-empty, fix the names before calling
[`commaData()`](https://carl-stone.github.io/commaKit/reference/commaData.md).
Additional sample metadata columns are allowed and are preserved in
`sampleInfo(object)`.

## Genome Size Problems

Genome information is strongly recommended because it becomes the
`Seqinfo` attached to methylation sites. Use one of these forms:

- a named integer or numeric vector of chromosome sizes
- a FASTA file path
- a named
  [`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
- a whole `BSgenome` object

``` r

# Named vector: simplest and robust for one bacterial chromosome
genome <- c(NC_000913 = 4641652L)

# FASTA path: names come from FASTA sequence headers
genome <- "E_coli_K12.fa"

# Whole BSgenome object: pass the object, not one chromosome extracted with $
genome <- BSgenome.Ecoli.NCBI.20080805::BSgenome.Ecoli.NCBI.20080805
```

Do not pass a single unnamed sequence such as
`BSgenomeObject$NC_000913`. That creates a `DNAString` without a
chromosome name, so commaKit cannot attach the sequence length to
methylation sites. Use the whole `BSgenome` object or a named vector
instead.

Chromosome names must match the methylation files. If the BED file uses
`chr1` but the genome vector is named `NC_000913`, rename one side
before import. commaKit stops with a clear error when methylation data
contain chromosomes that are absent from `genome`, because it cannot
attach valid `Seqinfo` for those sites. Extra chromosomes in the genome
input are dropped with a message when they are not present in the
methylation data.

``` r

bed_chroms <- unique(read.table("ctrl_1.modkit.bed", nrows = 1000)[[1]])
names(genome)
setdiff(bed_chroms, names(genome))
```

## Annotation File Problems

Annotation is optional. Use `annotation = NULL` while debugging
methylation file import, then add annotation once the `commaData` object
can be constructed.

commaKit accepts annotation as a `GRanges` object or as a GFF, GFF3, or
BED file path. Loading annotation files requires the Bioconductor
package `rtracklayer`.

``` r

if (requireNamespace("rtracklayer", quietly = TRUE)) {
  genes <- loadAnnotation("genes.gff3", feature_types = "gene")
  obj <- commaData(files, col_data,
    genome = c(chr1 = 4641652L),
    annotation = genes
  )
}
```

Common symptoms and fixes:

| Symptom | Likely cause | Fix |
|----|----|----|
| `Package 'rtracklayer' is required` | The optional annotation dependency is not installed. | Install with `BiocManager::install("rtracklayer")`, or set `annotation = NULL`. |
| `Annotation file not found` | The path is wrong relative to the R working directory. | Use [`normalizePath()`](https://rdrr.io/r/base/normalizePath.html) or an absolute path to confirm the file. |
| Unsupported extension | The file extension is not `.gff`, `.gff3`, or `.bed` after removing compression suffixes. | Rename or convert the file, or import it yourself as `GRanges`. |
| No requested feature types | `feature_types` does not match the standardized feature type column. | Inspect `table(mcols(loadAnnotation(file))$feature_type)`. |

## Quick Checklist

Before opening an issue or debugging further, collect these facts:

``` r

file.exists(files)
names(files)
col_data$sample_name

first_rows <- read.table(files[[1]], nrows = 5, fill = TRUE)
dim(first_rows)
first_rows[, seq_len(min(18, ncol(first_rows)))]

genome
```

Then try a minimal constructor call without annotation or motif
scanning:

``` r

obj <- commaData(
  files = files,
  colData = col_data,
  genome = c(chr1 = 4641652L),
  annotation = NULL,
  motif = NULL,
  min_coverage = 1L
)
```

If this works, add filters, annotation, motif scanning, and higher
`min_coverage` one at a time. That usually isolates whether the failure
is in file parsing, sample metadata, genome information, or optional
annotation.

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
#> [1] BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39       desc_1.4.3          R6_2.6.1           
#>  [4] bookdown_0.47       fastmap_1.2.0       xfun_0.60          
#>  [7] cachem_1.1.0        knitr_1.51          htmltools_0.5.9    
#> [10] rmarkdown_2.31      lifecycle_1.0.5     cli_3.6.6          
#> [13] sass_0.4.10         pkgdown_2.2.1       textshaping_1.0.5  
#> [16] jquerylib_0.1.4     systemfonts_1.3.2   compiler_4.5.2     
#> [19] tools_4.5.2         ragg_1.5.2          bslib_0.11.0       
#> [22] evaluate_1.0.5      yaml_2.3.12         BiocManager_1.30.27
#> [25] otel_0.2.0          jsonlite_2.0.0      rlang_1.3.0        
#> [28] fs_2.1.0            htmlwidgets_1.6.4
```
