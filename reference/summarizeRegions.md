# Summarize methylation counts over genomic regions

Aggregates site-level count evidence over user-supplied genomic regions.
This is a descriptive regional summary, not a DMR caller, smoothing
method, or region-level statistical test.

## Usage

``` r
summarizeRegions(
  object,
  regions,
  min_sites = 1L,
  mod_type = NULL,
  motif = NULL,
  mod_context = NULL
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object with count evidence assays.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object defining regions to summarize.

- min_sites:

  Integer. Minimum number of overlapping sites with usable count
  evidence and positive valid coverage required before
  `region_methylation` is reported. Regions with fewer usable sites
  return `NA` methylation for that sample. Default: `1L`.

- mod_type, motif, mod_context:

  Optional site filters. If provided, only sites matching these
  modification annotations are summarized. `mod_context` values use
  commaKit's `"<mod_type>_<motif>"` format, e.g. `"6mA_GATC"`. Requested
  values are validated and misspellings are errors.

## Value

A tidy `data.frame` with one row per region × sample containing region
coordinates/metadata, `sample_name`, `n_sites`, summed count evidence,
and `region_methylation` computed as
`sum(mod_counts) / sum(valid_coverage)`. Regions with no overlaps are
retained with zero counts and `NA` methylation; empty `regions` input
returns a typed 0-row `data.frame`.

## Details

The primary regional methylation statistic is count-based:
\$\$region\\methylation = sum(mod\\counts) / sum(valid\\coverage)\$\$
where `valid_coverage` is the `coverage` assay for sites with
non-missing modified counts and positive coverage. For modkit-style
inputs this corresponds to
`Nvalid_cov = Nmod + Nother_mod + Ncanonical`.

`summarizeRegions()` intentionally does not average beta values by
default. Beta-only/reconstructed-count import paths should record
provenance and message when counts are reconstructed; regional summaries
consume the count assays present in the `commaData` object.

## Examples

``` r
data(comma_example_data)
regions <- GenomicRanges::GRanges(
  seqnames = "chr_sim",
  ranges = IRanges::IRanges(start = 1, end = 5000)
)
summarizeRegions(comma_example_data, regions, mod_type = "6mA")
#>   region_id seqnames start  end width strand sample_name n_sites
#> 1  region_1  chr_sim     1 5000  5000      *      ctrl_1      16
#> 2  region_1  chr_sim     1 5000  5000      *      ctrl_2      16
#> 3  region_1  chr_sim     1 5000  5000      *      ctrl_3      16
#> 4  region_1  chr_sim     1 5000  5000      *     treat_1      16
#> 5  region_1  chr_sim     1 5000  5000      *     treat_2      16
#> 6  region_1  chr_sim     1 5000  5000      *     treat_3      16
#>   total_mod_counts total_valid_coverage region_methylation
#> 1             1195                 1329          0.8991723
#> 2             1150                 1284          0.8956386
#> 3              926                 1026          0.9025341
#> 4             1155                 1473          0.7841141
#> 5              924                 1172          0.7883959
#> 6              875                 1220          0.7172131
#>   total_canonical_counts
#> 1                    134
#> 2                    134
#> 3                    100
#> 4                    318
#> 5                    248
#> 6                    345
```
