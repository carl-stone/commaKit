# Annotate methylation sites relative to genomic features

Assigns genomic feature annotations to methylation sites stored in a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object using
[`findOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).

## Usage

``` r
annotateSites(
  object,
  features = NULL,
  feature_col = "feature_type",
  name_col = "name",
  window = 50L,
  keep = c("all", "overlap", "proximity", "metagene"),
  metadata_cols = NULL
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object.

- features:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of genomic features to annotate against. If `NULL` (default), the
  annotation stored in `object` via
  [`annotation`](https://rdrr.io/pkg/BiocGenerics/man/annotation.html)`(object)`
  is used. Must have mcols columns named by `feature_col` and
  `name_col`.

- feature_col:

  Character string. Name of the `mcols` column in `features` that
  contains the feature type. Default: `"feature_type"`.

- name_col:

  Character string. Name of the `mcols` column in `features` that
  contains the feature name. Default: `"name"`.

- window:

  Integer. Search window in base pairs. All features within this
  distance of each site are returned. Default: `50L`.

- keep:

  Character string controlling which output columns are retained:

  `"all"`

  : (default) Return all four columns for all associations — including
    sites within the window but not inside a feature.

  `"overlap"`

  : Subset each site's associations to features where
    `rel_position == 0` (inside the feature). Drop the `rel_position`
    and `frac_position` columns.

  `"proximity"`

  : Keep all associations. Drop the `frac_position` column.

  `"metagene"`

  : Subset to `rel_position == 0`. Drop the `rel_position` column.
    Retain `frac_position`.

- metadata_cols:

  Character vector or `NULL`. Names of additional `mcols` columns in
  `features` to pass through as parallel list columns in `rowData`. Each
  column `X` is stored as `X_values` (`CharacterList`). Default: `NULL`.

## Value

A
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object identical to `object` except that `rowData` has been extended
with annotation list columns. With `keep = "all"` (default):

- `feature_types`:

  CharacterList. Feature type for each association per site.

- `feature_names`:

  CharacterList. Feature name for each association per site.

- `rel_position`:

  IntegerList. Signed relative position (bp): 0 inside, negative
  upstream, positive downstream.

- `frac_position`:

  NumericList. Fractional position in \\\[0, 1\]\\ inside features; `NA`
  outside.

Intergenic sites (no features within `window`) receive length-0 list
elements in all columns.

## Details

The function searches within `window` bp of each site and returns every
feature found. Four parallel list columns are added to `rowData`:
feature types, feature names, signed relative position, and fractional
position within the feature.

Signed `rel_position` convention (strand-aware):

- **0** — site is *inside* the feature.

- **negative** — site is upstream (before the feature in the direction
  of transcription).

- **positive** — site is downstream (after the feature).

`frac_position` is in \\\[0, 1\]\\ when the site is inside a feature
(TSS = 0, TTS = 1, strand-aware) and `NA` when the site is outside.

All list columns are strictly parallel: element \\j\\ in
`feature_types[[i]]` corresponds to element \\j\\ in
`feature_names[[i]]`, `rel_position[[i]]`, and `frac_position[[i]]`.
This invariant is preserved by any metadata columns added via
`metadata_cols`.

## Examples

``` r
data(comma_example_data)
# Default: unified output with all four columns
annotated <- annotateSites(comma_example_data)
si <- siteInfo(annotated)
# Sites inside at least one feature (rel_position == 0):
sum(sapply(as.list(si$rel_position), function(x) any(x == 0L)))
#> [1] 17
# Fractional position of first inside site:
inside_idx <- which(lengths(si$frac_position) > 0)[1]
si$frac_position[[inside_idx]]
#> [1] 0.8857715

# Backward-compatible overlap output (no rel_position / frac_position):
ann_ov <- annotateSites(comma_example_data, keep = "overlap")
siteInfo(ann_ov)$feature_names[[1]]
#> [1] "geneA"
```
