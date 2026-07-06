# Genome browser-style methylation track plot

Plots methylation beta values for individual sites along a chromosome
region in a genome browser-style layout, with one panel per sample.
Optionally overlays genomic feature annotations as colored rectangles in
a separate track below the methylation data.

## Usage

``` r
plot_genome_track(
  object,
  chromosome,
  start = NULL,
  end = NULL,
  mod_type = NULL,
  motif = NULL,
  mod_context = NULL,
  annotation = NULL
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  object.

- chromosome:

  Character string. The chromosome (sequence name) to plot. Must be
  present in `names(genomeSizes(object))`.

- start:

  Integer or `NULL`. Start position of the region to display (1-based,
  inclusive). If `NULL`, the plot begins at position 1.

- end:

  Integer or `NULL`. End position of the region to display (1-based,
  inclusive). If `NULL`, the plot extends to the end of the chromosome.

- mod_type:

  Character vector of modification types to display (e.g., `"6mA"`,
  `c("6mA", "5mC")`). If `NULL` (default), all modification types are
  shown, colored differently.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are displayed (e.g., `"GATC"`). If `NULL`
  (default), all motifs are included.

- mod_context:

  Character vector or `NULL`. If provided, only sites whose
  `mod_context` rowData column matches one of the supplied values are
  included. If `NULL` (default), all modification contexts are included.

- annotation:

  A `GRanges` object of genomic features to display in the annotation
  track, `NULL` (default, uses `annotation(object)` if available), or
  `FALSE` to suppress the annotation track entirely.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The methylation track shows individual sites as points (x =
genomic position, y = beta value 0–1) colored by modification type,
faceted by sample. If annotation features are present on the selected
chromosome (and `annotation` is not `FALSE`), they are displayed as
colored rectangles in a separate annotation panel below.

## See also

[`annotateSites`](https://carl-stone.github.io/commaKit/reference/annotateSites.md),
[`plot_metagene`](https://carl-stone.github.io/commaKit/reference/plot_metagene.md)

## Examples

``` r
data(comma_example_data)
plot_genome_track(comma_example_data, chromosome = "chr_sim")


# Restrict to a region
plot_genome_track(comma_example_data,
  chromosome = "chr_sim",
  start = 1000, end = 50000
)


# One modification type, no annotation
plot_genome_track(comma_example_data,
  chromosome = "chr_sim",
  mod_type = "6mA", annotation = FALSE
)

```
