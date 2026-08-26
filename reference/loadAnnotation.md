# Load genomic feature annotations from a GFF3 or BED file

Reads a GFF3 or BED annotation file and returns a
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object with standardized metadata columns. The result can be passed
directly to the `annotation` argument of
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md).

## Usage

``` r
loadAnnotation(file, feature_types = NULL, ...)
```

## Arguments

- file:

  Character string. Path to a GFF3 (`.gff`, `.gff3`, `.gff.gz`,
  `.gff3.gz`) or BED (`.bed`) file.

- feature_types:

  Character vector or `NULL`. If provided, only features with a matching
  `type` (GFF3) or `name` (BED) are retained. Common GFF3 types include
  `"gene"`, `"CDS"`, `"rRNA"`, `"tRNA"`. `NULL` retains all features.

- ...:

  Additional arguments passed to `rtracklayer::import()`.

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object. The `mcols` always include:

- `feature_type`:

  Character. The feature type (from GFF3 `type` column, or `"region"`
  for BED).

- `name`:

  Character. Feature name or identifier (from GFF3 `Name` or `ID`
  attribute, or BED `name` column).

GFF3 files may also produce:

- `feature_subtype`:

  Character. The `feature_type` GFF3 *attribute* (column 9), if present.
  For EcoCyc annotations this encodes sigma factor identity
  (`"Sigma70"`, `"Sigma24"`, etc.) for
  `transcription_factor_binding_site` features.

- `transcription_unit`:

  Character. The regulated transcription unit(s) for protein and RNA
  binding sites (EcoCyc format).

All other metadata columns from the source file are preserved.

## Details

This function requires the `rtracklayer` package (Bioconductor):

      BiocManager::install("rtracklayer")

Both NCBI-style GFF3 (with `gene_biotype`, `product` attributes) and
Ensembl-style GFF3 are supported.

## Examples

``` r
# Load the bundled example GFF3 annotation
gff_file <- system.file("extdata", "example.gff3", package = "commaKit")
if (requireNamespace("rtracklayer", quietly = TRUE)) {
  ann <- loadAnnotation(gff_file)
  ann <- loadAnnotation(gff_file, feature_types = "gene")
}

if (FALSE) { # \dontrun{
# Load only genes and CDS from your own file
ann <- loadAnnotation("my_genome.gff3", feature_types = c("gene", "CDS"))
} # }
```
