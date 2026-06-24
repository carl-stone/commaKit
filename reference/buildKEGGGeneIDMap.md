# Build a KEGG gene ID map for symbol translation

Fetches the complete NCBI Gene ID \\\leftrightarrow\\ KEGG gene ID
correspondence for an organism in a single API call
([`KEGGREST::keggConv`](https://rdrr.io/pkg/KEGGREST/man/keggConv.html)),
then joins it against a gene symbol table supplied either as a
Bioconductor `OrgDb` object or as a plain `data.frame`. The result maps
each gene symbol to the organism's KEGG-internal gene identifier (e.g.
b-numbers for *E. coli* K-12), and can be passed to
[`buildKEGGTermGene()`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md)
to translate the `gene` column of `term2gene` from KEGG IDs to symbols
that match your annotation data.

## Usage

``` r
buildKEGGGeneIDMap(
  organism,
  OrgDb = NULL,
  entrez2symbol = NULL,
  keys_col = "SYMBOL",
  id_col = "ENTREZID",
  file = NULL
)
```

## Arguments

- organism:

  Character string; KEGG organism code (e.g., `"eco"` for *Escherichia
  coli* K-12). See <https://www.genome.jp/kegg/catalog/org_list.html>.

- OrgDb:

  A Bioconductor `OrgDb` object for automatic symbol lookup. `NULL`
  (default) if `entrez2symbol` is provided instead.

- entrez2symbol:

  A `data.frame` with columns `entrez_id` and `symbol`. `NULL` (default)
  if `OrgDb` is provided.

- keys_col:

  Character string; the `OrgDb` column containing gene symbols. Default
  `"SYMBOL"`. Ignored when `entrez2symbol` is provided.

- id_col:

  Character string; the `OrgDb` column containing NCBI Gene IDs. Default
  `"ENTREZID"`. Ignored when `entrez2symbol` is provided.

- file:

  Character string or `NULL`. Path to an RDS cache file.

  - If `file` already exists, it is loaded and returned without any API
    calls.

  - If `file` does not exist, the mapping is built and saved.

  - `NULL` (default) disables caching.

  A warning is issued when the cache is older than 90 days.

## Value

A `data.frame` with columns:

- `symbol`:

  Gene symbol matching identifiers in your annotation data (e.g.
  `"lacZ"`, `"rpoD"`).

- `kegg_id`:

  Corresponding KEGG gene identifier after stripping the organism prefix
  (e.g. `"b0344"` for *E. coli* K-12).

Pass this to
[`buildKEGGTermGene()`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md)
as `id_map = ...` to translate the `gene` column of `term2gene` from
KEGG IDs to symbols.

## Why this is needed

[`buildKEGGTermGene()`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md)
returns pathway-gene associations where genes are identified by
KEGG-internal IDs (b-numbers for *E. coli*, Entrez IDs for human, etc.).
Annotation data produced by
[`loadAnnotation()`](https://carl-stone.github.io/commaKit/reference/loadAnnotation.md)
uses gene symbols from the GFF3 `Name=` attribute. Without translation,
no genes will match during enrichment.

## Input modes

Exactly one of `OrgDb` or `entrez2symbol` must be supplied:

- `OrgDb`:

  A Bioconductor `OrgDb` annotation object (e.g., `org.EcK12.eg.db`).
  NCBI Gene IDs (`id_col`) and gene symbols (`keys_col`) are extracted
  via
  [`AnnotationDbi::select()`](https://rdrr.io/pkg/AnnotationDbi/man/AnnotationDb-class.html).
  Requires internet access only for the single `keggConv` call.

- `entrez2symbol`:

  A two-column `data.frame` with columns `entrez_id` (character NCBI
  Gene IDs) and `symbol` (gene symbols). Use this when no `OrgDb`
  package is available, e.g. when you have downloaded a mapping from
  NCBI's `gene_info` file.

## See also

[`buildKEGGTermGene`](https://carl-stone.github.io/commaKit/reference/buildKEGGTermGene.md),
[`enrichMethylation`](https://carl-stone.github.io/commaKit/reference/enrichMethylation.md)

## Examples

``` r
# \donttest{
if (requireNamespace("KEGGREST", quietly = TRUE) &&
    requireNamespace("org.EcK12.eg.db", quietly = TRUE)) {
  id_file  <- file.path(tempdir(), "eco_id_map.rds")
  id_map <- buildKEGGGeneIDMap("eco",
                               OrgDb = org.EcK12.eg.db::org.EcK12.eg.db,
                               file  = id_file)

  # Apply when building pathway map:
  kegg_file <- file.path(tempdir(), "eco_kegg.rds")
  kegg <- buildKEGGTermGene("eco", file = kegg_file, id_map = id_map)
  # kegg$term2gene$gene now contains symbols instead of b-numbers
}
#> Fetching KEGG gene ID map for organism 'eco' ...
#> 'select()' returned 1:1 mapping between keys and columns
#> KEGG ID map cached to: /tmp/RtmpJlNtDA/eco_id_map.rds
#> Done. 4634 gene symbols mapped to KEGG IDs.
#> Fetching KEGG pathway data for organism 'eco' ...
#> KEGG data cached to: /tmp/RtmpJlNtDA/eco_kegg.rds
#> Done. 4984 gene-pathway associations across 138 pathways.

# Manual table alternative:
if (requireNamespace("KEGGREST", quietly = TRUE)) {
  ent2sym <- data.frame(
    entrez_id = c("945076", "945803"),
    symbol    = c("lacZ",   "lacY"),
    stringsAsFactors = FALSE
  )
  id_map <- buildKEGGGeneIDMap("eco", entrez2symbol = ent2sym)
}
#> Fetching KEGG gene ID map for organism 'eco' ...
#> Done. 2 gene symbols mapped to KEGG IDs.
# }
```
