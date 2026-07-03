# Build a KEGG term-to-gene mapping for use with enrichMethylation()

Fetches all KEGG pathway-gene associations for an organism using only
two API calls (`keggLink` and `keggList`), then returns the result as a
pair of data frames that can be passed directly to
[`enrichMethylation()`](https://carl-stone.github.io/commaKit/reference/enrichMethylation.md)
via the `kegg_term2gene` and `kegg_term2name` arguments. Optionally
caches the result to an RDS file so that subsequent calls are fully
offline.

## Usage

``` r
buildKEGGTermGene(organism, file = NULL, strip_prefix = TRUE, id_map = NULL)
```

## Arguments

- organism:

  Character string; KEGG organism code (e.g., `"eco"` for *Escherichia
  coli* K-12, `"hsa"` for *Homo sapiens*). Browse organism codes at
  <https://www.genome.jp/kegg/catalog/org_list.html>.

- file:

  Character string or `NULL`. Path to an RDS cache file.

  - If `file` already exists, the cached object is loaded and returned
    immediately — no network access occurs.

  - If `file` does not exist, KEGG is queried and the result is saved to
    `file`.

  - `NULL` (default) disables caching.

  A warning is issued when a cache file is older than 90 days,
  suggesting the user refresh it to pick up pathway database updates.

- strip_prefix:

  Logical; if `TRUE` (default), the organism prefix is stripped from
  gene IDs (e.g., `"eco:b0001"` becomes `"b0001"`) and the `"path:"`
  prefix is stripped from pathway IDs (e.g., `"path:eco00010"` becomes
  `"eco00010"`). The trailing organism qualifier is also stripped from
  pathway names (e.g., `"... - Escherichia coli K-12"` becomes `"..."`).

- id_map:

  A `data.frame` with columns `symbol` and `kegg_id` as returned by
  [`buildKEGGGeneIDMap()`](https://carl-stone.github.io/commaKit/reference/buildKEGGGeneIDMap.md).
  When provided, the `gene` column of `term2gene` is translated from
  KEGG-internal IDs (e.g. b-numbers) to gene symbols before returning,
  so that genes match the identifiers in your annotation data. KEGG IDs
  with no matching symbol are preserved as-is; no pathway genes are
  dropped. `NULL` (default) leaves identifiers unchanged.

## Value

A named list with two elements:

- `term2gene`:

  A `data.frame` with columns `term` (KEGG pathway ID) and `gene` (gene
  ID matching the identifiers in your annotation data).

- `term2name`:

  A `data.frame` with columns `term` (KEGG pathway ID) and `name`
  (human-readable pathway description).

Pass these to
[`enrichMethylation()`](https://carl-stone.github.io/commaKit/reference/enrichMethylation.md)
as
`kegg_term2gene = result$term2gene, kegg_term2name = result$term2name`.

## Details

This function exists because
[`enrichKEGG`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html)
and [`gseKEGG`](https://rdrr.io/pkg/clusterProfiler/man/gseKEGG.html)
fire one HTTP request per pathway when fetching gene lists, which
quickly exceeds the KEGG API rate limit for organisms with many
pathways. `buildKEGGTermGene` retrieves the same data in two bulk calls
— one for gene-pathway links and one for pathway names.

## See also

[`buildKEGGGeneIDMap`](https://carl-stone.github.io/commaKit/reference/buildKEGGGeneIDMap.md),
[`enrichMethylation`](https://carl-stone.github.io/commaKit/reference/enrichMethylation.md)

## Examples

``` r
# \donttest{
if (requireNamespace("KEGGREST", quietly = TRUE)) {
  # Fetch once and cache to disk (use tempdir() so check directory stays clean)
  kegg_file <- file.path(tempdir(), "eco_kegg.rds")
  kegg <- buildKEGGTermGene("eco", file = kegg_file)

  # Load from cache on subsequent calls (no network)
  kegg <- buildKEGGTermGene("eco", file = kegg_file)

  # Translate b-numbers to symbols with an id_map:
  # id_map <- buildKEGGGeneIDMap("eco", OrgDb = org.EcK12.eg.db)
  # kegg   <- buildKEGGTermGene("eco", file = kegg_file, id_map = id_map)
  # res    <- enrichMethylation(obj,
  #             kegg_term2gene = kegg$term2gene,
  #             kegg_term2name = kegg$term2name)
}
#> Loading KEGG data from cache: /tmp/RtmpRS3aNN/eco_kegg.rds
#> Loading KEGG data from cache: /tmp/RtmpRS3aNN/eco_kegg.rds
# }
```
