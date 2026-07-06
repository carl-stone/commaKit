# Create a commaData object from methylation calling output files

Constructor for the
[`commaData-class`](https://carl-stone.github.io/commaKit/reference/commaData-class.md)
S4 class. Parses one or more methylation calling output files (modkit,
Megalodon, or Dorado), merges them into a sites × samples matrix
representation, and optionally loads genomic annotation and motif site
positions.

## Usage

``` r
commaData(
  files,
  colData,
  genome = NULL,
  annotation = NULL,
  mod_type = NULL,
  motif = NULL,
  expected_mod_contexts = NULL,
  min_coverage = 5L,
  caller = "modkit"
)
```

## Arguments

- files:

  Named character vector mapping sample names to file paths. Names must
  match `colData$sample_name`. Example:
  `c(ctrl_1 = "/path/to/ctrl_1.bed", treat_1 = "/path/to/treat_1.bed")`.

- colData:

  A `data.frame` with one row per sample. Must contain columns
  `sample_name` and `replicate`. A `condition` column is optional and is
  used by design-aware functions such as
  [`diffMethyl`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
  when requested in their formula. Additional columns (e.g.,
  `file_path`, `batch`) are preserved.

- genome:

  Genome size information: a named integer vector of chromosome sizes
  (e.g., `c(NC_000913 = 4641652L)`), a path to a FASTA file, a
  `DNAStringSet` (Biostrings), or a `BSgenome` object. For
  single-chromosome genomes pass the `BSgenome` object directly or a
  named integer vector — do not index into the BSgenome with `$` (e.g.,
  `BSgenome.Ecoli.NCBI.20080805$NC_000913`) as that yields a `DNAString`
  which has no chromosome name and cannot be used. Set to `NULL` to omit
  genome information (not recommended). When a multi-sequence source is
  provided, genome info is automatically restricted to chromosomes
  present in the data. Chromosomes are assumed circular by default when
  genome information is attached, matching the package's bacterial
  genome default; edit `seqinfo(rowRanges(object))` after construction
  if a chromosome should be treated as linear.

- annotation:

  Optional. Path to a GFF3 or BED annotation file, or a pre-loaded
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object. If `NULL`, the annotation slot is left empty.

- mod_type:

  Optional character vector specifying which modification types to
  retain (e.g., `"6mA"` or `c("6mA", "5mC")`). If `NULL`, all
  modification types detected in the files are kept.

- motif:

  Optional character string. A DNA sequence motif (e.g., `"GATC"`) to
  locate in the genome using
  [`findMotifSites`](https://carl-stone.github.io/commaKit/reference/findMotifSites.md).
  The results are stored in the `motifSites` slot as a genome-wide
  `GRanges` of all motif instances. Requires `genome` to be a FASTA
  path, `DNAStringSet`, or `BSgenome` object (not a named integer
  vector). If `NULL`, the `motifSites` slot is left empty. *Note:* this
  argument is distinct from `rowData(object)$motif`, which stores the
  per-site sequence context extracted automatically from the modkit
  `mod_code` field (e.g., `"a,GATC,1"` → `motif = "GATC"`) and is `NA`
  for Dorado and Megalodon callers.

- expected_mod_contexts:

  Named list or `NULL`. If provided, specifies which modification type /
  sequence motif combinations to retain. Names must be modification type
  strings (e.g., `"6mA"`, `"5mC"`). Values are character vectors of
  motif strings (e.g., `c("GATC", "ACCACC")`). Sites whose `mod_context`
  (`paste(mod_type, motif, sep = "_")`) does not match any name–value
  pair are dropped before the object is assembled. A message is emitted
  reporting the number of sites dropped per modification type. Use
  `NULL` (default) to retain all sites. Example:
  `list("6mA" = "GATC", "5mC" = c("CCWGG", "CCGG"))`. *Note:* for
  Dorado/Megalodon callers where `motif` is `NA`, the `mod_context`
  falls back to just `mod_type` (e.g., `"6mA"`), so those sites are only
  retained if you include `NA` in the motif vector for that type (e.g.,
  `list("6mA" = NA)`).

- min_coverage:

  Integer. Minimum read depth to include a site. Sites present in a
  sample with coverage below this threshold have their beta value set to
  `NA`. Sites absent from a sample entirely are also `NA`. Default `5`.

- caller:

  Character string specifying the methylation caller that produced the
  input files. One of `"modkit"` (default), `"megalodon"`, or
  `"dorado"`.

## Value

A valid `commaData` object.

## Details

The constructor uses a parse-then-merge strategy:

1.  Each file is parsed independently using the appropriate parser.

2.  Sites are identified by their genomic position (chromosome,
    position, strand) plus modification type and motif context.

3.  The union of all sites across all samples is taken, using
    `findOverlaps()` for alignment.

4.  Beta values and coverage are arranged into sites \\\times\\ samples
    matrices, with `NA` for samples that do not cover a given site.

5.  Observed modified, canonical, and non-target modified read counts
    are preserved as `mod_counts`, `canonical_counts`, and
    `other_mod_counts` assays when reported by the caller;
    probability-only callers store `NA` in those assays.

6.  Assay-layer provenance and default roles are recorded in
    `metadata(object)$assay_provenance` and
    `metadata(object)$assay_defaults`.

7.  Sites where coverage is below `min_coverage` in a sample have their
    beta value set to `NA` (but coverage is preserved).

## See also

[`commaData-class`](https://carl-stone.github.io/commaKit/reference/commaData-class.md),
[`methylation`](https://carl-stone.github.io/commaKit/reference/methylation.md),
[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`modCounts`](https://carl-stone.github.io/commaKit/reference/modCounts.md),
[`canonicalCounts`](https://carl-stone.github.io/commaKit/reference/canonicalCounts.md),
[`assayLayers`](https://carl-stone.github.io/commaKit/reference/assayLayers.md),
[`assayProvenance`](https://carl-stone.github.io/commaKit/reference/assayProvenance.md),
[`sampleInfo`](https://carl-stone.github.io/commaKit/reference/sampleInfo.md),
[`siteInfo`](https://carl-stone.github.io/commaKit/reference/siteInfo.md),
[`modTypes`](https://carl-stone.github.io/commaKit/reference/modTypes.md),
[`loadAnnotation`](https://carl-stone.github.io/commaKit/reference/loadAnnotation.md),
[`findMotifSites`](https://carl-stone.github.io/commaKit/reference/findMotifSites.md)

## Examples

``` r
# Construct a commaData object from built-in example data
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
#> caller: modkit
#> min_coverage: 5

if (FALSE) { # \dontrun{
# Load two modkit BED files (requires user-provided files)
cd <- commaData(
  files = c(
    ctrl_1  = "ctrl_1_modkit.bed",
    treat_1 = "treat_1_modkit.bed"
  ),
  colData = data.frame(
    sample_name = c("ctrl_1", "treat_1"),
    condition   = c("control", "treatment"),
    replicate   = c(1L, 1L)
  ),
  genome = c(chr1 = 4641652L),
  annotation = "MG1655.gff3",
  caller = "modkit"
)
cd
} # }
```
