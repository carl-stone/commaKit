# Parse a Dorado BAM file with MM/ML base modification tags

Reads a Dorado-aligned BAM file containing MM (base modification) and ML
(modification likelihood) tags, and aggregates per-read modification
probabilities into per-site beta values. This is an internal function
called by
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
when `caller = "dorado"`.

## Usage

``` r
.parseDorado(file, sample_name, mod_type = NULL, min_coverage = 5L)
```

## Arguments

- file:

  Character string. Path to the Dorado-aligned BAM file. An accompanying
  `.bai` index file must exist (required by `Rsamtools`).

- sample_name:

  Character string. Sample name (used in messages only).

- mod_type:

  Character vector or `NULL`. If provided, only sites with a matching
  modification type are returned. `NULL` retains all types.

- min_coverage:

  Integer. Minimum read depth required to retain a site. Default `5L`.

## Value

A `data.frame` with columns:

- `chrom`:

  Chromosome name (character).

- `position`:

  1-based genomic position (integer).

- `strand`:

  Strand, `"+"` or `"-"` (character).

- `mod_type`:

  Modification type: `"6mA"`, `"5mC"`, or `"4mC"` (character).

- `motif`:

  Always `NA` — Dorado BAM MM/ML tags do not encode sequence context
  motifs. Motif information can be added post-hoc with
  [`findMotifSites`](https://carl-stone.github.io/commaKit/reference/findMotifSites.md)
  if needed.

- `beta`:

  Proportion of reads called methylated, range 0–1 (numeric).

- `coverage`:

  Total read depth at this site (integer).

- `mod_counts`:

  Reads called modified using the ML probability threshold.

- `canonical_counts`:

  Reads not called modified using the ML probability threshold.

- `other_mod_counts`:

  Always `NA`; direct Dorado BAM parsing does not decompose non-target
  modified-read classes.

## Details

The function reads the BAM using
[`scanBam`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html) and
processes the MM/ML tags from each aligned read. The MM tag encodes
which bases carry modifications and their positions within the read
sequence (as inter-base offsets); the parallel ML tag provides
modification probabilities (0-255 scaled to 0-1).

A base is called modified when its ML probability is \\\> 0.5\\.
Per-site statistics are aggregated by counting modified reads
(\\n\_{\text{mod}}\\) and total reads at each genomic position
(\\n\_{\text{total}}\\). The beta value is \\\beta = n\_{\text{mod}} /
n\_{\text{total}}\\.

CIGAR operations are used to map read positions to reference
coordinates. Soft-clipped (`S`) and hard-clipped (`H`) bases are
excluded; insertions (`I`) consume read positions without advancing the
reference.

Two kinds of information loss can occur and should be distinguished:

- Read-level skips:

  An entire read is skipped when its CIGAR string is malformed or its
  MM/ML tags are missing or truncated. No modified-base calls from that
  read contribute to any site.

- Call-level drops:

  Individual modified-base calls whose read positions cannot be mapped
  to reference sites (inserted or soft-clipped positions) are dropped
  before site aggregation. The rest of the read is still used; only the
  unmappable call is discarded.

**Recommended workflow:** For most users, it is simpler to first run
`modkit pileup` on the Dorado BAM, then load the resulting BED file with
`caller = "modkit"`. Direct BAM parsing is provided for users who prefer
to avoid the modkit step.
