# Parse a modkit pileup BED file into a tidy per-site data frame

Reads a single-sample modkit `pileup` output file and returns a tidy
data frame of per-site methylation values. This is an internal function
called by
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md).

## Usage

``` r
.parseModkit(file, sample_name, mod_type = NULL, min_coverage = 5L)
```

## Arguments

- file:

  Character string. Path to the modkit pileup BED file.

- sample_name:

  Character string. Name for this sample (used in messages only; not
  added to the returned data frame).

- mod_type:

  Character vector or `NULL`. If provided, only sites with a matching
  modification type (e.g., `"6mA"`, `"5mC"`, `"4mC"`) are returned.
  `NULL` retains all modification types.

- min_coverage:

  Integer. Minimum read depth required to retain a site. Sites with
  `coverage < min_coverage` are dropped. Default `5`.

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

  Sequence context motif extracted from the modkit `mod_code` field
  (e.g., `"GATC"`, `"CCWGG"`). `NA` if the `mod_code` field does not
  contain motif information (older modkit formats).

- `beta`:

  Proportion of reads called methylated, computed as `Nmod / Nvalid_cov`
  from authoritative count fields. Range 0-1 for positive-coverage rows;
  0 for zero-coverage rows.

- `coverage`:

  Total read depth at this site (integer).

- `mod_counts`:

  Observed reads called as this modification.

- `canonical_counts`:

  Observed reads called canonical/unmodified.

- `other_mod_counts`:

  Observed reads called as a non-target modification. Together,
  `mod_counts + canonical_counts + other_mod_counts` is the modkit
  `Nvalid_cov` denominator used for `coverage`.
