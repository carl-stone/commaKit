# Validate and coerce genome input to a named integer vector

Validate and coerce genome input to a named integer vector

## Usage

``` r
.validateGenomeInfo(genome)
```

## Arguments

- genome:

  A BSgenome object, path to a FASTA file, named `DNAStringSet`, or
  named integer vector of chromosome sizes (e.g., `c(chr1 = 1000000L)`).

## Value

Named integer vector of chromosome sizes, or `NULL` if `genome` is
`NULL`.
