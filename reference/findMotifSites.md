# Find all instances of a sequence motif in a genome

Searches a genome for all occurrences of a DNA sequence motif on both
strands and returns their positions as a
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object. This is the bacterial equivalent of CpG island finding — for
example, locating all GATC (Dam methyltransferase) sites in an *E. coli*
genome.

## Usage

``` r
findMotifSites(genome, motif, ...)
```

## Arguments

- genome:

  A `DNAStringSet` object (Biostrings), a `BSgenome` object, or a
  character string giving the path to a FASTA file containing the genome
  sequence.

- motif:

  A character string specifying the DNA sequence motif to search for
  (e.g., `"GATC"` or `"CCWGG"`). Standard IUPAC ambiguity codes are
  supported (e.g., W = A or T). The motif is matched exactly on the
  forward strand; its reverse complement is searched on the minus
  strand.

- ...:

  Additional arguments (currently unused; reserved for future use).

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object with one range per motif instance. Ranges are 1-based and have
width equal to `nchar(motif)`. The `mcols` column `motif` records the
search motif string.

## Details

This function requires the `Biostrings` package (Bioconductor). If
`genome` is a `BSgenome` object, `BSgenome` is also required. Both can
be installed with:

      BiocManager::install(c("Biostrings", "BSgenome"))

Palindromic motifs (e.g., GATC) will have both forward and reverse
strand hits at each location. Non-palindromic motifs will have distinct
forward and reverse hits.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find all GATC sites in a genome from a FASTA file
gatc_sites <- findMotifSites(genome = "MG1655.fa", motif = "GATC")
} # }

if (FALSE) { # \dontrun{
# With a BSgenome object (requires BSgenome.Ecoli.NCBI.20080805)
# gatc_sites <- findMotifSites(genome = BSgenome.Ecoli.NCBI.20080805, motif = "GATC")
} # }
```
