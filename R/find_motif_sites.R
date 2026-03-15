#' @importFrom GenomicRanges GRanges mcols<-
#' @importFrom IRanges IRanges
NULL

#' Find all instances of a sequence motif in a genome
#'
#' Searches a genome for all occurrences of a DNA sequence motif on both
#' strands and returns their positions as a \code{\link[GenomicRanges]{GRanges}}
#' object. This is the bacterial equivalent of CpG island finding — for
#' example, locating all GATC (Dam methyltransferase) sites in an \emph{E. coli}
#' genome.
#'
#' @param genome A \code{BSgenome} object, or a character string giving the
#'   path to a FASTA file containing the genome sequence.
#' @param motif A character string specifying the DNA sequence motif to search
#'   for (e.g., \code{"GATC"} or \code{"CCWGG"}). Standard IUPAC ambiguity
#'   codes are supported (e.g., W = A or T). The motif is matched exactly on
#'   the forward strand; its reverse complement is searched on the minus strand.
#' @param ... Additional arguments (currently unused; reserved for future use).
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object with one range per
#'   motif instance. Ranges are 1-based and have width equal to
#'   \code{nchar(motif)}. The \code{mcols} column \code{motif} records the
#'   search motif string.
#'
#' @details
#' This function requires the \code{Biostrings} package (Bioconductor). If
#' \code{genome} is a \code{BSgenome} object, \code{BSgenome} is also required.
#' Both can be installed with:
#' \preformatted{
#'   BiocManager::install(c("Biostrings", "BSgenome"))
#' }
#' Palindromic motifs (e.g., GATC) will have both forward and reverse strand
#' hits at each location. Non-palindromic motifs will have distinct forward and
#' reverse hits.
#'
#' @examples
#' \donttest{
#' # Find all GATC sites in a genome from a FASTA file
#' gatc_sites <- findMotifSites(genome = "MG1655.fa", motif = "GATC")
#'
#' # With a BSgenome object
#' library(BSgenome.Ecoli.NCBI.20080805)
#' gatc_sites <- findMotifSites(genome = BSgenome.Ecoli.NCBI.20080805, motif = "GATC")
#' }
#'
#' @export
findMotifSites <- function(genome, motif, ...) {
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
        stop(
            "Package 'Biostrings' is required for findMotifSites(). ",
            "Install it with: BiocManager::install('Biostrings')"
        )
    }

    if (missing(genome) || is.null(genome)) {
        stop("genome must be provided (BSgenome object or path to FASTA file)")
    }
    if (!is.character(motif) || length(motif) != 1 || nchar(motif) == 0) {
        stop("motif must be a non-empty character string (e.g., 'GATC')")
    }

    # ── Load genome sequences ───────────────────────────────────────────────
    seqs <- .loadGenomeSequences(genome)

    # ── Search each chromosome ──────────────────────────────────────────────
    motif_dna    <- Biostrings::DNAString(motif)
    motif_rc     <- Biostrings::reverseComplement(motif_dna)
    is_palindrome <- as.character(motif_dna) == as.character(motif_rc)

    all_ranges <- lapply(seq_along(seqs), function(i) {
        chr_name <- names(seqs)[i]
        seq_i    <- seqs[[i]]

        # Forward strand hits
        fwd_hits <- Biostrings::matchPattern(motif_dna, seq_i, fixed = FALSE)
        fwd_starts <- BiocGenerics::start(fwd_hits)
        fwd_gr   <- GenomicRanges::GRanges(
            seqnames = rep(chr_name, length(fwd_starts)),
            ranges   = IRanges::IRanges(start = fwd_starts,
                                        end   = BiocGenerics::end(fwd_hits)),
            strand   = rep("+", length(fwd_starts))
        )

        # Reverse strand hits (skip if palindrome to avoid duplicates at same pos)
        if (!is_palindrome) {
            rev_hits <- Biostrings::matchPattern(motif_rc, seq_i, fixed = FALSE)
            rev_starts <- BiocGenerics::start(rev_hits)
            rev_gr   <- GenomicRanges::GRanges(
                seqnames = rep(chr_name, length(rev_starts)),
                ranges   = IRanges::IRanges(start = rev_starts,
                                            end   = BiocGenerics::end(rev_hits)),
                strand   = rep("-", length(rev_starts))
            )
            c(fwd_gr, rev_gr)
        } else {
            # Palindrome: assign both strands to the same positions
            if (length(fwd_gr) > 0L) {
                rev_gr <- fwd_gr
                GenomicRanges::strand(rev_gr) <- "-"
                c(fwd_gr, rev_gr)
            } else {
                fwd_gr
            }
        }
    })

    result <- do.call(c, all_ranges)
    GenomicRanges::mcols(result)$motif <- rep(motif, length(result))
    GenomicRanges::sort(result)
}

#' Load genome sequences from a FASTA path or BSgenome object
#' @return A \code{DNAStringSet} with one element per chromosome.
#' @keywords internal
.loadGenomeSequences <- function(genome) {
    if (is.character(genome)) {
        if (!file.exists(genome)) {
            stop("FASTA genome file not found: ", genome)
        }
        Biostrings::readDNAStringSet(genome)
    } else if (is(genome, "BSgenome")) {
        if (!requireNamespace("BSgenome", quietly = TRUE)) {
            stop("Package 'BSgenome' is required to use BSgenome objects.")
        }
        chr_names <- BSgenome::seqnames(genome)
        seqs      <- BSgenome::getSeq(genome, chr_names)
        names(seqs) <- chr_names
        seqs
    } else {
        stop(
            "genome must be a path to a FASTA file or a BSgenome object. Got: ",
            class(genome)
        )
    }
}
