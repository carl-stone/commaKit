#' @importFrom methods is
NULL

# Internal genome utility functions used by the commaData constructor and
# sliding window analyses.

#' Validate and coerce genome input to a named integer vector of chromosome sizes
#'
#' @param genome A BSgenome object, path to a FASTA file, or named integer
#'   vector of chromosome sizes (e.g., \code{c(chr1 = 1000000L)}).
#'
#' @return Named integer vector of chromosome sizes, or \code{NULL} if
#'   \code{genome} is \code{NULL}.
#'
#' @keywords internal
.validateGenomeInfo <- function(genome) {
    if (is.null(genome)) {
        return(NULL)
    }

    if (is.integer(genome) && !is.null(names(genome))) {
        return(genome)
    }

    if (is.numeric(genome) && !is.null(names(genome))) {
        result <- as.integer(genome)
        names(result) <- names(genome)
        return(result)
    }

    if (is.character(genome) && length(genome) == 1) {
        if (!file.exists(genome)) {
            stop("genome: file not found: ", genome)
        }
        if (!requireNamespace("Biostrings", quietly = TRUE)) {
            stop(
                "Package 'Biostrings' is required to read FASTA genome files. ",
                "Install it with: BiocManager::install('Biostrings')"
            )
        }
        seqs <- Biostrings::readDNAStringSet(genome)
        sizes <- as.integer(Biostrings::width(seqs))
        names(sizes) <- names(seqs)
        return(sizes)
    }

    # BSgenome object
    if (is(genome, "BSgenome")) {
        if (!requireNamespace("BSgenome", quietly = TRUE)) {
            stop(
                "Package 'BSgenome' is required to use BSgenome objects. ",
                "Install it with: BiocManager::install('BSgenome')"
            )
        }
        chr_names <- BSgenome::seqnames(genome)
        chr_sizes <- as.integer(GenomeInfoDb::seqlengths(genome))
        names(chr_sizes) <- chr_names
        return(chr_sizes)
    }

    stop(
        "genome must be a named integer vector of chromosome sizes, ",
        "a path to a FASTA file, or a BSgenome object. Got: ",
        class(genome)
    )
}

#' Vectorized circular genome index
#'
#' Converts a position that may wrap past the end (or before the start) of a
#' circular genome into a valid 1-based position.
#'
#' @param position Integer vector of genomic positions (1-based).
#' @param genome_size Integer. Total size of the circular genome in base pairs.
#'
#' @return Integer vector of positions in the range \code{[1, genome_size]}.
#'
#' @keywords internal
.circularIndex <- function(position, genome_size) {
    stopifnot(is.numeric(genome_size), length(genome_size) == 1, genome_size > 0)
    ((as.integer(position) - 1L) %% as.integer(genome_size)) + 1L
}

#' Build a Seqinfo object from a named integer vector of chromosome sizes
#'
#' @param genome_info Named integer vector of chromosome sizes, as returned by
#'   \code{.validateGenomeInfo()}.
#' @param genome_name Optional character string for the genome name.
#'
#' @return A \code{GenomeInfoDb::Seqinfo} object, or \code{NULL} if
#'   \code{genome_info} is \code{NULL}.
#'
#' @importFrom GenomicRanges seqinfo
#' @keywords internal
.makeSeqinfo <- function(genome_info, genome_name = NA_character_) {
    if (is.null(genome_info)) {
        return(NULL)
    }
    GenomeInfoDb::Seqinfo(
        seqnames  = names(genome_info),
        seqlengths = genome_info,
        isCircular = rep(NA, length(genome_info)),
        genome     = genome_name
    )
}
