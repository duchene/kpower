# Alignment I/O: read, write, subset, and concatenate alignments

#' Read a FASTA alignment
#'
#' @param path Path to a FASTA file.
#' @return Named character vector (names = taxa, values = sequences).
read_fasta <- function(path) {
  lines <- readLines(path)
  headers <- which(startsWith(lines, ">"))
  if (length(headers) == 0) stop("No FASTA headers found in: ", path)

  taxa <- sub("^>\\s*", "", lines[headers])
  taxa <- sub("\\s+.*", "", taxa)  # first word only

  seqs <- character(length(headers))
  for (i in seq_along(headers)) {
    start <- headers[i] + 1
    end   <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    seqs[i] <- paste(lines[start:end][nzchar(trimws(lines[start:end]))],
                     collapse = "")
  }
  names(seqs) <- taxa
  seqs
}

#' Read a sequential PHYLIP alignment
#'
#' @param path Path to a PHYLIP file.
#' @return Named character vector (names = taxa, values = sequences).
read_phylip <- function(path) {
  lines <- readLines(path)
  nonempty <- lines[nzchar(trimws(lines))]
  header <- trimws(nonempty[1])
  parts  <- strsplit(header, "\\s+")[[1]]
  ntaxa  <- as.integer(parts[1])

  data_lines <- nonempty[-1]
  if (length(data_lines) < ntaxa) {
    stop("Expected ", ntaxa, " sequence lines in PHYLIP file: ", path)
  }

  taxa <- character(ntaxa)
  seqs <- character(ntaxa)
  for (i in seq_len(ntaxa)) {
    parts   <- strsplit(trimws(data_lines[i]), "\\s+", perl = TRUE)[[1]]
    taxa[i] <- parts[1]
    seqs[i] <- paste(parts[-1], collapse = "")
  }
  names(seqs) <- taxa
  seqs
}

#' Read a FASTA or PHYLIP alignment
#'
#' Detects format from the first non-blank line: `>` = FASTA, otherwise PHYLIP.
#'
#' @param path Path to the alignment file.
#' @return Named character vector (names = taxa, values = sequences).
read_alignment <- function(path) {
  lines <- readLines(path, n = 5)
  first <- trimws(lines[nzchar(trimws(lines))][1])
  if (startsWith(first, ">")) read_fasta(path) else read_phylip(path)
}

#' Write a sequential PHYLIP alignment
#'
#' @param seqs Named character vector (names = taxa, values = sequences).
#' @param path Output file path.
write_phylip <- function(seqs, path) {
  ntaxa  <- length(seqs)
  nsites <- nchar(seqs[1])
  out    <- character(ntaxa + 1)
  out[1] <- paste(ntaxa, nsites)
  for (i in seq_along(seqs)) {
    out[i + 1] <- paste(names(seqs)[i], seqs[i])
  }
  writeLines(out, path)
}

#' Extract a contiguous range of sites from an alignment
#'
#' @param seqs Named character vector of sequences.
#' @param start First site (1-indexed).
#' @param end Last site (inclusive).
#' @return Named character vector with the sub-sequences.
subset_sites <- function(seqs, start, end) {
  vapply(seqs, function(s) substr(s, start, end), character(1),
         USE.NAMES = TRUE)
}

#' Concatenate multiple alignments by sites
#'
#' All alignments must share the same taxa in the same order.
#'
#' @param aln_list List of named character vectors (each an alignment).
#' @return Named character vector of the concatenated alignment.
concatenate_alignments <- function(aln_list) {
  taxa   <- names(aln_list[[1]])
  result <- rep("", length(taxa))
  names(result) <- taxa
  for (aln in aln_list) {
    result <- paste0(result, aln[taxa])
  }
  result
}
