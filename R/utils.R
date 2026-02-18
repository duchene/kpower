# Utility functions for kpower

#' Find the IQ-TREE executable
#'
#' Searches PATH and then the user option `kpower.iqtree_path`.
#'
#' @return Path to the IQ-TREE binary as a character string.
#' @export
find_iqtree <- function() {
  # Check user-set option first
  opt <- getOption("kpower.iqtree_path")
  if (!is.null(opt) && file.exists(opt)) return(opt)

  # Search PATH for iqtree3, iqtree2, then iqtree
  for (bin in c("iqtree3", "iqtree2", "iqtree")) {
    path <- Sys.which(bin)
    if (nzchar(path)) return(path)
  }

  stop(
    "IQ-TREE executable not found. Install IQ-TREE and ensure it is on PATH, ",
    "or set options(kpower.iqtree_path = '/path/to/iqtree2')."
  )
}

#' Check the IQ-TREE installation and report its version
#'
#' Calls `iqtree2 --version` and prints the version string. Useful for
#' verifying the installation before running a full analysis.
#'
#' @param iqtree_bin Path to the IQ-TREE executable (auto-detected if NULL).
#' @return Invisibly, the version string.
#' @export
check_iqtree <- function(iqtree_bin = NULL) {
  if (is.null(iqtree_bin)) iqtree_bin <- find_iqtree()
  result <- processx::run(iqtree_bin, "--version", error_on_status = FALSE)
  version_line <- strsplit(result$stdout, "\n")[[1]][1]
  message("IQ-TREE found: ", iqtree_bin)
  message(version_line)
  invisible(version_line)
}

#' Get alignment length from a FASTA file
#'
#' @param fasta Path to a FASTA alignment file.
#' @return Integer number of sites.
alignment_length <- function(fasta) {
  lines <- readLines(fasta)
  seq_lines <- lines[!startsWith(lines, ">")]
  nchar(paste(seq_lines[seq_lines != ""], collapse = ""))
}

#' Build a unique run prefix inside a directory
#'
#' @param outdir Base output directory.
#' @param label Short label string (e.g. "empirical_K4" or "sim_003_K2").
#' @return Full path prefix string suitable for IQ-TREE --prefix.
make_prefix <- function(outdir, label) {
  dir.create(file.path(outdir, label), showWarnings = FALSE, recursive = TRUE)
  file.path(outdir, label, label)
}

#' Zero-pad an integer for consistent file naming
#'
#' @param i Integer.
#' @param width Total character width (default 4).
#' @return Zero-padded character string.
pad_int <- function(i, width = 4) {
  formatC(i, width = width, flag = "0")
}
