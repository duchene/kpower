# Functions to parse IQ-TREE output files

#' Parse an IQ-TREE report file (.iqtree)
#'
#' Extracts log-likelihood, AIC, AICc, BIC, and number of free parameters.
#'
#' @param iqtree_file Path to a `.iqtree` report file.
#' @return Named list with elements: lnL, df, AIC, AICc, BIC.
parse_iqtree_report <- function(iqtree_file) {
  lines <- readLines(iqtree_file)

  extract <- function(pattern) {
    hit <- grep(pattern, lines, value = TRUE)
    if (length(hit) == 0) return(NA_real_)
    as.numeric(sub(paste0(".*", pattern, "\\s*([0-9eE.+-]+).*"), "\\1", hit[1]))
  }

  list(
    lnL  = extract("Log-likelihood of the tree:"),
    df   = extract("Number of free parameters.*:"),
    AIC  = extract("Akaike information criterion \\(AIC\\) score:"),
    AICc = extract("Corrected Akaike information criterion \\(AICc\\) score:"),
    BIC  = extract("Bayesian information criterion \\(BIC\\) score:")
  )
}

#' Extract the plain AliSim command string from an IQ-TREE report file
#'
#' IQ-TREE appends an "ALISIM COMMAND" section at the end of every `.iqtree`
#' report file. This function extracts the plain simulation command (starting
#' with `--alisim simulated_MSA`), which contains the fully-specified fitted
#' model string, tree path, and alignment length. `simulate_alignments()` then
#' blends this with gap mimicking by replacing the output prefix and appending
#' `-s original_alignment`.
#'
#' @param iqtree_file Path to a `.iqtree` report file.
#' @return Character string of the plain AliSim argument line, or NULL if the
#'   section is not found (caller falls back to constructing the command).
parse_alisim_string <- function(iqtree_file) {
  if (!file.exists(iqtree_file)) return(NULL)

  lines <- readLines(iqtree_file)

  section_idx <- grep("^ALISIM COMMAND", lines)
  if (length(section_idx) == 0) return(NULL)

  after <- lines[seq(section_idx[1] + 1, length(lines))]

  # The plain command starts directly with "--alisim simulated_MSA ..."
  # (no binary prefix) and carries the full -m and -t specification.
  target_idx <- grep("^--alisim\\s+simulated", after)
  if (length(target_idx) == 0) target_idx <- grep("^--alisim", after)
  if (length(target_idx) == 0) return(NULL)

  raw <- trimws(after[target_idx[1]])
  if (!nzchar(raw)) return(NULL)
  raw
}
