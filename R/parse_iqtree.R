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

#' Extract the AliSim command string from an IQ-TREE report file (.iqtree)
#'
#' IQ-TREE appends an "ALISIM COMMAND" section at the end of every `.iqtree`
#' report file containing two ready-to-use invocations: one for a plain
#' simulation and one that mimics the empirical alignment's gap pattern. This
#' function locates that section and returns the appropriate command string,
#' avoiding the need to reconstruct the model specification manually (which is
#' error-prone for complex mixture model strings).
#'
#' @param iqtree_file Path to a `.iqtree` report file.
#' @param mimic_gaps Logical; if TRUE (default) return the command that mimics
#'   the empirical gap pattern (`--alisim mimicked_MSA ...`). If FALSE return
#'   the plain simulation command (`--alisim simulated_MSA ...`).
#' @return Character string of the full AliSim argument line, or NULL if the
#'   section is not found (caller falls back to constructing the command).
parse_alisim_string <- function(iqtree_file, mimic_gaps = TRUE) {
  if (!file.exists(iqtree_file)) return(NULL)

  lines <- readLines(iqtree_file)

  # Locate the ALISIM COMMAND section header
  section_idx <- grep("^ALISIM COMMAND", lines)
  if (length(section_idx) == 0) return(NULL)

  # Work only with lines after the section header
  after <- lines[seq(section_idx[1] + 1, length(lines))]

  # Find all lines that contain "--alisim" within this section
  hit_idx <- grep("--alisim", after)
  if (length(hit_idx) == 0) return(NULL)

  # IQ-TREE writes two --alisim lines:
  #   1. "--alisim simulated_MSA ..."  — plain simulation
  #   2. "--alisim mimicked_MSA ..."   — gap-mimicking simulation
  # Select the appropriate one based on mimic_gaps.
  if (mimic_gaps) {
    target_idx <- grep("--alisim\\s+mimicked", after)
    if (length(target_idx) == 0) target_idx <- hit_idx  # fall back to first
  } else {
    target_idx <- grep("--alisim\\s+simulated", after)
    if (length(target_idx) == 0) target_idx <- hit_idx
  }

  raw <- trimws(after[target_idx[1]])
  if (!nzchar(raw)) return(NULL)
  raw
}
