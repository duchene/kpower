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

#' Parse tree weights from a MAST IQ-TREE report
#'
#' Extracts the `Tree weights:` line from a `.iqtree` report file produced
#' by a MAST run.
#'
#' @param iqtree_file Path to a `.iqtree` report file.
#' @return Numeric vector of tree weights, or NULL if not found.
parse_tree_weights <- function(iqtree_file) {
  lines <- readLines(iqtree_file)
  hit <- grep("^Tree weights:", lines, value = TRUE)
  if (length(hit) == 0) return(NULL)
  wstr <- sub("Tree weights:\\s*", "", hit[1])
  as.numeric(strsplit(wstr, ",\\s*")[[1]])
}

#' Parse MAST model parameters for AliSim reconstruction
#'
#' Extracts the substitution rate parameters, state frequencies, and
#' FreeRate / Gamma categories from a MAST `.iqtree` report.  These are the
#' *linked* (shared) model parameters used across all tree classes.
#'
#' @param iqtree_file Path to a `.iqtree` report file.
#' @return Named list with elements:
#'   \describe{
#'     \item{rate_params}{Numeric vector of 6 exchangeability parameters
#'       (A-C, A-G, A-T, C-G, C-T, G-T).}
#'     \item{freqs}{Numeric vector of 4 state frequencies
#'       (pi(A), pi(C), pi(G), pi(T)).}
#'     \item{rate_cats}{List of lists, each with \code{proportion} and
#'       \code{rate} (or NULL if no rate heterogeneity).}
#'   }
parse_mast_model_params <- function(iqtree_file) {
  lines <- readLines(iqtree_file)

  # --- Exchangeability parameters ---
  rate_names  <- c("A-C", "A-G", "A-T", "C-G", "C-T", "G-T")
  rate_params <- numeric(6)
  for (i in seq_along(rate_names)) {
    pat <- paste0("^\\s*", rate_names[i], ":\\s+([0-9.eE+-]+)")
    hit <- grep(pat, lines, value = TRUE)
    if (length(hit) > 0)
      rate_params[i] <- as.numeric(regmatches(hit[1], regexpr("[0-9.eE+-]+$", hit[1])))
  }

  # --- State frequencies ---
  freq_labels <- c("pi\\(A\\)", "pi\\(C\\)", "pi\\(G\\)", "pi\\(T\\)")
  freqs <- numeric(4)
  for (i in seq_along(freq_labels)) {
    pat <- paste0("^\\s*", freq_labels[i], "\\s*=\\s*([0-9.eE+-]+)")
    hit <- grep(pat, lines, value = TRUE)
    if (length(hit) > 0)
      freqs[i] <- as.numeric(sub(".*=\\s*", "", hit[1]))
  }

  # --- Rate heterogeneity categories ---
  # "Site proportion and rates:  (0.5368,0.007065) (0.2633,0.2372) ..."
  prop_line <- grep("^Site proportion and rates:", lines, value = TRUE)
  rate_cats <- NULL
  if (length(prop_line) > 0) {
    pairs <- regmatches(prop_line[1],
                        gregexpr("\\([0-9.eE+-]+,[0-9.eE+-]+\\)", prop_line[1]))[[1]]
    rate_cats <- lapply(pairs, function(p) {
      nums <- as.numeric(strsplit(gsub("[()]", "", p), ",")[[1]])
      list(proportion = nums[1], rate = nums[2])
    })
  }

  list(rate_params = rate_params, freqs = freqs, rate_cats = rate_cats)
}

#' Build a fully parameterised AliSim model string from parsed parameters
#'
#' Constructs a string like
#' `GTR{1.886/4.762/0.823/0.722/6.05/1}+FO{0.399/0.169/0.218/0.214}+R4{0.537/0.007/0.263/0.237/...}`.
#'
#' @param base_model Base substitution model (e.g. `"GTR"`).
#' @param params List returned by `parse_mast_model_params()`.
#' @return Character string suitable for `iqtree --alisim -m`.
build_alisim_model_string <- function(base_model, params) {
  ms <- paste0(base_model, "{", paste(params$rate_params, collapse = "/"), "}")
  ms <- paste0(ms, "+FU{", paste(params$freqs, collapse = "/"), "}")

  if (!is.null(params$rate_cats) && length(params$rate_cats) > 0) {
    K <- length(params$rate_cats)
    cat_strs <- vapply(params$rate_cats, function(cc) {
      paste(cc$proportion, cc$rate, sep = "/")
    }, character(1))
    ms <- paste0(ms, "+R", K, "{", paste(cat_strs, collapse = "/"), "}")
  }
  ms
}
