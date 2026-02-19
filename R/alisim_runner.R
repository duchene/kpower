# AliSim simulation wrapper

#' Simulate replicate alignments under the best-fitting K model
#'
#' Extracts the plain AliSim command from the ALISIM COMMAND section of the
#' `.iqtree` report (which carries the full fitted model string and tree path),
#' then blends it with gap mimicking by appending `-s original_alignment`.
#' This combines the fully-specified model from the plain command with the
#' gap-pattern reproduction from the mimicked approach. Falls back to
#' constructing the command from scratch if the section is not found.
#'
#' @param fit_result The list returned by `fit_model()` for K_best.
#' @param alignment Path to the original empirical alignment, used for gap
#'   mimicking via `-s`.
#' @param n_sites Integer number of sites (matches the empirical alignment).
#' @param B Integer number of replicate alignments to simulate.
#' @param outdir Directory in which to write simulated alignments.
#' @param seed Integer random seed for reproducibility.
#' @param iqtree_bin Path to the IQ-TREE executable.
#' @param threads Number of threads (default `"1"`).
#' @param timeout Seconds before the run is killed (default 7200).
#' @return Character vector of length B giving paths to the simulated files.
simulate_alignments <- function(fit_result, alignment, n_sites, B = 1000,
                                outdir = tempdir(), seed = 1,
                                iqtree_bin = find_iqtree(),
                                threads = "1", timeout = 7200) {
  sim_dir    <- file.path(outdir, "simulations")
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  sim_prefix <- file.path(sim_dir, "sim")

  alisim_string <- parse_alisim_string(fit_result$iqtree_file)

  if (!is.null(alisim_string)) {
    message("Using AliSim command from ALISIM COMMAND section of .iqtree report.")
    args <- build_args_from_alisim_string(
      alisim_string, sim_prefix, alignment, B, seed, n_sites, threads
    )
  } else {
    message("ALISIM COMMAND section not found; constructing from fit parameters.")
    args <- build_alisim_args(
      fit_result, alignment, sim_prefix, n_sites, B, seed, threads
    )
  }

  run_iqtree(iqtree_bin, args, timeout = timeout)

  # AliSim writes files as sim_1.phy, sim_2.phy etc. (or .fa / .fasta)
  sim_files <- sort(list.files(
    sim_dir,
    pattern     = "^sim_[0-9]+\\.(phy|fa|fasta)$",
    full.names  = TRUE
  ))

  if (length(sim_files) == 0) {
    stop("AliSim produced no output files in: ", sim_dir)
  }
  if (length(sim_files) != B) {
    warning("Expected ", B, " simulated files but found ", length(sim_files), ".")
  }

  sim_files
}

# Build AliSim args from the plain command in the .iqtree ALISIM COMMAND
# section, replacing the prefix with sim_prefix and appending -s alignment
# for gap mimicking.
build_args_from_alisim_string <- function(alisim_string, sim_prefix, alignment,
                                          B, seed, n_sites, threads) {
  # Tokenise respecting double-quoted strings (e.g. the -m model string).
  # processx passes arguments directly to the binary without a shell, so
  # quotes must be stripped rather than preserved.
  tokens <- tokenise_command(alisim_string)

  # Drop binary name if present
  if (grepl("^iqtree", tokens[1], ignore.case = TRUE)) tokens <- tokens[-1]

  # Replace --alisim value (simulated_MSA) with our output prefix
  tokens <- replace_or_append(tokens, "--alisim",          sim_prefix)
  tokens <- replace_or_append(tokens, "--num-alignments",  as.character(B))
  tokens <- replace_or_append(tokens, "--seed",            as.character(seed))
  tokens <- replace_or_append(tokens, "--length",          as.character(n_sites))
  tokens <- replace_or_append(tokens, "-T",                as.character(threads))
  tokens <- replace_or_append(tokens, "--redo",            NULL, flag_only = TRUE)

  # Add the original alignment for gap mimicking (-s)
  tokens <- replace_or_append(tokens, "-s",                alignment)

  tokens
}

# Build AliSim args from scratch when no ALISIM COMMAND section is found.
build_alisim_args <- function(fit_result, alignment, sim_prefix, n_sites, B,
                              seed, threads) {
  c(
    "--alisim",         sim_prefix,
    "-m",               fit_result$model_string,
    "-t",               fit_result$treefile,
    "-s",               alignment,
    "--length",         as.character(n_sites),
    "--num-alignments", as.character(B),
    "--site-rate",      "SAMPLING",
    "--seed",           as.character(seed),
    "-T",               as.character(threads),
    "--redo"
  )
}

# Tokenise a command string respecting double-quoted substrings.
# Quoted tokens have their surrounding double quotes stripped so that
# processx can pass them directly to the binary without shell interpretation.
tokenise_command <- function(s) {
  tokens  <- character(0)
  current <- ""
  in_quote <- FALSE
  chars <- strsplit(trimws(s), "")[[1]]

  for (ch in chars) {
    if (ch == '"') {
      in_quote <- !in_quote   # toggle; do not append the quote character
    } else if (ch == " " && !in_quote) {
      if (nzchar(current)) {
        tokens  <- c(tokens, current)
        current <- ""
      }
    } else {
      current <- paste0(current, ch)
    }
  }
  if (nzchar(current)) tokens <- c(tokens, current)
  tokens
}

# Helper: find a flag in a token vector and replace its next value,
# or append flag+value if absent. If flag_only = TRUE, just ensure the
# flag itself is present with no value.
replace_or_append <- function(tokens, flag, value, flag_only = FALSE) {
  idx <- which(tokens == flag)
  if (length(idx) > 0) {
    if (!flag_only && !is.null(value)) tokens[idx[1] + 1] <- value
  } else {
    if (flag_only) {
      tokens <- c(tokens, flag)
    } else if (!is.null(value)) {
      tokens <- c(tokens, flag, value)
    }
  }
  tokens
}
