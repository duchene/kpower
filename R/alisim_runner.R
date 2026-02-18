# AliSim simulation wrapper

#' Simulate replicate alignments under the best-fitting K model
#'
#' Attempts to extract the AliSim command directly from the ALISIM COMMAND
#' section of the `.iqtree` report file produced during the K_best empirical
#' fit. If that fails, constructs the AliSim call from the fitted model
#' parameters.
#'
#' @param fit_result The list returned by `fit_model()` for K_best.
#' @param n_sites Integer number of sites (should match the empirical
#'   alignment length).
#' @param B Integer number of replicate alignments to simulate.
#' @param outdir Directory in which to write simulated alignments.
#' @param seed Integer random seed for reproducibility.
#' @param mimic_gaps Logical; if TRUE (default) use the gap-mimicking AliSim
#'   command (`--alisim mimicked_MSA`) to reproduce the empirical gap pattern.
#'   If FALSE use the plain simulation command.
#' @param iqtree_bin Path to the IQ-TREE executable.
#' @param threads Number of threads (default `"AUTO"`).
#' @param timeout Seconds before the run is killed (default 7200).
#' @return Character vector of length B giving paths to the simulated FASTA
#'   files.
simulate_alignments <- function(fit_result, n_sites, B = 1000,
                                outdir = tempdir(), seed = 1,
                                mimic_gaps = TRUE,
                                iqtree_bin = find_iqtree(),
                                threads = "AUTO", timeout = 7200) {
  sim_dir    <- file.path(outdir, "simulations")
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  sim_prefix <- file.path(sim_dir, "sim")

  # --- Attempt to reuse the AliSim string from the .iqtree report ----------
  alisim_string <- parse_alisim_string(fit_result$iqtree_file,
                                       mimic_gaps = mimic_gaps)

  if (!is.null(alisim_string)) {
    message(
      "Using AliSim command from ALISIM COMMAND section of .iqtree report",
      if (mimic_gaps) " (gap-mimicking)." else " (plain simulation)."
    )
    args <- build_args_from_log_string(
      alisim_string, sim_prefix, B, seed, n_sites, threads
    )
  } else {
    message("ALISIM COMMAND section not found; constructing from fit parameters.")
    args <- build_alisim_args(
      fit_result, sim_prefix, n_sites, B, seed, threads
    )
  }

  run_iqtree(iqtree_bin, args, timeout = timeout)

  # AliSim writes files as sim_prefix_1.fa, sim_prefix_2.fa, ...
  # (exact suffix depends on IQ-TREE version; try .fa and .fasta)
  sim_files <- c(
    list.files(sim_dir, pattern = "^sim_[0-9]+\\.fa$",    full.names = TRUE),
    list.files(sim_dir, pattern = "^sim_[0-9]+\\.fasta$", full.names = TRUE)
  )
  sim_files <- sort(sim_files)

  if (length(sim_files) == 0) {
    stop("AliSim produced no output files in: ", sim_dir)
  }
  if (length(sim_files) != B) {
    warning("Expected ", B, " simulated files but found ", length(sim_files), ".")
  }

  sim_files
}

# Build AliSim args by modifying the string extracted from the log file.
# Replaces the original prefix, seed, and num-alignments with new values.
build_args_from_log_string <- function(alisim_string, sim_prefix, B, seed,
                                       n_sites, threads) {
  # Tokenise the command string (naive split on whitespace; handles most cases)
  tokens <- strsplit(trimws(alisim_string), "\\s+")[[1]]

  # Drop the binary name if present (first token often "iqtree2")
  if (grepl("iqtree", tokens[1], ignore.case = TRUE)) tokens <- tokens[-1]

  # Replace or append key flags with our values
  tokens <- replace_or_append(tokens, "--alisim",          sim_prefix)
  tokens <- replace_or_append(tokens, "--num-alignments",  as.character(B))
  tokens <- replace_or_append(tokens, "--seed",            as.character(seed))
  tokens <- replace_or_append(tokens, "--length",          as.character(n_sites))
  tokens <- replace_or_append(tokens, "-T",                as.character(threads))
  tokens <- replace_or_append(tokens, "--redo",            NULL, flag_only = TRUE)

  tokens
}

# Build AliSim args from scratch using the fit_result fields.
build_alisim_args <- function(fit_result, sim_prefix, n_sites, B, seed,
                              threads) {
  c(
    "--alisim",         sim_prefix,
    "-m",               fit_result$model_string,
    "-t",               fit_result$treefile,
    "--length",         as.character(n_sites),
    "--num-alignments", as.character(B),
    "--site-rate",      "SAMPLING",
    "--seed",           as.character(seed),
    "-T",               as.character(threads),
    "--redo"
  )
}

# Helper: in a token vector, find a flag and replace its next value,
# or append the flag+value if absent. If flag_only = TRUE, just ensure
# the flag is present with no value.
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
