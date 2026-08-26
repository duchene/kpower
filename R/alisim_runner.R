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

  # --site-rate defaults to MEAN, which hands each site its shrunken posterior
  # mean rate from -s rather than drawing from the fitted rate distribution.
  # That under-produces constant sites, so replicates are harder to fit than the
  # empirical alignment and the bootstrap is no longer a valid null.
  tokens <- replace_or_append(tokens, "--site-rate",       "SAMPLING")

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

# ---------------------------------------------------------------------------
# MAST simulation: per-class AliSim + concatenation
# ---------------------------------------------------------------------------

#' Copy an alignment's gap pattern onto a simulated alignment
#'
#' AliSim's `-s` gap mimicking is used on the +R/+H path, but the MAST path
#' simulates each class at its own length, so a single `-s` template cannot line
#' up. Applying the mask positionally after concatenation reproduces the source
#' gaps exactly instead, and each window inherits its real counterpart's gap
#' profile.
#'
#' @param sim Named character vector: the simulated alignment.
#' @param gap_mask Logical matrix (taxa x sites) from `build_gap_mask()`.
#' @return `sim` with gaps written in.
#' @keywords internal
apply_gap_mask <- function(sim, gap_mask) {
  shared <- intersect(names(sim), rownames(gap_mask))
  if (length(shared) == 0L) return(sim)
  for (nm in shared) {
    ch <- strsplit(sim[[nm]], "")[[1]]
    g  <- gap_mask[nm, ]
    if (length(g) != length(ch)) next
    ch[g] <- "-"
    sim[[nm]] <- paste(ch, collapse = "")
  }
  sim
}

#' Build a taxa x sites logical gap mask from an alignment
#'
#' @param alignment Path to the source alignment.
#' @return Logical matrix with taxon names as rownames.
#' @keywords internal
build_gap_mask <- function(alignment) {
  a <- read_alignment(alignment)
  m <- do.call(rbind, strsplit(unname(a), ""))
  rownames(m) <- names(a)
  matrix(m %in% c("-", "?", "N", "n"), nrow = nrow(m), dimnames = dimnames(m))
}

#' Simulate replicate alignments under a fitted MAST model
#'
#' For each tree class in the MAST model, runs AliSim to generate B
#' full-length alignments. Each replicate then draws its own multinomial
#' split of the sites across classes (see `draw_class_counts()`) and takes
#' that many sites from each class, concatenated in class order, to produce B
#' full-length simulated alignments. Classes stay in contiguous blocks; only
#' the block sizes vary between replicates.
#'
#' @param mast_result Named list from the K_best MAST fit, containing at
#'   least: `treefile`, `iqtree_file`, `tree_weights`, `model_string`.
#' @param base_model Base substitution model (e.g. `"GTR"`).
#' @param n_sites Integer total number of sites in the empirical alignment.
#' @param B Integer number of replicate alignments.
#' @param seed Integer random seed (each class is offset by +k).
#' @param outdir Output directory.
#' @param iqtree_bin Path to IQ-TREE.
#' @param threads Number of threads.
#' @param timeout AliSim timeout in seconds.
#' @return Character vector of length B giving paths to the combined simulated
#'   alignment files.
simulate_mast_alignments <- function(mast_result, base_model, n_sites, B,
                                     seed, outdir, iqtree_bin, threads,
                                     timeout = 7200, alignment = NULL) {
  sim_dir <- file.path(outdir, "mast_simulations")
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

  weights <- mast_result$tree_weights
  K       <- length(weights)

  # Per-replicate class site counts, drawn from the fitted weights. Each
  # class is simulated at full length so any replicate's draw can be served
  # from it; the replicate takes the first n sites of each class.
  class_counts <- draw_class_counts(B, n_sites, weights, seed)   # K x B

  # Parse shared model parameters and build AliSim model string
  params    <- parse_mast_model_params(mast_result$iqtree_file)
  model_str <- build_alisim_model_string(base_model, params)

  # Per-class trees from the MAST treefile (one tree per line)
  all_trees <- readLines(mast_result$treefile)
  if (length(all_trees) != K) {
    stop("MAST treefile has ", length(all_trees),
         " trees but tree_weights has ", K, " entries.")
  }

  # --- Run AliSim once per class (each producing B alignments) ---------------
  class_dirs <- character(K)
  for (k in seq_len(K)) {
    class_dir <- file.path(sim_dir, paste0("class_", k))
    dir.create(class_dir, showWarnings = FALSE, recursive = TRUE)

    tree_file <- file.path(class_dir, "tree.newick")
    writeLines(all_trees[k], tree_file)

    sim_prefix <- file.path(class_dir, "sim")
    args <- c(
      "--alisim",         sim_prefix,
      "-m",               model_str,
      "-t",               tree_file,
      "--length",         as.character(n_sites),
      "--num-alignments", as.character(B),
      "--seed",           as.character(seed + k),
      "-T",               threads,
      "--redo"
    )
    message("  Simulating tree class ", k, " (weight ",
            round(weights[k], 4), ", ~",
            round(weights[k] * n_sites), " sites per replicate) ...")
    run_iqtree(iqtree_bin, args, timeout = timeout)
    class_dirs[k] <- class_dir
  }

  # --- Concatenate per-class replicates into full alignments -----------------
  combined_dir <- file.path(sim_dir, "combined")
  dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

  gap_mask <- if (is.null(alignment)) NULL else build_gap_mask(alignment)

  combined_files <- character(B)
  for (b in seq_len(B)) {
    counts <- class_counts[, b]
    class_alns <- lapply(seq_len(K), function(k) {
      # AliSim writes sim_1.phy, sim_2.phy, etc.
      f <- file.path(class_dirs[k], paste0("sim_", b, ".phy"))
      if (!file.exists(f)) f <- file.path(class_dirs[k], paste0("sim_", b, ".fa"))
      if (!file.exists(f)) stop("Missing simulated file for class ", k,
                                ", replicate ", b, ": ", f)
      # This replicate's share of class k: its first counts[k] sites. A class
      # can draw zero sites, which yields an empty (dropped) block.
      subset_sites(read_alignment(f), 1L, counts[k])
    })

    combined <- concatenate_alignments(class_alns)
    if (!is.null(gap_mask)) combined <- apply_gap_mask(combined, gap_mask)
    combined_files[b] <- file.path(combined_dir,
                                   paste0("sim_", pad_int(b), ".phy"))
    write_phylip(combined, combined_files[b])
  }

  combined_files
}
