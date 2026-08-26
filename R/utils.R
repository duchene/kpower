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

#' Get alignment length from a FASTA or PHYLIP alignment file
#'
#' @param alignment Path to a FASTA or PHYLIP alignment file.
#' @return Integer number of sites.
alignment_length <- function(alignment) {
  seqs <- read_alignment(alignment)
  nchar(seqs[1])
}

# Null-coalesce operator (base R does not have %||% before 4.3)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Run a function over an index in parallel via future.apply
#'
#' Replaces `parallel::mclapply()`. Forking with `mclapply()` installs a
#' SIGCHLD handler in the main process that conflicts with `processx::run()`
#' (used for every IQ-TREE call, including from the main process itself),
#' leaking zombie IQ-TREE/R processes and occasionally hanging `processx`
#' calls that rely on SIGCHLD-based exit detection. `future::multisession`
#' uses separate background R processes instead of forking, so it doesn't
#' touch the main process's signal handlers.
#'
#' @param X Vector or list to iterate over.
#' @param FUN Function of one argument, applied to each element of `X`.
#'   Should not throw for expected failure modes; wrap in `tryCatch()`
#'   upstream and return the condition object instead, if per-element
#'   failures need to be detected downstream (as `assess_power()` and
#'   `assess_mast_power()` do).
#' @param n_cores Number of parallel workers. Runs sequentially via
#'   `lapply()` when `n_cores <= 1`.
#' @return A list of results, one per element of `X`.
run_parallel <- function(X, FUN, n_cores = 1) {
  if (n_cores > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    old_plan <- future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(old_plan), add = TRUE)
    future.apply::future_lapply(X, FUN)
  } else {
    lapply(X, FUN)
  }
}

# Format an elapsed time (seconds) as a compact human-readable string,
# e.g. 4.2s, 3m 12s, 1h 05m. Used for survey timing output.
format_elapsed <- function(seconds) {
  if (is.na(seconds)) return("NA")
  if (seconds < 60) return(sprintf("%.1fs", seconds))
  if (seconds < 3600) {
    m <- floor(seconds / 60)
    s <- round(seconds - 60 * m)
    return(sprintf("%dm %02ds", m, s))
  }
  h <- floor(seconds / 3600)
  m <- round((seconds - 3600 * h) / 60)
  sprintf("%dh %02dm", h, m)
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

#' Build an IQ-TREE model string for a given K and mixture type
#'
#' For K = 1, returns the base model alone.  For K > 1, appends the mixture
#' suffix.  The `*R` type is special: IQ-TREE implements it via
#' `MIX{base+F,...}*R{K}` (each rate class gets its own substitution model).
#'
#' @param base_model Base substitution model (e.g. `"GTR"`).
#' @param mix_type Mixture type: `"+R"`, `"*R"`, `"+H"`, or `"*H"`.
#' @param K Integer number of mixture categories.
#' @return Character string suitable for IQ-TREE `-m`.
build_model_str <- function(base_model, mix_type, K) {
  if (K == 1) return(base_model)
  paste0(base_model, mix_type, K)
}

#' Build an IQ-TREE model string for MAST tree-mixture models
#'
#' For linked MAST (`+T`): `base+FO+rate+T`.
#' For unlinked MAST (`*T`): `MIX{base+FO,...,base+FO}+rate+T`, giving each
#' tree class its own substitution parameters.
#'
#' @param base_model Base substitution model.
#' @param rate_model Rate heterogeneity string (e.g. `"+R3"`).
#' @param K Number of tree classes.
#' @param unlinked Logical; if TRUE, use `MIX{...}` syntax.
#' @return Character string for IQ-TREE `-m`.
build_mast_model_str <- function(base_model, rate_model, K, unlinked = FALSE) {
  if (K <= 1) {
    return(paste0(base_model, "+FO", rate_model))
  }
  if (unlinked) {
    component <- paste0(base_model, "+FO")
    mix_part  <- paste0("MIX{", paste(rep(component, K), collapse = ","), "}")
    paste0(mix_part, rate_model, "+T")
  } else {
    paste0(base_model, "+FO", rate_model, "+T")
  }
}

#' Draw per-replicate class site counts from MAST tree weights
#'
#' A single deterministic `round(weights * n_sites)` gives every replicate the
#' same class sizes; with all other simulation inputs fixed the B replicates
#' then come out identical, i.e. a bootstrap carrying no sampling variation.
#' Each replicate instead draws its own multinomial split of the sites across
#' classes. Reproducible given `seed`, and the caller's RNG stream is restored
#' on exit.
#'
#' @param B Number of replicates.
#' @param n_sites Total sites per replicate; every column sums to this.
#' @param weights Numeric vector of K class weights.
#' @param seed Integer seed.
#' @return Integer matrix, K rows x B columns.
#' @keywords internal
draw_class_counts <- function(B, n_sites, weights, seed) {
  had_seed <- exists(".Random.seed", envir = globalenv())
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv()))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)

  set.seed(seed)
  stats::rmultinom(B, size = n_sites, prob = weights)
}
