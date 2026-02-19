# Main entry point

#' Assess alignment power to discriminate phylogenetic mixture model categories
#'
#' Fits mixture models with K = 1, 2, ..., K_max categories to an empirical
#' alignment, identifies the best K via information criteria, simulates B
#' replicate alignments from that model using IQ-TREE's AliSim, refits all K
#' values on each replicate, and returns a power estimate alongside a
#' publication-ready figure.
#'
#' @param alignment Path to the input alignment (FASTA or PHYLIP).
#' @param K_max Integer; maximum number of mixture categories to evaluate.
#' @param K_min Integer; minimum K (default 1).
#' @param base_model Base substitution model string (e.g. `"GTR"`, `"LG"`).
#'   Default `"GTR"`.
#' @param mix_type Mixture family: `"+R"` (FreeRate, default), `"+H"` (GHOST).
#' @param ic Information criterion used to select K_best and compute power:
#'   `"AIC"`, `"AICc"`, or `"BIC"` (default `"BIC"`).
#' @param fixed_tree Tree handling: `"NJ"` (per-K BioNJ, default), a path to
#'   a fixed tree file, or `NULL` (full heuristic search).
#' @param B Integer number of parametric bootstrap replicates (default 1000).
#' @param seed Integer random seed passed to AliSim (default 1).
#' @param outdir Directory for all IQ-TREE output files. Defaults to a
#'   temporary directory.
#' @param iqtree_bin Path to the IQ-TREE executable. Detected automatically
#'   if not supplied.
#' @param n_cores Number of cores controlling both R-level parallel bootstrap
#'   refits and IQ-TREE's `-T` thread count (default 1).
#' @param timeout Per-run timeout in seconds (default 3600).
#'
#' @return An object of class `kpower_result`, a list containing:
#'   \describe{
#'     \item{empirical}{Data frame of K, lnL, df, AIC, AICc, BIC for the
#'       empirical alignment.}
#'     \item{sim_ic}{Long-format data frame of replicate, K, and IC scores
#'       from all bootstrap fits.}
#'     \item{K_best}{Integer; K selected from empirical data.}
#'     \item{power}{Numeric; proportion of simulations that recover K_best.}
#'     \item{ic}{Character; the IC used.}
#'     \item{plot}{A ggplot2 object (IC profile figure).}
#'   }
#' @export
kpower <- function(alignment,
                   K_max,
                   K_min      = 1L,
                   base_model = "GTR",
                   mix_type   = "+R",
                   ic         = "BIC",
                   fixed_tree = "NJ",
                   B          = 1000L,
                   seed       = 1L,
                   outdir     = tempdir(),
                   iqtree_bin = find_iqtree(),
                   n_cores    = 1L,
                   timeout    = 3600L) {

  ic      <- match.arg(ic, c("AIC", "AICc", "BIC"))
  K_values <- seq.int(K_min, K_max)
  threads  <- as.character(n_cores)   # n_cores drives IQ-TREE -T and mclapply

  emp_outdir <- file.path(outdir, "empirical")
  dir.create(emp_outdir, showWarnings = FALSE, recursive = TRUE)

  # --- Step 1: Fit all K to empirical data ----------------------------------
  message("Fitting K = ", K_min, " to ", K_max, " on empirical alignment ...")
  empirical_ic <- fit_all_K(
    alignment    = alignment,
    K_values     = K_values,
    base_model   = base_model,
    mix_type     = mix_type,
    fixed_tree   = fixed_tree,
    outdir       = emp_outdir,
    label_prefix = "empirical_",
    iqtree_bin   = iqtree_bin,
    threads      = threads,
    timeout      = timeout
  )

  # --- Step 2: Select K_best ------------------------------------------------
  K_best     <- K_values[which.min(empirical_ic[[ic]])]
  best_label  <- paste0("empirical_K", K_best)
  best_prefix <- make_prefix(emp_outdir, best_label)
  best_fit <- list(
    K            = K_best,
    model_string = if (K_best == 1) base_model else
                     paste0(base_model, mix_type, K_best),
    treefile     = paste0(best_prefix, ".treefile"),
    logfile      = paste0(best_prefix, ".log"),
    iqtree_file  = paste0(best_prefix, ".iqtree")
  )

  message("K_best = ", K_best, " (selected by ", ic, ")")

  # --- Step 3: Simulate B alignments under K_best ---------------------------
  n_sites <- alignment_length(alignment)
  message("Simulating ", B, " alignments of ", n_sites, " sites via AliSim ...")
  sim_files <- simulate_alignments(
    fit_result = best_fit,
    alignment  = alignment,
    n_sites    = n_sites,
    B          = B,
    outdir     = outdir,
    seed       = seed,
    iqtree_bin = iqtree_bin,
    threads    = threads
  )

  # --- Step 4: Refit all K on each simulated alignment ---------------------
  message("Refitting K = ", K_min, " to ", K_max,
          " on ", B, " simulated alignments ...")
  power_result <- assess_power(
    sim_files  = sim_files,
    K_values   = K_values,
    K_best     = K_best,
    ic         = ic,
    base_model = base_model,
    mix_type   = mix_type,
    fixed_tree = fixed_tree,
    outdir     = outdir,
    iqtree_bin = iqtree_bin,
    threads    = threads,
    n_cores    = n_cores,
    timeout    = timeout
  )

  message(sprintf(
    "Power: %.1f%% of simulations recover K_best = %d under %s",
    power_result$power * 100, K_best, ic
  ))

  # --- Step 5: Build figure -------------------------------------------------
  fig <- plot_kpower(
    empirical_ic = empirical_ic,
    sim_ic       = power_result$sim_ic,
    K_best       = K_best,
    power        = power_result$power,
    ic           = ic
  )

  structure(
    list(
      empirical = empirical_ic,
      sim_ic    = power_result$sim_ic,
      K_best    = K_best,
      power     = power_result$power,
      ic        = ic,
      plot      = fig
    ),
    class = "kpower_result"
  )
}

#' Print method for kpower_result
#' @export
print.kpower_result <- function(x, ...) {
  cat("kpower result\n")
  cat("  IC used  :", x$ic, "\n")
  cat("  K_best   :", x$K_best, "\n")
  cat(sprintf("  Power    : %.1f%%\n", x$power * 100))
  invisible(x)
}
