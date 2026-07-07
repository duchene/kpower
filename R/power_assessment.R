# Power assessment: compare empirical IC profiles to simulation distributions

#' Fit all K models to a single alignment and return an IC table
#'
#' Used for both empirical and simulated alignments.
#'
#' @param alignment Path to the alignment file.
#' @param K_values Integer vector of K values to fit.
#' @param base_model Base substitution model string.
#' @param mix_type Mixture type suffix (e.g. `"+R"`).
#' @param fixed_tree Tree handling passed to `fit_model()`: `"NJ"`, a file
#'   path, or `NULL`.
#' @param outdir Output directory for IQ-TREE files.
#' @param label_prefix String prepended to the per-K label for `--prefix`.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param threads Number of threads.
#' @param timeout Per-run timeout in seconds.
#' @return Data frame with columns: K, lnL, df, AIC, AICc, BIC.
fit_all_K <- function(alignment, K_values, base_model = "GTR",
                      mix_type = "+R", fixed_tree = "NJ", outdir = tempdir(),
                      label_prefix = "", iqtree_bin = find_iqtree(),
                      threads = "1", timeout = Inf) {
  results <- lapply(K_values, function(K) {
    label <- paste0(label_prefix, "K", K)
    fit   <- fit_model(
      alignment  = alignment,
      K          = K,
      base_model = base_model,
      mix_type   = mix_type,
      fixed_tree = fixed_tree,
      outdir     = outdir,
      label      = label,
      iqtree_bin = iqtree_bin,
      threads    = threads,
      timeout    = timeout
    )
    data.frame(
      K    = fit$K,
      lnL  = fit$lnL,
      df   = fit$df,
      AIC  = fit$AIC,
      AICc = fit$AICc,
      BIC  = fit$BIC,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

#' Run the full parametric bootstrap power assessment
#'
#' For each of B simulated alignments, fits all K models and records IC
#' scores. Returns a long-format data frame of all simulation results,
#' alongside the power estimate.
#'
#' @param sim_files Character vector of paths to simulated alignments (from
#'   `simulate_alignments()`).
#' @param K_values Integer vector of K values to fit on each simulation.
#' @param K_best Integer; the K selected from empirical data.
#' @param ic Character; which IC to use for power calculation: `"AIC"`,
#'   `"AICc"`, or `"BIC"` (default `"BIC"`).
#' @param base_model Base substitution model string.
#' @param mix_type Mixture type suffix.
#' @param fixed_tree Tree handling for simulation refits: `"NJ"`, a file
#'   path, or `NULL`. Default `"NJ"`.
#' @param outdir Output directory for IQ-TREE files.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param threads Number of threads.
#' @param n_cores Number of cores for parallel simulation refits (default 1).
#' @param timeout Per-run timeout in seconds.
#' @return List with:
#'   \describe{
#'     \item{sim_ic}{Long-format data frame: replicate, K, lnL, AIC, AICc,
#'       BIC}
#'     \item{power}{Proportion of replicates that select K_best under `ic`}
#'   }
assess_power <- function(sim_files, K_values, K_best, ic = "BIC",
                         base_model = "GTR", mix_type = "+R",
                         fixed_tree = "NJ", outdir = tempdir(),
                         iqtree_bin = find_iqtree(), threads = "1",
                         n_cores = 1, timeout = Inf) {
  sim_outdir <- file.path(outdir, "sim_fits")
  dir.create(sim_outdir, showWarnings = FALSE, recursive = TRUE)

  run_one <- function(b) {
    label_prefix <- paste0("sim", pad_int(b), "_")
    tbl <- fit_all_K(
      alignment    = sim_files[b],
      K_values     = K_values,
      base_model   = base_model,
      mix_type     = mix_type,
      fixed_tree   = fixed_tree,
      outdir       = sim_outdir,
      label_prefix = label_prefix,
      iqtree_bin   = iqtree_bin,
      threads      = threads,
      timeout      = timeout
    )
    tbl$replicate <- b
    tbl
  }
  run_one_safe <- function(b) tryCatch(run_one(b), error = function(e) e)

  sim_results <- run_parallel(
    seq_along(sim_files), run_one_safe, n_cores = n_cores
  )

  # failed replicates come back as error conditions; detect and report them
  failed <- vapply(sim_results, inherits, logical(1), "error")
  if (any(failed)) {
    idx <- which(failed)
    stop(
      "Bootstrap refit(s) failed for replicate(s): ", paste(idx, collapse = ", "),
      "\nFirst error: ", conditionMessage(sim_results[[idx[1]]])
    )
  }

  sim_ic <- do.call(rbind, sim_results)

  # For each replicate, which K minimises the chosen IC?
  best_per_rep <- tapply(
    sim_ic[[ic]], sim_ic$replicate,
    function(x) K_values[which.min(x)]
  )
  power <- mean(best_per_rep == K_best)

  list(sim_ic = sim_ic, power = power)
}

#' Run the parametric bootstrap power assessment for MAST models
#'
#' For each simulated alignment, fits K = 1 (single-tree BioNJ) and
#' K = 2..K_max (MAST with top-K ranked trees).
#'
#' @param sim_files Character vector of paths to simulated alignments.
#' @param K_values Integer vector of K values.
#' @param K_best Integer; K selected from empirical data.
#' @param ic Character; IC for power calculation.
#' @param base_model Base substitution model string.
#' @param rate_model Rate heterogeneity string (e.g. `"+R3"`).
#' @param tree_files Named list of tree file paths per K (from
#'   `build_mast_tree_files()`).
#' @param unlinked Logical; if TRUE, use MIX syntax for unlinked per-tree
#'   substitution parameters.
#' @param outdir Output directory.
#' @param iqtree_bin Path to IQ-TREE.
#' @param threads Number of threads.
#' @param n_cores Number of cores for parallel replicates.
#' @param timeout Per-run timeout in seconds.
#' @return List with `sim_ic` (data frame) and `power` (numeric).
assess_mast_power <- function(sim_files, K_values, K_best, ic = "BIC",
                              base_model, rate_model, tree_files,
                              unlinked = FALSE, fixed_tree = "NJ",
                              outdir, iqtree_bin, threads,
                              n_cores = 1, timeout = Inf) {
  sim_outdir <- file.path(outdir, "mast_sim_fits")
  dir.create(sim_outdir, showWarnings = FALSE, recursive = TRUE)

  run_one <- function(b) {
    label_prefix <- paste0("sim", pad_int(b), "_")
    tbl <- fit_mast_all_K(
      alignment    = sim_files[b],
      K_values     = K_values,
      base_model   = base_model,
      rate_model   = rate_model,
      tree_files   = tree_files,
      unlinked     = unlinked,
      fixed_tree   = fixed_tree,
      outdir       = sim_outdir,
      label_prefix = label_prefix,
      iqtree_bin   = iqtree_bin,
      threads      = threads,
      timeout      = timeout
    )
    tbl$replicate <- b
    tbl
  }
  run_one_safe <- function(b) tryCatch(run_one(b), error = function(e) e)

  sim_results <- run_parallel(
    seq_along(sim_files), run_one_safe, n_cores = n_cores
  )

  failed <- vapply(sim_results, inherits, logical(1), "error")
  if (any(failed)) {
    idx <- which(failed)
    stop(
      "MAST bootstrap refit(s) failed for replicate(s): ",
      paste(idx, collapse = ", "),
      "\nFirst error: ", conditionMessage(sim_results[[idx[1]]])
    )
  }

  sim_ic <- do.call(rbind, sim_results)

  best_per_rep <- tapply(
    sim_ic[[ic]], sim_ic$replicate,
    function(x) K_values[which.min(x)]
  )
  power <- mean(best_per_rep == K_best)

  list(sim_ic = sim_ic, power = power)
}

#' Compute K_best and power under all three information criteria
#'
#' For each of AIC, AICc, and BIC: selects K_best from the empirical IC
#' profile, then computes power as the proportion of bootstrap replicates
#' that recover that K_best.
#'
#' @param empirical_ic Data frame with columns K, AIC, AICc, BIC from
#'   empirical fits.
#' @param sim_ic Long-format data frame with replicate, K, AIC, AICc, BIC.
#' @param K_values Integer vector of K values.
#' @return Named list with elements AIC, AICc, BIC, each containing
#'   `K_best` (integer) and `power` (numeric between 0 and 1).
compute_power_all_ic <- function(empirical_ic, sim_ic, K_values) {
  ics <- c("AIC", "AICc", "BIC")
  result <- lapply(stats::setNames(ics, ics), function(crit) {
    K_best <- K_values[which.min(empirical_ic[[crit]])]
    best_per_rep <- tapply(
      sim_ic[[crit]], sim_ic$replicate,
      function(x) K_values[which.min(x)]
    )
    list(K_best = K_best, power = mean(best_per_rep == K_best))
  })
  result
}
