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

#' K minimising an IC within each replicate
#'
#' Reads K from each replicate's own `K` column rather than indexing a shared
#' `K_values` vector positionally, because a replicate's K range can be
#' shorter than the requested one (see `mast_candidates()`, which reduces
#' K_max to the number of usable candidate trees).
#'
#' @param sim_ic Long-format data frame with `replicate`, `K` and IC columns.
#' @param ic Name of the IC column.
#' @return Named numeric vector of best K per replicate; `NA` for a replicate
#'   with no finite IC value.
best_K_per_replicate <- function(sim_ic, ic) {
  vapply(
    split(seq_len(nrow(sim_ic)), sim_ic$replicate),
    function(idx) {
      vals <- sim_ic[[ic]][idx]
      ok   <- is.finite(vals)
      if (!any(ok)) return(NA_real_)
      as.numeric(sim_ic$K[idx][ok][which.min(vals[ok])])
    },
    numeric(1)
  )
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

  # failed replicates come back as error conditions; drop them (with a warning)
  # rather than aborting the whole family, so a single bad refit does not discard
  # every successful replicate. Only error if nothing survives.
  failed <- vapply(sim_results, inherits, logical(1), "error")
  if (any(failed)) {
    idx <- which(failed)
    warning(
      "Bootstrap refit(s) failed for ", sum(failed), "/", length(failed),
      " replicate(s) (dropped): ", paste(idx, collapse = ", "),
      "\nFirst error: ", conditionMessage(sim_results[[idx[1]]])
    )
  }
  sim_results <- sim_results[!failed]
  if (length(sim_results) == 0) {
    stop("All bootstrap replicates failed; cannot compute power.")
  }

  sim_ic <- do.call(rbind, sim_results)

  # For each replicate, which K minimises the chosen IC?
  best_per_rep <- best_K_per_replicate(sim_ic, ic)
  power <- mean(best_per_rep == K_best, na.rm = TRUE)

  list(sim_ic = sim_ic, power = power)
}

#' Run the parametric bootstrap power assessment for MAST models
#'
#' Each replicate derives its own candidate trees from its own alignment via
#' `mast_candidates()`, then fits K = 1 (single-tree BioNJ) and K >= 2 (MAST
#' with its own top-K ranked trees).
#'
#' The empirical tree sets are deliberately not passed in. Reusing them would
#' hand each replicate the exact topologies it was simulated from, pinning the
#' IC minimum at K_best and forcing power to 100%.
#'
#' @param sim_files Character vector of paths to simulated alignments.
#' @param K_values Integer vector of K values.
#' @param K_best Integer; K selected from empirical data.
#' @param ic Character; IC for power calculation.
#' @param base_model Base substitution model string.
#' @param rate_model Rate heterogeneity string (e.g. `"+R3"`), or `NULL` to
#'   re-derive it per replicate from that replicate's windows.
#' @param unlinked Logical; if TRUE, use MIX syntax for unlinked per-tree
#'   substitution parameters.
#' @param window_method Window tree estimator: `"NJ"` or `"fast"`.
#' @param filter_trees Logical; passed to `mast_candidates()`. With `TRUE`
#'   (default) a replicate whose windows yield fewer usable trees fits a
#'   shorter K range than requested, and is recorded as such in `sim_ic`.
#'   No floor is applied: a replicate that cannot reach the empirical K_best
#'   is scored as a miss.
#' @param fixed_tree Tree handling for the K = 1 fit.
#' @param outdir Output directory.
#' @param iqtree_bin Path to IQ-TREE.
#' @param threads Number of threads.
#' @param n_cores Number of cores for parallel replicates.
#' @param seed Optional base seed for per-replicate window tree searches.
#' @param timeout Per-run timeout in seconds.
#' @return List with `sim_ic` (data frame) and `power` (numeric).
assess_mast_power <- function(sim_files, K_values, K_best, ic = "BIC",
                              base_model, rate_model = NULL,
                              unlinked = FALSE, window_method = "NJ",
                              filter_trees = TRUE, fixed_tree = "NJ",
                              outdir, iqtree_bin, threads,
                              n_cores = 1, seed = NULL, timeout = Inf) {
  sim_outdir <- file.path(outdir, "mast_sim_fits")
  dir.create(sim_outdir, showWarnings = FALSE, recursive = TRUE)

  run_one <- function(b) {
    rep_dir <- file.path(sim_outdir, paste0("sim", pad_int(b)))
    dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)

    cand <- mast_candidates(
      alignment     = sim_files[b],
      K_values      = K_values,
      base_model    = base_model,
      rate_model    = rate_model,
      unlinked      = unlinked,
      window_method = window_method,
      filter_trees  = filter_trees,
      outdir        = rep_dir,
      iqtree_bin    = iqtree_bin,
      threads       = threads,
      timeout       = timeout,
      seed          = if (is.null(seed)) NULL else seed + 1000L * b
    )

    tbl <- fit_mast_all_K(
      alignment    = sim_files[b],
      K_values     = cand$K_values_effective,
      base_model   = base_model,
      rate_model   = cand$rate_model,
      tree_files   = cand$tree_files,
      unlinked     = unlinked,
      fixed_tree   = fixed_tree,
      outdir       = rep_dir,
      label_prefix = "fit_",
      iqtree_bin   = iqtree_bin,
      threads      = threads,
      timeout      = timeout,
      mast_max_fit = cand$mast_max
    )
    tbl$replicate  <- b
    tbl$rate_model <- cand$rate_model
    # Per-replicate audit trail: how far the candidate-tree filter cut this
    # replicate back, and whether that alone put K_best out of its reach.
    tbl$K_max_effective    <- cand$K_max_effective
    tbl$n_dropped          <- sum(cand$drop_counts)
    tbl$drop_failed        <- cand$drop_counts[["failed"]]
    tbl$drop_unparsable    <- cand$drop_counts[["unparsable"]]
    tbl$drop_taxon_set     <- cand$drop_counts[["taxon_set"]]
    tbl$drop_uninformative <- cand$drop_counts[["uninformative"]]
    tbl$n_distinct         <- cand$n_distinct
    tbl$below_K_best       <- cand$K_max_effective < K_best
    tbl
  }
  run_one_safe <- function(b) tryCatch(run_one(b), error = function(e) e)

  sim_results <- run_parallel(
    seq_along(sim_files), run_one_safe, n_cores = n_cores
  )

  # Drop failed replicates (with a warning) rather than aborting the whole
  # family; only error if nothing survives. See assess_power() for rationale.
  failed <- vapply(sim_results, inherits, logical(1), "error")
  if (any(failed)) {
    idx <- which(failed)
    warning(
      "MAST bootstrap refit(s) failed for ", sum(failed), "/", length(failed),
      " replicate(s) (dropped): ", paste(idx, collapse = ", "),
      "\nFirst error: ", conditionMessage(sim_results[[idx[1]]])
    )
  }
  sim_results <- sim_results[!failed]
  if (length(sim_results) == 0) {
    stop("All MAST bootstrap replicates failed; cannot compute power.")
  }

  sim_ic <- do.call(rbind, sim_results)

  best_per_rep <- best_K_per_replicate(sim_ic, ic)
  power <- mean(best_per_rep == K_best, na.rm = TRUE)

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
#' @param K_values Unused, kept for call compatibility: K now comes from the
#'   `K` column of each table, since empirical and replicate K ranges can
#'   differ from the requested range.
#' @return Named list with elements AIC, AICc, BIC, each containing
#'   `K_best` (integer) and `power` (numeric between 0 and 1).
compute_power_all_ic <- function(empirical_ic, sim_ic, K_values = NULL) {
  ics <- c("AIC", "AICc", "BIC")
  result <- lapply(stats::setNames(ics, ics), function(crit) {
    K_best <- empirical_ic$K[which.min(empirical_ic[[crit]])]
    best_per_rep <- best_K_per_replicate(sim_ic, crit)
    list(K_best = K_best,
         power  = mean(best_per_rep == K_best, na.rm = TRUE))
  })
  result
}
