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
#' @param mix_type Mixture family: `"+R"` (FreeRate, default), `"*R"` (FreeRate
#'   unlinked), `"+H"` (GHOST linked), `"*H"` (GHOST unlinked), `"+T"` (MAST
#'   tree mixtures linked), or `"*T"` (MAST tree mixtures unlinked).
#' @param ic Information criterion used to select K_best and compute power:
#'   `"AIC"`, `"AICc"`, or `"BIC"` (default `"BIC"`).
#' @param fixed_tree Tree handling: `"NJ"` (per-K BioNJ, default), a path to
#'   a fixed tree file, or `NULL` (heuristic search with `--fast`). Ignored when
#'   `mix_type = "+T"` (MAST derives trees from alignment windows).
#' @param rate_model Within-class rate heterogeneity for `+T`/`*T`
#'   (e.g. `"+R4"`). When `NULL` (default) it is derived from the per-window
#'   model selections of whichever alignment is being analysed, so each
#'   bootstrap replicate derives its own.
#' @param fast_trees Deprecated alias for `tree_search = "fast"`.
#' @param filter_trees Logical, `+T`/`*T` only; default `TRUE`. Drop candidate
#'   trees MAST cannot use — a window that produced no tree, a tree on a
#'   gap-trimmed subset of the taxa, or a star tree with no bipartition — and
#'   reduce `K_max` to the number of survivors, reporting it as
#'   `K_max_effective`. Windows are not re-split and trees are not
#'   re-inferred. Duplicate topologies are kept (they are a real
#'   per-tree-branch-length model comparison, which the IC judges on its
#'   own) and only counted. `FALSE` restores the previous behaviour, where an
#'   unusable window aborts the analysis.
#' @param tree_search Tree estimation for `+T`/`*T`: `"NJ"` (BioNJ, default)
#'   or `"fast"` (`--fast` heuristic search). Applies to both the window trees
#'   and the K = 1 tree, so the two cannot disagree, and is applied identically
#'   to the empirical alignment and to every replicate. Overrides `fixed_tree`
#'   on the MAST pathway.
#' @param B Integer number of parametric bootstrap replicates (default 1000).
#' @param seed Integer random seed passed to AliSim (default 1).
#' @param outdir Directory for all IQ-TREE output files. Defaults to a
#'   temporary directory.
#' @param iqtree_bin Path to the IQ-TREE executable. Detected automatically
#'   if not supplied.
#' @param n_cores Number of parallel R workers for bootstrap refits via
#'   `future.apply::future_lapply()` (default 1).
#' @param threads Number of threads for each IQ-TREE run (`-T`). Default 1.
#'   Independent of `n_cores`, so e.g. `n_cores = 4, threads = 2` runs 4
#'   parallel refits each using 2 IQ-TREE threads (8 CPUs total).
#' @param timeout Per-run IQ-TREE timeout in seconds. Default `10 * 3600`
#'   (10 hours): a safety net so a hung fit cannot stall a run indefinitely.
#'   Set to `Inf` to disable, or lower it for shorter jobs.
#'
#' @return An object of class `kpower_result`, a list containing:
#'   \describe{
#'     \item{empirical}{Data frame of K, lnL, df, AIC, AICc, BIC for the
#'       empirical alignment.}
#'     \item{sim_ic}{Long-format data frame of replicate, K, and IC scores
#'       from all bootstrap fits.}
#'     \item{K_best}{Integer; K selected from empirical data under `ic`.}
#'     \item{power}{Numeric; proportion of simulations that recover K_best
#'       under `ic`.}
#'     \item{power_all}{Named list with elements AIC, AICc, BIC, each a list
#'       of K_best and power under that criterion.}
#'     \item{ic}{Character; the primary IC used.}
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
                   rate_model = NULL,
                   tree_search = "NJ",
                   fast_trees = FALSE,
                   filter_trees = TRUE,
                   B          = 1000L,
                   seed       = 1L,
                   outdir     = tempdir(),
                   iqtree_bin = find_iqtree(),
                   n_cores    = 1L,
                   threads    = 1L,
                   timeout    = 10 * 3600) {

  ic      <- match.arg(ic, c("AIC", "AICc", "BIC"))
  threads <- as.character(threads)
  if (isTRUE(fast_trees)) tree_search <- "fast"   # legacy alias

  # ---- MAST pathway (+T or *T) ----------------------------------------------
  if (mix_type %in% c("+T", "*T")) {
    return(kpower_mast(
      alignment  = alignment,
      K_max      = K_max,
      K_min      = K_min,
      base_model = base_model,
      mix_type   = mix_type,
      ic         = ic,
      fixed_tree = fixed_tree,
      rate_model = rate_model,
      tree_search = tree_search,
      filter_trees = filter_trees,
      B          = B,
      seed       = seed,
      outdir     = outdir,
      iqtree_bin = iqtree_bin,
      n_cores    = n_cores,
      threads    = threads,
      timeout    = timeout
    ))
  }

  # ---- +R / *R / +H / *H pathway ------------------------------------------
  K_values <- seq.int(K_min, K_max)

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
  K_best     <- empirical_ic$K[which.min(empirical_ic[[ic]])]
  best_label  <- paste0("empirical_K", K_best)
  best_prefix <- make_prefix(emp_outdir, best_label)
  best_fit <- list(
    K            = K_best,
    model_string = build_model_str(base_model, mix_type, K_best),
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

  # --- Step 5: Compute power under all ICs -----------------------------------
  power_all <- compute_power_all_ic(
    empirical_ic = empirical_ic,
    sim_ic       = power_result$sim_ic,
    K_values     = K_values
  )

  message(sprintf(
    "Power: %.1f%% of simulations recover K_best = %d under %s",
    power_all[[ic]]$power * 100, K_best, ic
  ))

  # --- Step 6: Build figure -------------------------------------------------
  fig <- plot_kpower(
    empirical_ic = empirical_ic,
    sim_ic       = power_result$sim_ic,
    K_best       = K_best,
    power        = power_all[[ic]]$power,
    ic           = ic
  )

  structure(
    list(
      empirical = empirical_ic,
      sim_ic    = power_result$sim_ic,
      K_best    = K_best,
      power     = power_all[[ic]]$power,
      power_all = power_all,
      ic        = ic,
      mix_type  = mix_type,
      plot      = fig
    ),
    class = "kpower_result"
  )
}

#' Print method for kpower_result
#' @param x A `kpower_result` object.
#' @param ... Ignored.
#' @export
print.kpower_result <- function(x, ...) {
  cat("kpower result\n")
  cat("  mix_type :", x$mix_type %||% "n/a", "\n")
  cat("  IC used  :", x$ic, "\n")
  cat("  K_best   :", x$K_best, "\n")
  cat(sprintf("  Power    : %.1f%%\n", x$power * 100))
  if (!is.null(x$rate_model))
    cat("  Rate het :", x$rate_model, "\n")
  if (!is.null(x$power_all)) {
    cat("  ---\n")
    for (crit in c("AIC", "AICc", "BIC")) {
      pa <- x$power_all[[crit]]
      cat(sprintf("  %-4s : K_best = %d, Power = %.1f%%\n",
                  crit, pa$K_best, pa$power * 100))
    }
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# MAST (tree-mixture) pathway
# ---------------------------------------------------------------------------

#' Power assessment for MAST tree-mixture models
#'
#' Internal workhorse called by `kpower(mix_type = "+T")`.
#'
#' Workflow:
#' 1. Split alignment into K_max windows.
#' 2. Estimate a candidate tree per window (BioNJ or `--fast`, per
#'    `tree_search`).
#' 3. Determine within-class rate heterogeneity from window model selections.
#' 3b. Drop candidate trees MAST cannot use and reduce K_max to the survivors
#'    (`filter_candidate_trees()`; skipped when `filter_trees = FALSE`).
#' 4. Fit MAST with all surviving candidate trees.
#' 5. Rank trees by weight; build tree sets for K = 1..K_max_effective.
#' 6. Fit all K values (K = 1 is a standard BioNJ single-tree fit).
#' 7. Select K_best by IC.
#' 8. Simulate B alignments by per-class AliSim + concatenation.
#' 9. Refit all K on each replicate.
#' 10. Report power and IC profile figure.
#'
#' @inheritParams kpower
#' @param mix_type Either `"+T"` (linked) or `"*T"` (unlinked per-tree
#'   substitution parameters via `MIX{...}+T`).
#' @return Object of class `kpower_result`.
#' @keywords internal
kpower_mast <- function(alignment, K_max, K_min = 1L, base_model = "GTR",
                        mix_type = "+T", ic = "BIC", fixed_tree = "NJ",
                        rate_model = NULL, tree_search = "NJ",
                        filter_trees = TRUE,
                        B = 1000L, seed = 1L,
                        outdir = tempdir(), iqtree_bin = find_iqtree(),
                        n_cores = 1L, threads = "1", timeout = 10 * 3600) {

  K_values <- seq.int(K_min, K_max)
  unlinked <- (mix_type == "*T")
  n_sites  <- alignment_length(alignment)

  # Window trees and the K = 1 tree must use the same estimator.
  ts            <- resolve_tree_search(tree_search)
  window_method <- ts$window_method
  fixed_tree    <- ts$fixed_tree

  # --- Steps 1-6: candidate trees, ranking, nested tree sets ---------------
  # Same function is run on every bootstrap replicate, so no empirical tree
  # can leak into a refit.
  message("Deriving candidate trees (", window_method, ") from ", K_max,
          " windows ...")
  emp <- mast_candidates(
    alignment     = alignment,
    K_values      = K_values,
    base_model    = base_model,
    rate_model    = rate_model,
    unlinked      = unlinked,
    window_method = window_method,
    filter_trees  = filter_trees,
    outdir        = file.path(outdir, "empirical"),
    iqtree_bin    = iqtree_bin,
    threads       = threads,
    timeout       = timeout,
    seed          = seed
  )
  rate_model   <- emp$rate_model
  all_trees    <- emp$all_trees
  ranked       <- emp$ranked
  tree_files   <- emp$tree_files
  mast_max     <- emp$mast_max
  mast_max_dir <- emp$mast_dir

  # The filter may have cut K_max back; everything downstream uses the
  # effective range, and K_max_effective is reported so a K_best sitting on a
  # reduced ceiling is detectable rather than silent.
  K_values        <- emp$K_values_effective
  K_max_effective <- emp$K_max_effective

  message("Rate heterogeneity: ", rate_model)
  message("Tree ranking by weight: ", paste(ranked, collapse = ", "))

  message("Fitting K = ", min(K_values), " to ", max(K_values),
          " on empirical alignment ...")
  empirical_ic <- fit_mast_all_K(
    alignment    = alignment,
    K_values     = K_values,
    base_model   = base_model,
    rate_model   = rate_model,
    tree_files   = tree_files,
    unlinked     = unlinked,
    fixed_tree   = fixed_tree,
    outdir       = mast_max_dir,
    label_prefix = "empirical_",
    iqtree_bin   = iqtree_bin,
    threads      = threads,
    timeout      = timeout,
    mast_max_fit = mast_max
  )

  # --- Step 7: Select K_best ------------------------------------------------
  K_best <- empirical_ic$K[which.min(empirical_ic[[ic]])]
  message("K_best = ", K_best, " (selected by ", ic, ")")
  if (K_best == K_max_effective && K_max_effective < K_max)
    message("Note: K_best sits on the reduced ceiling (K_max_effective = ",
            K_max_effective, " < K_max = ", K_max,
            "), so it reads as \"at least ", K_best, "\".")

  # Retrieve the K_best fit result for simulation
  if (K_best == 1) {
    # K=1 has no MAST tree weights; skip MAST simulation
    # Fall back to standard single-tree simulation
    best_label  <- paste0("empirical_K1")
    best_prefix <- make_prefix(mast_max_dir, best_label)
    best_fit <- list(
      K            = 1,
      model_string = paste0(base_model, "+FO", rate_model),
      treefile     = paste0(best_prefix, ".treefile"),
      iqtree_file  = paste0(best_prefix, ".iqtree")
    )
    use_mast_sim <- FALSE
  } else {
    # Re-use the MAST fit for K_best (or the effective K_max if it is that)
    if (K_best == K_max_effective) {
      best_mast <- mast_max
    } else {
      message("Refitting MAST with K_best = ", K_best, " trees ...")
      kbest_model_str <- build_mast_model_str(base_model, rate_model, K_best,
                                              unlinked)
      best_mast <- fit_mast_model(
        alignment  = alignment,
        tree_file  = tree_files[[as.character(K_best)]],
        model_str  = kbest_model_str,
        outdir     = mast_max_dir,
        label      = paste0("mast_Kbest_", K_best),
        iqtree_bin = iqtree_bin,
        threads    = threads,
        timeout    = timeout
      )
    }
    use_mast_sim <- TRUE
  }

  # --- Step 8: Simulate B alignments ----------------------------------------
  message("Simulating ", B, " alignments of ", n_sites, " sites ...")
  if (use_mast_sim) {
    sim_files <- simulate_mast_alignments(
      mast_result = best_mast,
      base_model  = base_model,
      n_sites     = n_sites,
      B           = B,
      seed        = seed,
      outdir      = outdir,
      iqtree_bin  = iqtree_bin,
      threads     = threads,
      timeout     = timeout,
      alignment   = alignment
    )
  } else {
    # K_best = 1: standard single-tree AliSim
    sim_files <- simulate_alignments(
      fit_result = best_fit,
      alignment  = alignment,
      n_sites    = n_sites,
      B          = B,
      outdir     = outdir,
      seed       = seed,
      iqtree_bin = iqtree_bin,
      threads    = threads,
      timeout    = timeout
    )
  }

  # --- Step 9: Refit all K on each replicate --------------------------------
  message("Refitting K = ", K_min, " to ", K_max,
          " on ", B, " simulated alignments ...")
  # Replicates are given the *requested* K range, not the empirical
  # alignment's effective one: each derives its own candidate trees and so its
  # own effective ceiling. No floor is applied at the empirical K_best.
  power_result <- assess_mast_power(
    sim_files  = sim_files,
    K_values   = seq.int(K_min, K_max),
    K_best     = K_best,
    ic         = ic,
    base_model = base_model,
    rate_model = rate_model,
    unlinked   = unlinked,
    window_method = window_method,
    filter_trees = filter_trees,
    fixed_tree = fixed_tree,
    outdir     = outdir,
    iqtree_bin = iqtree_bin,
    threads    = threads,
    n_cores    = n_cores,
    seed       = seed,
    timeout    = timeout
  )

  # --- Step 10: Power under all ICs and figure ------------------------------
  power_all <- compute_power_all_ic(
    empirical_ic = empirical_ic,
    sim_ic       = power_result$sim_ic,
    K_values     = K_values
  )

  message(sprintf(
    "Power: %.1f%% of simulations recover K_best = %d under %s",
    power_all[[ic]]$power * 100, K_best, ic
  ))

  fig <- plot_kpower(
    empirical_ic = empirical_ic,
    sim_ic       = power_result$sim_ic,
    K_best       = K_best,
    power        = power_all[[ic]]$power,
    ic           = ic
  )

  structure(
    list(
      empirical    = empirical_ic,
      sim_ic       = power_result$sim_ic,
      K_best       = K_best,
      power        = power_all[[ic]]$power,
      power_all    = power_all,
      ic           = ic,
      mix_type     = mix_type,
      rate_model   = rate_model,
      tree_weights = if (use_mast_sim) best_mast$tree_weights else NULL,
      ranked_trees = ranked,
      K_max_requested = K_max,
      K_max_effective = K_max_effective,
      ceiling_limited = (K_best == K_max_effective && K_max_effective < K_max),
      drop_counts     = emp$drop_counts,
      drop_reasons    = emp$drop_reasons,
      n_distinct      = emp$n_distinct,
      plot         = fig
    ),
    class = "kpower_result"
  )
}
