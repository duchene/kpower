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
#'   a fixed tree file, or `NULL` (full heuristic search). Ignored when
#'   `mix_type = "+T"` (MAST derives trees from alignment windows).
#' @param B Integer number of parametric bootstrap replicates (default 1000).
#' @param seed Integer random seed passed to AliSim (default 1).
#' @param outdir Directory for all IQ-TREE output files. Defaults to a
#'   temporary directory.
#' @param iqtree_bin Path to the IQ-TREE executable. Detected automatically
#'   if not supplied.
#' @param n_cores Number of parallel R workers for bootstrap refits via
#'   `parallel::mclapply()` (default 1).
#' @param threads Number of threads for each IQ-TREE run (`-T`). Default 1.
#'   Independent of `n_cores`, so e.g. `n_cores = 4, threads = 2` runs 4
#'   parallel refits each using 2 IQ-TREE threads (8 CPUs total).
#' @param timeout Per-run timeout in seconds (default 3600).
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
                   fast_trees = FALSE,
                   B          = 1000L,
                   seed       = 1L,
                   outdir     = tempdir(),
                   iqtree_bin = find_iqtree(),
                   n_cores    = 1L,
                   threads    = 1L,
                   timeout    = 3600L) {

  ic      <- match.arg(ic, c("AIC", "AICc", "BIC"))
  threads <- as.character(threads)

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
      fast_trees = fast_trees,
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
  K_best     <- K_values[which.min(empirical_ic[[ic]])]
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
#' 2. Run MFP on each window to estimate a candidate tree.
#' 3. Determine within-class rate heterogeneity from window model selections.
#' 4. Fit MAST with all K_max candidate trees.
#' 5. Rank trees by weight; build tree sets for K = 1..K_max.
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
                        fast_trees = FALSE,
                        B = 1000L, seed = 1L,
                        outdir = tempdir(), iqtree_bin = find_iqtree(),
                        n_cores = 1L, threads = "1", timeout = 3600L) {

  K_values <- seq.int(K_min, K_max)
  unlinked <- (mix_type == "*T")
  n_sites  <- alignment_length(alignment)

  # --- Step 1: Split alignment into K_max windows ---------------------------
  message("Splitting alignment into ", K_max, " windows ...")
  windows <- split_alignment_windows(alignment, K_max, outdir)

  # --- Step 2: Estimate tree per window (MFP + full search) -----------------
  message("Estimating trees for ", K_max, " windows via ModelFinder ...")
  window_results <- estimate_window_trees(
    windows, outdir, iqtree_bin, threads, timeout,
    fast_trees = fast_trees
  )

  # --- Step 3: Determine rate heterogeneity ---------------------------------
  rate_model <- determine_rate_heterogeneity(window_results)
  message("Rate heterogeneity from windows: ", rate_model)

  # --- Step 4: Collect candidate trees and run MAST with K_max --------------
  tree_file <- collect_candidate_trees(
    window_results, file.path(outdir, "candidate_trees.newick")
  )
  all_trees      <- readLines(tree_file)
  mast_model_str <- build_mast_model_str(base_model, rate_model, K_max,
                                         unlinked)

  message("Fitting MAST (", mast_model_str, ") with ", K_max, " trees ...")
  mast_max_dir <- file.path(outdir, "mast_empirical")
  dir.create(mast_max_dir, showWarnings = FALSE, recursive = TRUE)

  mast_max <- fit_mast_model(
    alignment  = alignment,
    tree_file  = tree_file,
    model_str  = mast_model_str,
    outdir     = mast_max_dir,
    label      = paste0("mast_K", K_max),
    iqtree_bin = iqtree_bin,
    threads    = threads,
    timeout    = timeout
  )

  # --- Step 5: Rank trees by weight -----------------------------------------
  ranked <- rank_trees_by_weight(mast_max$tree_weights)
  message("Tree ranking by weight: ", paste(ranked, collapse = ", "))

  # --- Step 6: Build tree files for each K and fit all K --------------------
  tree_files <- build_mast_tree_files(all_trees, ranked, K_values, outdir)

  message("Fitting K = ", K_min, " to ", K_max, " on empirical alignment ...")
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
    timeout      = timeout
  )

  # --- Step 7: Select K_best ------------------------------------------------
  K_best <- K_values[which.min(empirical_ic[[ic]])]
  message("K_best = ", K_best, " (selected by ", ic, ")")

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
    # Re-use the MAST fit for K_best (or K_max if K_best == K_max)
    if (K_best == K_max) {
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
      timeout     = timeout
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
  power_result <- assess_mast_power(
    sim_files  = sim_files,
    K_values   = K_values,
    K_best     = K_best,
    ic         = ic,
    base_model = base_model,
    rate_model = rate_model,
    tree_files = tree_files,
    unlinked   = unlinked,
    fixed_tree = fixed_tree,
    outdir     = outdir,
    iqtree_bin = iqtree_bin,
    threads    = threads,
    n_cores    = n_cores,
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
      plot         = fig
    ),
    class = "kpower_result"
  )
}
