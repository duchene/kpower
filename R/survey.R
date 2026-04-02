# Cross-model survey: two-phase pipeline
#
# Phase 1: Fit empirical IC profiles for every requested mixture family.
# Phase 2: For each family's K_best, run the parametric bootstrap to
#           estimate power.  The two phases share infrastructure so the
#           bootstrap is not a separate "optional" step but the natural
#           second half of the survey.

#' Survey model complexity across mixture families with power assessment
#'
#' **Phase 1 — Empirical fitting.** For each mixture family in `mix_types`,
#' fit K = K_min..K_max to the empirical alignment and record IC profiles.
#' Identify K_best per family under the chosen IC.
#'
#' **Phase 2 — Parametric bootstrap.** For each family, simulate B
#' alignments from the K_best model, refit the full K range on each
#' replicate, and compute power (proportion recovering K_best).
#'
#' The result is a comparative table of all families ranked by IC at their
#' K_best, each annotated with the bootstrap power estimate.  This directly
#' answers: *which model family and complexity does the alignment best
#' support, and how reliably?*
#'
#' @param alignment Path to the input alignment (FASTA or PHYLIP).
#' @param K_max Integer; maximum number of mixture categories.
#' @param K_min Integer; minimum K (default 1).
#' @param base_model Base substitution model (default `"GTR"`).
#' @param mix_types Character vector of mixture families to survey.
#'   Default: `c("+R", "*R", "+H", "*H", "+T", "*T")`.
#' @param ic Primary information criterion: `"AIC"`, `"AICc"`, or `"BIC"`
#'   (default `"BIC"`).
#' @param fixed_tree Tree handling for rate/GHOST models (default `"NJ"`).
#'   Ignored for `+T`/`*T`.
#' @param B Integer number of parametric bootstrap replicates per family
#'   (default 100).
#' @param seed Integer random seed (default 1).
#' @param outdir Output directory.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param n_cores Number of parallel R workers for bootstrap refits.
#' @param threads Number of threads per IQ-TREE run.
#' @param timeout Per-run timeout in seconds.
#'
#' @return An object of class `kpower_survey`, a list containing:
#'   \describe{
#'     \item{families}{Named list; for each mix_type, a list with
#'       `empirical` (IC table), `K_best`, `sim_ic` (bootstrap results),
#'       `power_all` (power under all three ICs).}
#'     \item{comparison}{Data frame ranking families by IC at K_best,
#'       with power under each criterion.}
#'     \item{best}{Named list: `mix_type`, `K`, `ic_value` for the overall
#'       winner under the primary IC.}
#'     \item{ic}{Character; the primary IC.}
#'     \item{plots}{Named list of ggplot2 IC-profile figures, one per
#'       mix_type.}
#'   }
#' @export
kpower_survey <- function(alignment,
                          K_max,
                          K_min      = 1L,
                          base_model = "GTR",
                          mix_types  = c("+R", "*R", "+H", "*H", "+T", "*T"),
                          ic         = "BIC",
                          fixed_tree = "NJ",
                          fast_trees = FALSE,
                          B          = 100L,
                          seed       = 1L,
                          outdir     = tempdir(),
                          iqtree_bin = find_iqtree(),
                          n_cores    = 1L,
                          threads    = 1L,
                          timeout    = 3600L) {

  ic       <- match.arg(ic, c("AIC", "AICc", "BIC"))
  K_values <- seq.int(K_min, K_max)
  threads  <- as.character(threads)
  n_sites  <- alignment_length(alignment)

  families <- list()
  plots    <- list()

  for (mt in mix_types) {
    mt_dir <- file.path(outdir, survey_dir_label(mt))
    dir.create(mt_dir, showWarnings = FALSE, recursive = TRUE)

    message("\n===== ", mt, " =====")

    family <- tryCatch(
      survey_one_family(
        alignment  = alignment,
        K_values   = K_values,
        base_model = base_model,
        mix_type   = mt,
        ic         = ic,
        fixed_tree = fixed_tree,
        fast_trees = fast_trees,
        n_sites    = n_sites,
        B          = B,
        seed       = seed,
        outdir     = mt_dir,
        iqtree_bin = iqtree_bin,
        n_cores    = n_cores,
        threads    = threads,
        timeout    = timeout
      ),
      error = function(e) {
        warning("Survey failed for ", mt, ": ", conditionMessage(e))
        NULL
      }
    )

    families[[mt]] <- family

    if (!is.null(family)) {
      pa <- family$power_all[[ic]]
      plots[[mt]] <- plot_kpower(
        empirical_ic = family$empirical,
        sim_ic       = family$sim_ic,
        K_best       = family$K_best,
        power        = pa$power,
        ic           = ic
      )
    }
  }

  # --- Comparison table and best model --------------------------------------
  comparison <- build_survey_comparison(families, ic)
  best       <- identify_survey_best(comparison, ic)

  message("\n===== Survey complete =====")
  if (!is.na(best$mix_type)) {
    message("Best model: ", best$mix_type, " K=", best$K,
            " (", ic, " = ", round(best$ic_value, 2), ")")
  }

  structure(
    list(
      families   = families,
      comparison = comparison,
      best       = best,
      ic         = ic,
      plots      = plots
    ),
    class = "kpower_survey"
  )
}


# ---------------------------------------------------------------------------
# Internal: run both phases for a single mixture family
# ---------------------------------------------------------------------------

#' Phase 1 + Phase 2 for one mixture family
#'
#' @keywords internal
survey_one_family <- function(alignment, K_values, base_model, mix_type,
                              ic, fixed_tree, fast_trees = FALSE,
                              n_sites, B, seed, outdir,
                              iqtree_bin, n_cores, threads, timeout) {

  is_mast  <- mix_type %in% c("+T", "*T")
  unlinked <- mix_type %in% c("*R", "*T", "*H")

  # ========== Phase 1: Empirical IC profiles ================================

  if (is_mast) {
    phase1 <- survey_phase1_mast(
      alignment  = alignment,
      K_values   = K_values,
      base_model = base_model,
      mix_type   = mix_type,
      ic         = ic,
      fixed_tree = fixed_tree,
      fast_trees = fast_trees,
      outdir     = outdir,
      iqtree_bin = iqtree_bin,
      threads    = threads,
      timeout    = timeout
    )
    empirical_ic <- phase1$empirical_ic
    mast_ctx     <- phase1$mast_ctx
  } else {
    emp_dir <- file.path(outdir, "empirical")
    dir.create(emp_dir, showWarnings = FALSE, recursive = TRUE)

    message("Phase 1: Fitting K = ", min(K_values), "..", max(K_values),
            " (", mix_type, ") ...")
    empirical_ic <- fit_all_K(
      alignment    = alignment,
      K_values     = K_values,
      base_model   = base_model,
      mix_type     = mix_type,
      fixed_tree   = fixed_tree,
      outdir       = emp_dir,
      label_prefix = "empirical_",
      iqtree_bin   = iqtree_bin,
      threads      = threads,
      timeout      = timeout
    )
    mast_ctx <- NULL
  }

  K_best <- K_values[which.min(empirical_ic[[ic]])]
  message("  K_best = ", K_best, " (", ic, ")")

  # ========== Phase 2: Parametric bootstrap =================================

  message("Phase 2: Bootstrap (B = ", B, ") for ", mix_type,
          " K_best = ", K_best, " ...")

  if (is_mast) {
    phase2 <- survey_phase2_mast(
      alignment  = alignment,
      K_values   = K_values,
      K_best     = K_best,
      base_model = base_model,
      mix_type   = mix_type,
      ic         = ic,
      fixed_tree = fixed_tree,
      n_sites    = n_sites,
      B          = B,
      seed       = seed,
      mast_ctx   = mast_ctx,
      outdir     = outdir,
      iqtree_bin = iqtree_bin,
      n_cores    = n_cores,
      threads    = threads,
      timeout    = timeout
    )
  } else {
    phase2 <- survey_phase2_standard(
      alignment  = alignment,
      K_values   = K_values,
      K_best     = K_best,
      empirical_ic = empirical_ic,
      base_model = base_model,
      mix_type   = mix_type,
      fixed_tree = fixed_tree,
      n_sites    = n_sites,
      B          = B,
      seed       = seed,
      outdir     = outdir,
      iqtree_bin = iqtree_bin,
      n_cores    = n_cores,
      threads    = threads,
      timeout    = timeout
    )
  }

  power_all <- compute_power_all_ic(empirical_ic, phase2$sim_ic, K_values)

  message(sprintf("  Power (%s): %.1f%%",
                  ic, power_all[[ic]]$power * 100))

  list(
    empirical  = empirical_ic,
    K_best     = K_best,
    sim_ic     = phase2$sim_ic,
    power_all  = power_all,
    mix_type   = mix_type,
    mast_ctx   = mast_ctx
  )
}

# ---------------------------------------------------------------------------
# Phase 1 for MAST: window trees, ranking, empirical fits
# ---------------------------------------------------------------------------

survey_phase1_mast <- function(alignment, K_values, base_model, mix_type,
                               ic, fixed_tree = "NJ", fast_trees = FALSE,
                               outdir, iqtree_bin,
                               threads, timeout) {
  K_max    <- max(K_values)
  unlinked <- (mix_type == "*T")

  message("Phase 1: MAST windows and empirical fits ...")

  windows <- split_alignment_windows(alignment, K_max, outdir)
  window_results <- estimate_window_trees(
    windows, outdir, iqtree_bin, threads, timeout,
    fast_trees = fast_trees
  )
  rate_model <- determine_rate_heterogeneity(window_results)
  message("  Rate heterogeneity: ", rate_model)

  tree_file <- collect_candidate_trees(
    window_results, file.path(outdir, "candidate_trees.newick")
  )
  all_trees <- readLines(tree_file)

  # Fit MAST with K_max trees to get tree weights and ranking
  mast_max_dir <- file.path(outdir, "mast_empirical")
  dir.create(mast_max_dir, showWarnings = FALSE, recursive = TRUE)

  mast_model_str <- build_mast_model_str(base_model, rate_model, K_max,
                                         unlinked)
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

  ranked     <- rank_trees_by_weight(mast_max$tree_weights)
  tree_files <- build_mast_tree_files(all_trees, ranked, K_values, outdir)

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

  mast_ctx <- list(
    rate_model   = rate_model,
    tree_files   = tree_files,
    all_trees    = all_trees,
    ranked       = ranked,
    mast_max     = mast_max,
    mast_max_dir = mast_max_dir,
    unlinked     = unlinked,
    fixed_tree   = fixed_tree
  )

  list(empirical_ic = empirical_ic, mast_ctx = mast_ctx)
}

# ---------------------------------------------------------------------------
# Phase 2 for standard models: simulate + refit
# ---------------------------------------------------------------------------

survey_phase2_standard <- function(alignment, K_values, K_best,
                                   empirical_ic, base_model, mix_type,
                                   fixed_tree, n_sites, B, seed, outdir,
                                   iqtree_bin, n_cores, threads, timeout) {
  emp_dir     <- file.path(outdir, "empirical")
  best_label  <- paste0("empirical_K", K_best)
  best_prefix <- make_prefix(emp_dir, best_label)

  best_fit <- list(
    K            = K_best,
    model_string = build_model_str(base_model, mix_type, K_best),
    treefile     = paste0(best_prefix, ".treefile"),
    iqtree_file  = paste0(best_prefix, ".iqtree")
  )

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

  assess_power(
    sim_files  = sim_files,
    K_values   = K_values,
    K_best     = K_best,
    ic         = "BIC",     # placeholder — power_all computes all three
    base_model = base_model,
    mix_type   = mix_type,
    fixed_tree = fixed_tree,
    outdir     = outdir,
    iqtree_bin = iqtree_bin,
    threads    = threads,
    n_cores    = n_cores,
    timeout    = timeout
  )
}

# ---------------------------------------------------------------------------
# Phase 2 for MAST: simulate + refit
# ---------------------------------------------------------------------------

survey_phase2_mast <- function(alignment, K_values, K_best, base_model,
                               mix_type, ic, fixed_tree = "NJ", n_sites,
                               B, seed, mast_ctx,
                               outdir, iqtree_bin, n_cores, threads,
                               timeout) {
  unlinked     <- mast_ctx$unlinked
  rate_model   <- mast_ctx$rate_model
  tree_files   <- mast_ctx$tree_files
  mast_max_dir <- mast_ctx$mast_max_dir

  # Get the K_best MAST fit for simulation
  if (K_best == 1) {
    best_label  <- "empirical_K1"
    best_prefix <- make_prefix(mast_max_dir, best_label)
    best_fit <- list(
      K            = 1,
      model_string = build_mast_model_str(base_model, rate_model, 1, unlinked),
      treefile     = paste0(best_prefix, ".treefile"),
      iqtree_file  = paste0(best_prefix, ".iqtree")
    )
    use_mast_sim <- FALSE
  } else {
    K_max <- max(K_values)
    if (K_best == K_max) {
      best_mast <- mast_ctx$mast_max
    } else {
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
    best_fit     <- NULL
  }

  # Simulate
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
  }

  # Refit all K on each replicate
  assess_mast_power(
    sim_files  = sim_files,
    K_values   = K_values,
    K_best     = K_best,
    ic         = "BIC",     # placeholder — power_all computes all three
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
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Label for survey output subdirectory
#' @keywords internal
survey_dir_label <- function(mix_type) {
  base <- gsub("[+*]", "", mix_type)
  link <- if (startsWith(mix_type, "*")) "unlinked" else "linked"
  paste0("survey_", base, "_", link)
}

#' Build comparison data frame from survey family results
#' @keywords internal
build_survey_comparison <- function(families, ic) {
  rows <- lapply(names(families), function(mt) {
    fam <- families[[mt]]
    if (is.null(fam)) {
      return(data.frame(
        mix_type   = mt,
        K_best     = NA_integer_,
        lnL        = NA_real_,
        df         = NA_integer_,
        AIC        = NA_real_,
        AICc       = NA_real_,
        BIC        = NA_real_,
        power_AIC  = NA_real_,
        power_AICc = NA_real_,
        power_BIC  = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    emp      <- fam$empirical
    K_best   <- fam$K_best
    best_row <- emp[emp$K == K_best, ]
    pa       <- fam$power_all

    data.frame(
      mix_type   = mt,
      K_best     = K_best,
      lnL        = best_row$lnL,
      df         = best_row$df,
      AIC        = best_row$AIC,
      AICc       = best_row$AICc,
      BIC        = best_row$BIC,
      power_AIC  = pa$AIC$power,
      power_AICc = pa$AICc$power,
      power_BIC  = pa$BIC$power,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Identify the overall best model from the comparison table
#' @keywords internal
identify_survey_best <- function(comparison, ic) {
  successful <- comparison[!is.na(comparison[[ic]]), ]
  if (nrow(successful) == 0) {
    return(list(mix_type = NA, K = NA, ic_value = NA))
  }
  best_idx <- which.min(successful[[ic]])
  list(
    mix_type = successful$mix_type[best_idx],
    K        = successful$K_best[best_idx],
    ic_value = successful[[ic]][best_idx]
  )
}

#' Print method for kpower_survey
#' @param x A `kpower_survey` object.
#' @param ... Ignored.
#' @export
print.kpower_survey <- function(x, ...) {
  cat("kpower survey\n")
  cat("  Primary IC :", x$ic, "\n")
  if (!is.na(x$best$mix_type)) {
    cat("  Best model :", x$best$mix_type, " K =", x$best$K, "\n")
    cat(sprintf("  Best %-5s : %.2f\n", x$ic, x$best$ic_value))
  }
  cat("\n")

  comp <- x$comparison
  power_col <- paste0("power_", x$ic)
  disp <- data.frame(
    mix_type = comp$mix_type,
    K        = comp$K_best,
    df       = comp$df,
    lnL      = ifelse(is.na(comp$lnL), NA, round(comp$lnL, 1)),
    AIC      = ifelse(is.na(comp$AIC), NA, round(comp$AIC, 1)),
    AICc     = ifelse(is.na(comp$AICc), NA, round(comp$AICc, 1)),
    BIC      = ifelse(is.na(comp$BIC), NA, round(comp$BIC, 1)),
    power    = ifelse(is.na(comp[[power_col]]), "n/a",
                      sprintf("%.1f%%", comp[[power_col]] * 100)),
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)
  invisible(x)
}
