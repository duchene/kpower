# MAST tree-mixture workflow: window splitting, tree estimation, ranking

#' Split an alignment into K non-overlapping windows of near-equal size
#'
#' @param alignment Path to the input alignment.
#' @param K Number of windows.
#' @param outdir Output directory.
#' @return Character vector of paths to the K window alignments (PHYLIP format).
split_alignment_windows <- function(alignment, K, outdir) {
  seqs    <- read_alignment(alignment)
  n_sites <- nchar(seqs[1])
  breaks  <- round(seq(0, n_sites, length.out = K + 1))

  win_dir <- file.path(outdir, "windows")
  dir.create(win_dir, showWarnings = FALSE, recursive = TRUE)

  paths <- character(K)
  for (i in seq_len(K)) {
    start <- breaks[i] + 1L
    end   <- breaks[i + 1]
    win   <- subset_sites(seqs, start, end)
    paths[i] <- file.path(win_dir, paste0("window_", pad_int(i), ".phy"))
    write_phylip(win, paths[i])
  }
  paths
}

#' Estimate a tree for each alignment window using ModelFinder
#'
#' Runs `iqtree -m MFP` on each window to select the best model and infer a
#' tree under that model via full heuristic search.
#'
#' @param windows Character vector of window alignment paths.
#' @param outdir Output directory.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param threads Number of threads.
#' @param timeout Per-window timeout in seconds.
#' @return List of per-window results, each with: window, treefile, iqtree_file,
#'   best_model, tree (Newick string).
estimate_window_trees <- function(windows, outdir, iqtree_bin, threads,
                                  timeout) {
  tree_dir <- file.path(outdir, "window_trees")
  dir.create(tree_dir, showWarnings = FALSE, recursive = TRUE)

  lapply(seq_along(windows), function(i) {
    label  <- paste0("window_", pad_int(i))
    prefix <- make_prefix(tree_dir, label)

    args <- c(
      "-s", windows[i],
      "-m", "MFP",
      "--prefix", prefix,
      "-T", threads,
      "--redo"
    )
    run_iqtree(iqtree_bin, args, timeout = timeout)

    iqtree_file <- paste0(prefix, ".iqtree")
    treefile    <- paste0(prefix, ".treefile")
    best_model  <- parse_best_model(iqtree_file)

    list(
      window      = i,
      treefile    = treefile,
      iqtree_file = iqtree_file,
      best_model  = best_model,
      tree        = readLines(treefile)[1]
    )
  })
}

#' Parse the best-fit model from an IQ-TREE report (after MFP)
#'
#' @param iqtree_file Path to a `.iqtree` report file.
#' @return Best-fit model string, or NA if not found.
parse_best_model <- function(iqtree_file) {
  lines <- readLines(iqtree_file)
  hit <- grep("^Best-fit model according to", lines, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  sub(".*:\\s+", "", hit[1])
}

#' Determine within-class rate heterogeneity from window MFP results
#'
#' Examines the best-fit models from each window and summarises the rate
#' heterogeneity component.  Preference order: +R > +G.  The number of
#' categories is the most common value across windows.  Falls back to `+R4`
#' when parsing fails.
#'
#' @param window_results List of per-window results from
#'   `estimate_window_trees()`.
#' @return Character string like `"+R3"` or `"+G4"`.
determine_rate_heterogeneity <- function(window_results) {
  models <- vapply(window_results, function(r) r$best_model, character(1))

  # Try to extract +R{n} from each model
  r_hits <- regmatches(models, regexpr("\\+R[0-9]+", models))
  if (length(r_hits) > 0 && !all(is.na(r_hits))) {
    r_nums <- as.integer(sub("\\+R", "", r_hits))
    r_nums <- r_nums[!is.na(r_nums)]
    if (length(r_nums) > 0) {
      # Mode (most common); ties broken by taking the first
      tbl <- table(r_nums)
      n   <- as.integer(names(tbl)[which.max(tbl)])
      return(paste0("+R", n))
    }
  }

  # Fall back to +G{n}
  g_hits <- regmatches(models, regexpr("\\+G[0-9]*", models))
  if (length(g_hits) > 0 && !all(is.na(g_hits))) {
    g_nums <- as.integer(sub("\\+G", "", g_hits))
    g_nums[is.na(g_nums)] <- 4L  # IQ-TREE default
    tbl <- table(g_nums)
    n   <- as.integer(names(tbl)[which.max(tbl)])
    return(paste0("+G", n))
  }

  "+R4"
}

#' Write all candidate trees (one per line) to a single Newick file
#'
#' @param window_results List of per-window results.
#' @param outfile Path for the combined tree file.
#' @return `outfile` (invisibly).
collect_candidate_trees <- function(window_results, outfile) {
  trees <- vapply(window_results, function(r) r$tree, character(1))
  writeLines(trees, outfile)
  invisible(outfile)
}

#' Rank tree indices by descending weight
#'
#' @param tree_weights Numeric vector of MAST tree weights.
#' @return Integer vector of tree indices sorted by descending weight.
rank_trees_by_weight <- function(tree_weights) {
  order(tree_weights, decreasing = TRUE)
}

#' Build tree files for each K in K_min:K_max using ranked trees
#'
#' For K = 1 no tree file is needed (uses BioNJ). For K >= 2, writes the
#' top-K trees (by weight from the K_max run) to a Newick file.
#'
#' @param all_trees Character vector of Newick tree strings (one per candidate).
#' @param ranked_indices Integer vector of tree indices in weight-descending
#'   order.
#' @param K_values Integer vector of K values to prepare.
#' @param outdir Output directory.
#' @return Named list mapping K (as character) to tree file paths. K = 1 maps
#'   to `NULL`.
build_mast_tree_files <- function(all_trees, ranked_indices, K_values, outdir) {
  tree_dir <- file.path(outdir, "mast_tree_sets")
  dir.create(tree_dir, showWarnings = FALSE, recursive = TRUE)

  files <- list()
  for (K in K_values) {
    if (K <= 1) {
      files[["1"]] <- NULL
      next
    }
    tf <- file.path(tree_dir, paste0("trees_K", K, ".newick"))
    writeLines(all_trees[ranked_indices[seq_len(K)]], tf)
    files[[as.character(K)]] <- tf
  }
  files
}

#' Fit a MAST model to an alignment with a given set of candidate trees
#'
#' @param alignment Path to the input alignment.
#' @param tree_file Path to a Newick file with one tree per line.
#' @param model_str Full IQ-TREE model string including `+T`
#'   (e.g. `"GTR+FO+R3+T"`).
#' @param outdir Output directory.
#' @param label Label for the `--prefix`.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param threads Number of threads.
#' @param timeout Timeout in seconds.
#' @return Named list: model_string, prefix, treefile, iqtree_file, lnL, df,
#'   AIC, AICc, BIC, tree_weights.
fit_mast_model <- function(alignment, tree_file, model_str, outdir, label,
                           iqtree_bin, threads, timeout) {
  prefix <- make_prefix(outdir, label)

  args <- c(
    "-s", alignment,
    "-m", model_str,
    "-te", tree_file,
    "--prefix", prefix,
    "-T", threads,
    "-wspm",
    "--redo"
  )

  run_iqtree(iqtree_bin, args, timeout = timeout)

  iqtree_file <- paste0(prefix, ".iqtree")
  report      <- parse_iqtree_report(iqtree_file)
  weights     <- parse_tree_weights(iqtree_file)
  treefile    <- paste0(prefix, ".treefile")

  c(
    list(model_string  = model_str,
         prefix        = prefix,
         treefile      = treefile,
         iqtree_file   = iqtree_file,
         tree_weights  = weights),
    report
  )
}

#' Fit all K values for MAST analysis
#'
#' K = 1 is a standard single-tree BioNJ fit. K >= 2 uses MAST with the
#' top-K trees.
#'
#' @param alignment Path to the input alignment.
#' @param K_values Integer vector of K values.
#' @param base_model Base substitution model string.
#' @param rate_model Rate heterogeneity string (e.g. `"+R3"`).
#' @param tree_files Named list mapping K (character) to tree file paths,
#'   as returned by `build_mast_tree_files()`.
#' @param unlinked Logical; if TRUE, use MIX syntax for unlinked per-tree
#'   substitution parameters (*T mode).
#' @param outdir Output directory.
#' @param label_prefix Prefix for per-K labels.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param threads Number of threads.
#' @param timeout Per-run timeout in seconds.
#' @return Data frame with columns: K, lnL, df, AIC, AICc, BIC.
fit_mast_all_K <- function(alignment, K_values, base_model, rate_model,
                           tree_files, unlinked = FALSE, outdir,
                           label_prefix = "",
                           iqtree_bin, threads, timeout) {
  results <- lapply(K_values, function(K) {
    label <- paste0(label_prefix, "K", K)

    if (K == 1) {
      # Standard single-tree fit with BioNJ
      model_str <- build_mast_model_str(base_model, rate_model, 1, unlinked)
      fit <- fit_model(
        alignment  = alignment,
        K          = 1,
        base_model = model_str,
        mix_type   = "+R",       # dummy — won't be appended for K = 1
        fixed_tree = "NJ",
        outdir     = outdir,
        label      = label,
        iqtree_bin = iqtree_bin,
        threads    = threads,
        timeout    = timeout
      )
    } else {
      # MAST fit with top-K trees
      model_str <- build_mast_model_str(base_model, rate_model, K, unlinked)
      tf <- tree_files[[as.character(K)]]
      fit <- fit_mast_model(
        alignment  = alignment,
        tree_file  = tf,
        model_str  = model_str,
        outdir     = outdir,
        label      = label,
        iqtree_bin = iqtree_bin,
        threads    = threads,
        timeout    = timeout
      )
    }

    data.frame(
      K    = K,
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
