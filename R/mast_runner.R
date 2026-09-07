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

#' Estimate a tree for each alignment window
#'
#' `window_method` selects the estimator:
#' - `"NJ"`: BioNJ topology with `-t BIONJ --tree-fix` under `base_model+R`,
#'   letting IQ-TREE choose the FreeRate category count. Deterministic given
#'   the alignment, and cheap enough to repeat on every bootstrap replicate.
#' - `"fast"`: `--fast` heuristic search under `base_model+R`.
#'
#' Must match the estimator used for the K = 1 fit; see
#' `resolve_tree_search()`, which couples them.
#'
#' @param windows Character vector of window alignment paths.
#' @param outdir Output directory.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param threads Number of threads.
#' @param timeout Per-window timeout in seconds.
#' @param window_method Either `"NJ"` or `"fast"`.
#' @param base_model Base substitution model used by the NJ and fast methods.
#' @param seed Optional base seed; window i uses `seed + i`.
#' @return List of per-window results, each with: window, treefile,
#'   iqtree_file, best_model, tree (Newick string, or `NA` when the window
#'   produced none), ok (logical), error (message or `NA`).
#'
#' A window IQ-TREE cannot analyse — too short, too gappy, non-zero exit, or a
#' clean exit with no tree written — no longer aborts the whole alignment. It
#' comes back with `ok = FALSE`, and `filter_candidate_trees()` drops it.
estimate_window_trees <- function(windows, outdir, iqtree_bin, threads,
                                  timeout, window_method = "NJ",
                                  base_model = "GTR", seed = NULL) {
  window_method <- match.arg(window_method, c("NJ", "fast"))
  tree_dir <- file.path(outdir, "window_trees")
  dir.create(tree_dir, showWarnings = FALSE, recursive = TRUE)

  lapply(seq_along(windows), function(i) {
    label  <- paste0("window_", pad_int(i))
    prefix <- make_prefix(tree_dir, label)

    args <- switch(window_method,
      NJ = c("-s", windows[i], "-m", paste0(base_model, "+R"),
             "-t", "BIONJ", "--tree-fix",
             "--prefix", prefix, "-T", threads, "--redo"),
      fast = c("-s", windows[i], "-m", paste0(base_model, "+R"), "--fast",
               "--prefix", prefix, "-T", threads, "--redo")
    )
    if (!is.null(seed)) args <- c(args, "--seed", as.character(seed + i))

    iqtree_file <- paste0(prefix, ".iqtree")
    treefile    <- paste0(prefix, ".treefile")

    attempt <- tryCatch({
      run_iqtree(iqtree_bin, args, timeout = timeout)
      tree <- read_first_tree(treefile)
      if (is.na(tree))
        stop("IQ-TREE exited cleanly but wrote no usable tree")
      list(tree = tree, error = NA_character_)
    }, error = function(e) {
      list(tree = NA_character_, error = conditionMessage(e))
    })

    if (!is.na(attempt$error))
      message("  Window ", i, " unusable: ", attempt$error)

    list(
      window      = i,
      treefile    = treefile,
      iqtree_file = iqtree_file,
      best_model  = parse_best_model(iqtree_file),
      tree        = attempt$tree,
      ok          = !is.na(attempt$tree),
      error       = attempt$error
    )
  })
}

#' First usable Newick string in a tree file
#'
#' @param treefile Path to a `.treefile`.
#' @return The first non-blank line, or `NA_character_` if the file is
#'   missing, unreadable or empty.
read_first_tree <- function(treefile) {
  if (!file.exists(treefile)) return(NA_character_)
  lines <- tryCatch(readLines(treefile, warn = FALSE),
                    error = function(e) character(0))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(NA_character_)
  lines[1]
}

#' Parse the best-fit model from an IQ-TREE report (after MFP)
#'
#' @param iqtree_file Path to a `.iqtree` report file.
#' @return Best-fit model string, or NA if not found.
parse_best_model <- function(iqtree_file) {
  if (!file.exists(iqtree_file)) return(NA_character_)
  lines <- tryCatch(readLines(iqtree_file, warn = FALSE),
                    error = function(e) character(0))
  # ModelFinder writes "Best-fit model according to ...". Fixed-model runs
  # (e.g. GTR+R with -t BIONJ) do not, but still report the model they used.
  hit <- grep("^Best-fit model according to", lines, value = TRUE)
  if (length(hit) == 0) hit <- grep("^Model of substitution:", lines, value = TRUE)
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
  models <- models[!is.na(models)]
  if (length(models) == 0) return("+R4")

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
#' @param window_results Either a character vector of Newick strings (as
#'   returned by `filter_candidate_trees()`) or a list of per-window results.
#' @param outfile Path for the combined tree file.
#' @return `outfile` (invisibly).
collect_candidate_trees <- function(window_results, outfile) {
  trees <- if (is.character(window_results)) window_results
           else vapply(window_results, function(r) r$tree, character(1))
  if (anyNA(trees))
    stop("Refusing to write NA candidate trees to ", outfile)
  writeLines(trees, outfile)
  invisible(outfile)
}

# ---------------------------------------------------------------------------
# Candidate-tree filtering
# ---------------------------------------------------------------------------

# Branch lengths at or below this are treated as absent when asking whether a
# tree carries any bipartition. IQ-TREE's minimum branch length is 1e-6, and a
# branch sitting on that floor carries no signal.
ZERO_BRANCH_TOL <- 1e-6

#' Drop candidate trees that MAST cannot use, in window order
#'
#' Applied between window tree estimation and the MAST fit. Drop reasons:
#' - `failed` — the window produced no tree (`ok = FALSE`).
#' - `unparsable` — a tree string `ape::read.tree()` cannot read.
#' - `taxon_set` — the tree is not on the alignment's full taxon set. IQ-TREE
#'   removes all-gap sequences, so a very gappy window yields a tree on a
#'   subset of taxa; MAST's `-te` requires every tree on the same taxa, which
#'   makes this a hard drop and the concrete form "too many gaps" takes.
#' - `uninformative` — no internal branch survives collapsing near-zero
#'   branches, i.e. a star tree carrying no bipartition.
#'
#' Duplicate topologies are deliberately **not** dropped. In MAST every tree
#' carries its own branch lengths, so a repeated topology is a real model
#' comparison (37 branch lengths + 1 weight on 20 unrooted taxa), not a free
#' parameter, and the IC rejects it unaided: across the 480-run validation
#' suite a duplicate reached the selected tree set in 4 runs, all at 3000
#' sites. They are counted instead (`n_distinct`), so a K_best leaning on a
#' redundant topology is visible rather than silent.
#'
#' @param window_results List of per-window results from
#'   `estimate_window_trees()`.
#' @param taxa Character vector of the alignment's taxon names.
#' @param tol Branch lengths at or below this are treated as zero.
#' @return List with `trees` (kept Newick strings, window order), `kept`
#'   (the kept window results), `reasons` (per window), `counts` (named
#'   integer vector of drop counts) and `n_distinct` (distinct topologies
#'   among the kept trees).
filter_candidate_trees <- function(window_results, taxa,
                                   tol = ZERO_BRANCH_TOL) {
  reasons <- vapply(window_results, function(r) {
    if (!isTRUE(r$ok) || is.na(r$tree)) return("failed")
    tr <- tryCatch(suppressWarnings(ape::read.tree(text = r$tree)),
                   error = function(e) NULL)
    if (is.null(tr) || is.null(tr$tip.label)) return("unparsable")
    if (!setequal(tr$tip.label, taxa)) return("taxon_set")
    if (is_uninformative_tree(tr, tol)) return("uninformative")
    "kept"
  }, character(1))

  keep   <- reasons == "kept"
  trees  <- vapply(window_results[keep], function(r) r$tree, character(1))
  levels <- c("failed", "unparsable", "taxon_set", "uninformative")

  list(
    trees      = trees,
    kept       = window_results[keep],
    reasons    = reasons,
    counts     = vapply(stats::setNames(levels, levels),
                        function(l) sum(reasons == l), integer(1)),
    n_distinct = count_distinct_topologies(trees, tol)
  )
}

#' Does a tree carry no bipartition at all?
#'
#' @param tr An `ape` phylo object.
#' @param tol Branch lengths at or below this are treated as zero.
#' @return `TRUE` for a star tree (after collapsing near-zero branches).
is_uninformative_tree <- function(tr, tol = ZERO_BRANCH_TOL) {
  tr <- tryCatch(ape::unroot(tr), error = function(e) tr)
  if (!is.null(tr$edge.length))
    tr <- tryCatch(ape::di2multi(tr, tol = tol), error = function(e) tr)
  isTRUE(tr$Nnode <= 1L)
}

#' Count distinct unrooted topologies in a set of Newick strings
#'
#' Compared after unrooting and collapsing near-zero branches, so two trees
#' differing only by a zero-length branch count once. Diagnostic only —
#' nothing is dropped on the strength of it.
#'
#' @param trees Character vector of Newick strings.
#' @param tol Branch lengths at or below this are treated as zero.
#' @return Number of distinct topologies, or `NA_integer_` if none could be
#'   compared.
count_distinct_topologies <- function(trees, tol = ZERO_BRANCH_TOL) {
  if (length(trees) == 0) return(0L)
  if (length(trees) == 1) return(1L)

  topos <- lapply(trees, function(t) tryCatch({
    tr <- ape::di2multi(ape::unroot(ape::read.tree(text = t)), tol = tol)
    tr$edge.length <- NULL
    tr
  }, error = function(e) NULL))
  topos <- topos[!vapply(topos, is.null, logical(1))]
  if (length(topos) == 0) return(NA_integer_)

  reps <- list()
  for (tr in topos) {
    dup <- FALSE
    for (rp in reps) {
      d <- tryCatch(as.numeric(ape::dist.topo(rp, tr, method = "PH85")),
                    error = function(e) NA_real_)
      if (length(d) == 1 && !is.na(d) && d == 0) { dup <- TRUE; break }
    }
    if (!dup) reps[[length(reps) + 1L]] <- tr
  }
  length(reps)
}

#' Describe drop counts for a message
#'
#' @param counts Named integer vector from `filter_candidate_trees()`.
#' @return A string like `"1 failed, 2 taxon_set"`, or `"none"`.
describe_drops <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) == 0) return("none")
  paste(paste(counts, names(counts)), collapse = ", ")
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
#' @param fixed_tree Tree handling for the K=1 single-tree fit: `"NJ"`,
#'   a file path, or `NULL` (heuristic search). Default `"NJ"`.
#' @param outdir Output directory.
#' @param label_prefix Prefix for per-K labels.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param threads Number of threads.
#' @param timeout Per-run timeout in seconds.
#' @param mast_max_fit Optional fit for `K = max(K_values)` from
#'   `mast_candidates()`, reused instead of refitting the same model.
#' @return Data frame with columns: K, lnL, df, AIC, AICc, BIC.
fit_mast_all_K <- function(alignment, K_values, base_model, rate_model,
                           tree_files, unlinked = FALSE,
                           fixed_tree = "NJ", outdir,
                           label_prefix = "",
                           iqtree_bin, threads, timeout,
                           mast_max_fit = NULL) {
  K_max <- max(K_values)
  results <- lapply(K_values, function(K) {
    label <- paste0(label_prefix, "K", K)

    if (!is.null(mast_max_fit) && K == K_max && K > 1) {
      # mast_candidates() already fitted all K_max trees to derive the
      # weights; the ranked tree set is the same set in a different order,
      # so the fit is identical.
      fit <- mast_max_fit
    } else if (K == 1) {
      # Standard single-tree fit
      model_str <- build_mast_model_str(base_model, rate_model, 1, unlinked)
      fit <- fit_model(
        alignment  = alignment,
        K          = 1,
        base_model = model_str,
        mix_type   = "+R",       # dummy — won't be appended for K = 1
        fixed_tree = fixed_tree,
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

# ---------------------------------------------------------------------------
# Candidate-tree derivation (run identically on empirical and replicate data)
# ---------------------------------------------------------------------------

#' Derive MAST candidate trees and tree sets from a single alignment
#'
#' Runs the full candidate-tree pipeline on whatever alignment it is given:
#' window split, per-window tree estimation, filtering of trees MAST cannot
#' use, MAST fit at the effective K_max to obtain tree weights, ranking, and
#' nested top-K tree sets.
#'
#' Bootstrap replicates must call this on their own alignment. Passing the
#' empirical tree sets to a replicate hands it the topologies it was simulated
#' from, which pins the IC minimum at K_best and forces power to 100%.
#'
#' @param alignment Path to the alignment.
#' @param K_values Integer vector of K values.
#' @param base_model Base substitution model string.
#' @param rate_model Within-class rate heterogeneity (e.g. `"+R4"`). When
#'   `NULL`, it is derived from the per-window model selections.
#' @param unlinked Logical; `TRUE` for `*T` (MIX syntax).
#' @param window_method Passed to `estimate_window_trees()`.
#' @param filter_trees Logical; when `TRUE` (default) candidate trees MAST
#'   cannot use are dropped by `filter_candidate_trees()` and `K_max` falls to
#'   the number of survivors. Windows are *not* re-split and trees are *not*
#'   re-inferred — the surviving trees are the ones already estimated. When
#'   `FALSE`, a failed window is an error, as it was before the filter existed.
#' @param outdir Output directory for this alignment's files.
#' @param iqtree_bin Path to IQ-TREE executable.
#' @param threads Number of threads.
#' @param timeout Per-run timeout in seconds.
#' @param seed Optional base seed for window tree searches.
#' @return List with: rate_model, all_trees (the kept trees), ranked,
#'   tree_files, mast_max (`NULL` when fewer than two trees survive),
#'   mast_dir, K_values_effective, K_max_effective, K_max_requested,
#'   drop_reasons, drop_counts, n_distinct.
mast_candidates <- function(alignment, K_values, base_model,
                            rate_model = NULL, unlinked = FALSE,
                            window_method = "NJ", filter_trees = TRUE,
                            outdir, iqtree_bin, threads, timeout,
                            seed = NULL) {
  K_max <- max(K_values)
  K_min <- min(K_values)

  windows <- split_alignment_windows(alignment, K_max, outdir)
  window_results <- estimate_window_trees(
    windows, outdir, iqtree_bin, threads, timeout,
    window_method = window_method, base_model = base_model, seed = seed
  )

  if (is.null(rate_model)) rate_model <- determine_rate_heterogeneity(window_results)

  if (filter_trees) {
    filt <- filter_candidate_trees(
      window_results, taxa = names(read_alignment(alignment))
    )
  } else {
    bad <- !vapply(window_results, function(r) isTRUE(r$ok), logical(1))
    if (any(bad))
      stop("Window tree estimation failed for window(s) ",
           paste(which(bad), collapse = ", "),
           ". Set filter_trees = TRUE to drop unusable windows and reduce ",
           "K_max to the survivors.")
    filt <- list(
      trees      = vapply(window_results, function(r) r$tree, character(1)),
      reasons    = rep("kept", length(window_results)),
      counts     = c(failed = 0L, unparsable = 0L, taxon_set = 0L,
                     uninformative = 0L),
      n_distinct = NA_integer_
    )
  }

  all_trees       <- filt$trees
  K_max_effective <- length(all_trees)

  if (K_max_effective < K_max)
    message("  Candidate trees: ", K_max_effective, " of ", K_max,
            " usable (dropped ", describe_drops(filt$counts),
            "); K_max reduced to ", max(K_max_effective, 1L))

  mast_dir <- file.path(outdir, "mast_fits")
  dir.create(mast_dir, showWarnings = FALSE, recursive = TRUE)

  # Fewer than two usable trees leaves no mixture to fit: the alignment gets
  # the single-tree K = 1 fit and nothing else.
  if (K_max_effective < 2) {
    if (K_max >= 2)
      warning("Only ", K_max_effective, " usable candidate tree(s) for ",
              alignment, "; fitting K = 1 only.")
    return(list(
      rate_model         = rate_model,
      all_trees          = all_trees,
      ranked             = integer(0),
      tree_files         = list(),
      mast_max           = NULL,
      mast_dir           = mast_dir,
      K_values_effective = 1L,
      K_max_effective    = K_max_effective,
      K_max_requested    = K_max,
      drop_reasons       = filt$reasons,
      drop_counts        = filt$counts,
      n_distinct         = filt$n_distinct
    ))
  }

  K_values_effective <- K_values[K_values <= K_max_effective]
  if (length(K_values_effective) == 0) K_values_effective <- K_max_effective

  tree_file <- collect_candidate_trees(
    all_trees, file.path(outdir, "candidate_trees.newick")
  )

  mast_max <- fit_mast_model(
    alignment  = alignment,
    tree_file  = tree_file,
    model_str  = build_mast_model_str(base_model, rate_model,
                                      K_max_effective, unlinked),
    outdir     = mast_dir,
    label      = paste0("mast_K", K_max_effective),
    iqtree_bin = iqtree_bin,
    threads    = threads,
    timeout    = timeout
  )

  # A NULL/short weight vector used to slip through as order(NULL) ==
  # integer(0), which writes literal NA lines into the tree files. Error here
  # instead of corrupting every downstream fit.
  if (length(mast_max$tree_weights) != K_max_effective ||
      anyNA(mast_max$tree_weights))
    stop("Could not parse ", K_max_effective, " MAST tree weights from ",
         mast_max$iqtree_file)

  ranked <- rank_trees_by_weight(mast_max$tree_weights)

  list(
    rate_model         = rate_model,
    all_trees          = all_trees,
    ranked             = ranked,
    tree_files         = build_mast_tree_files(all_trees, ranked,
                                               K_values_effective, outdir),
    mast_max           = mast_max,
    mast_dir           = mast_dir,
    K_values_effective = K_values_effective,
    K_max_effective    = K_max_effective,
    K_max_requested    = K_max,
    drop_reasons       = filt$reasons,
    drop_counts        = filt$counts,
    n_distinct         = filt$n_distinct
  )
}


#' Resolve the `+T` tree-search setting into its two component estimators
#'
#' The window trees and the K = 1 tree must be estimated the same way, or the
#' K = 1 versus K >= 2 comparison is biased by tree quality rather than by the
#' number of tree classes. This couples them so they cannot be set apart.
#'
#' @param tree_search Either `"NJ"` (BioNJ throughout) or `"fast"`
#'   (`--fast` heuristic search throughout).
#' @return List with `window_method` and `fixed_tree`.
resolve_tree_search <- function(tree_search = "NJ") {
  tree_search <- match.arg(tree_search, c("NJ", "fast"))
  list(
    window_method = tree_search,
    fixed_tree    = if (tree_search == "NJ") "NJ" else NULL
  )
}
