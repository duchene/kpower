# Low-level IQ-TREE invocation

#' Run IQ-TREE with given arguments
#'
#' Thin wrapper around `processx::run()`. Does not invoke a shell.
#'
#' @param iqtree_bin Path to the IQ-TREE executable.
#' @param args Character vector of arguments passed directly to IQ-TREE.
#' @param timeout Maximum run time in seconds (default 3600).
#' @return Invisibly, the `processx` result object (with $status, $stdout,
#'   $stderr).
run_iqtree <- function(iqtree_bin, args, timeout = 3600) {
  result <- processx::run(
    command         = iqtree_bin,
    args            = args,
    timeout         = timeout,
    error_on_status = FALSE
  )
  if (result$status != 0) {
    stop("IQ-TREE exited with status ", result$status, ":\n", result$stderr)
  }
  invisible(result)
}

#' Fit a mixture model with K categories to an alignment
#'
#' Calls IQ-TREE with `--fast` and returns parsed IC scores. The tree is not
#' fixed so that topology is co-estimated with the model (empirical fits only;
#' simulation refits should pass a fixed tree via `tree`).
#'
#' @param alignment Path to the input alignment (FASTA or PHYLIP).
#' @param K Integer number of mixture categories.
#' @param base_model Base substitution model string (e.g. `"GTR"`, `"LG"`).
#' @param mix_type Mixture type suffix: `"+R"` (FreeRate), `"+H"` (GHOST),
#'   etc. Default `"+R"`.
#' @param tree Optional path to a fixed tree file. If supplied, `--tree-fix`
#'   is added.
#' @param outdir Directory in which to write IQ-TREE output files.
#' @param label Short label used to build the `--prefix` (default derived from
#'   K and mix_type).
#' @param iqtree_bin Path to the IQ-TREE executable.
#' @param threads Number of threads (passed to `-T`; default `"AUTO"`).
#' @param timeout Seconds before the run is killed (default 3600).
#' @return Named list: K, model_string, lnL, df, AIC, AICc, BIC, prefix,
#'   treefile, logfile, iqtree_file.
fit_model <- function(alignment, K, base_model = "GTR", mix_type = "+R",
                      tree = NULL, outdir = tempdir(), label = NULL,
                      iqtree_bin = find_iqtree(), threads = "AUTO",
                      timeout = 3600) {
  model_str <- if (K == 1) base_model else paste0(base_model, mix_type, K)
  if (is.null(label)) label <- paste0("K", K)
  prefix <- make_prefix(outdir, label)

  args <- c(
    "-s", alignment,
    "-m", model_str,
    "--fast",
    "--prefix", prefix,
    "-T", as.character(threads),
    "--redo"
  )

  if (!is.null(tree)) {
    args <- c(args, "-t", tree, "--tree-fix")
  }

  run_iqtree(iqtree_bin, args, timeout = timeout)

  iqtree_file <- paste0(prefix, ".iqtree")
  report      <- parse_iqtree_report(iqtree_file)
  treefile    <- paste0(prefix, ".treefile")
  logfile     <- paste0(prefix, ".log")

  c(
    list(K = K, model_string = model_str, prefix = prefix,
         treefile = treefile, logfile = logfile, iqtree_file = iqtree_file),
    report
  )
}
