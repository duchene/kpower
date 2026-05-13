# IC profile figure

# Avoid R CMD check NOTEs for ggplot2 non-standard evaluation
utils::globalVariables(c("K", ".data", "lab", "replicate", "mix_type"))

#' Plot IC profiles from empirical data and parametric bootstrap simulations
#'
#' Produces a figure with K on the x-axis and the chosen information criterion
#' on the y-axis. Each simulated replicate is drawn as a thin, semi-transparent
#' line; the empirical IC profile is overlaid as a thicker line.
#'
#' Two calling modes:
#'
#' \itemize{
#'   \item **Single family** (back-compat): pass a single data frame to
#'     `empirical_ic` and `sim_ic`, plus scalar `K_best` and `power`.
#'   \item **Multi-family**: pass named lists (names = mixture-family labels,
#'     e.g. "+R", "+H", "+T") to `empirical_ic` and `sim_ic`, with matching
#'     named integer vector `K_best` and named numeric vector `power`. The
#'     plot then colours lines by family and reports per-family power in
#'     the subtitle.
#' }
#'
#' @param empirical_ic Data frame (single-family) or named list of data
#'   frames (multi-family) with columns K and the IC column.
#' @param sim_ic       Data frame or named list of data frames with columns
#'   replicate, K, and the IC column.
#' @param K_best       Integer (single) or named integer vector (multi).
#' @param power        Numeric in [0,1] (single) or named numeric vector
#'   (multi).
#' @param ic Character; which IC column to plot: `"AIC"`, `"AICc"`, or `"BIC"`
#'   (default `"BIC"`).
#' @return A ggplot2 object.
#' @export
plot_kpower <- function(empirical_ic, sim_ic, K_best, power, ic = "BIC") {

  # Detect multi-family mode: empirical_ic is a named list of data frames
  multi <- is.list(empirical_ic) && !is.data.frame(empirical_ic)

  if (!multi) {
    return(plot_kpower_single(empirical_ic, sim_ic, K_best, power, ic))
  }

  fams <- names(empirical_ic)
  if (is.null(fams) || any(!nzchar(fams)))
    stop("In multi-family mode, `empirical_ic` must be a named list.")

  # Stack everything into long tibbles with a mix_type column
  emp_long <- do.call(rbind, lapply(fams, function(mt) {
    df <- empirical_ic[[mt]]
    df$mix_type <- mt
    df[, c("mix_type", "K", ic)]
  }))
  sim_long <- do.call(rbind, lapply(fams, function(mt) {
    df <- sim_ic[[mt]]
    if (is.null(df)) return(NULL)
    df$mix_type <- mt
    df[, c("mix_type", "replicate", "K", ic)]
  }))

  emp_long$mix_type <- factor(emp_long$mix_type, levels = fams)
  if (!is.null(sim_long))
    sim_long$mix_type <- factor(sim_long$mix_type, levels = fams)

  # Per-family K_best marker rows
  best_rows <- do.call(rbind, lapply(fams, function(mt) {
    df <- empirical_ic[[mt]]
    row <- df[df$K == K_best[[mt]], , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    row$mix_type <- factor(mt, levels = fams)
    row[, c("mix_type", "K", ic)]
  }))

  subtitle <- paste0(
    "Power (", ic, "): ",
    paste(vapply(fams, function(mt)
      sprintf("%s %.0f%%", mt, 100 * power[[mt]]),
      character(1)),
      collapse = "  |  ")
  )

  family_palette <- c(
    "+R" = "#d6604d", "+H" = "#4393c3", "+T" = "#2ca02c",
    "*R" = "#a50026", "*H" = "#2166ac", "*T" = "#1a7c1a"
  )
  pal <- family_palette[fams]
  n_missing <- sum(is.na(pal))
  if (n_missing > 0)
    pal[is.na(pal)] <- scales::hue_pal()(n_missing)
  names(pal) <- fams

  p <- ggplot2::ggplot()

  if (!is.null(sim_long)) {
    p <- p + ggplot2::geom_line(
      data    = sim_long,
      mapping = ggplot2::aes(
        x = K, y = .data[[ic]],
        colour = mix_type,
        group  = interaction(mix_type, replicate)
      ),
      alpha = 0.10, linewidth = 0.3
    )
  }

  p <- p +
    ggplot2::geom_line(
      data    = emp_long,
      mapping = ggplot2::aes(x = K, y = .data[[ic]], colour = mix_type,
                             group = mix_type),
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data    = emp_long,
      mapping = ggplot2::aes(x = K, y = .data[[ic]], colour = mix_type),
      size = 2.5
    ) +
    ggplot2::geom_point(
      data    = best_rows,
      mapping = ggplot2::aes(x = K, y = .data[[ic]], colour = mix_type),
      shape = 21, fill = "white", size = 4, stroke = 1.5
    ) +
    ggplot2::scale_colour_manual(name = "Family", values = pal) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(emp_long$K))
    ) +
    ggplot2::labs(
      x = "Number of mixture categories (K)",
      y = ic,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}


# ---- single-family plot (original behaviour) ------------------------------
plot_kpower_single <- function(empirical_ic, sim_ic, K_best, power,
                               ic = "BIC") {
  power_label <- paste0(
    "Power = ", round(power * 100, 1), "% of simulations select K = ", K_best
  )

  legend_df <- data.frame(
    K   = empirical_ic$K[1],
    ic  = empirical_ic[[ic]][1],
    lab = power_label
  )

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data    = sim_ic,
      mapping = ggplot2::aes(x = K, y = .data[[ic]], group = replicate),
      colour  = "#4393c3", alpha = 0.15, linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data    = empirical_ic,
      mapping = ggplot2::aes(x = K, y = .data[[ic]]),
      colour  = "#d6604d", linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data    = empirical_ic,
      mapping = ggplot2::aes(x = K, y = .data[[ic]]),
      colour  = "#d6604d", size = 2.5
    ) +
    ggplot2::geom_point(
      data    = empirical_ic[empirical_ic$K == K_best, ],
      mapping = ggplot2::aes(x = K, y = .data[[ic]]),
      colour  = "#d6604d", shape = 21, fill = "white",
      size = 4, stroke = 1.5
    ) +
    ggplot2::geom_line(
      data    = legend_df,
      mapping = ggplot2::aes(x = K, y = ic, colour = lab),
      alpha   = 0
    ) +
    ggplot2::scale_colour_manual(name = NULL, values = c("black"),
                                 labels = power_label) +
    ggplot2::scale_x_continuous(breaks = empirical_ic$K) +
    ggplot2::labs(x = "Number of mixture categories (K)", y = ic) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
}
