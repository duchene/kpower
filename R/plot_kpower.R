# IC profile figure

#' Plot IC profiles from empirical data and parametric bootstrap simulations
#'
#' Produces a figure with K on the x-axis and the chosen information criterion
#' on the y-axis. Each simulated replicate is drawn as a thin, semi-transparent
#' line; the empirical IC profile is overlaid as a thicker line. The
#' proportion of simulations that recover K_best is shown in the legend.
#'
#' @param empirical_ic Data frame with columns K and the IC column (output of
#'   `fit_all_K()` on the empirical alignment).
#' @param sim_ic Long-format data frame with columns replicate, K, and the IC
#'   column (output of `assess_power()$sim_ic`).
#' @param K_best Integer; the K selected from empirical data.
#' @param power Numeric in [0, 1]; proportion of simulations recovering K_best.
#' @param ic Character; which IC column to plot: `"AIC"`, `"AICc"`, or `"BIC"`
#'   (default `"BIC"`).
#' @return A ggplot2 object.
#' @export
plot_kpower <- function(empirical_ic, sim_ic, K_best, power, ic = "BIC") {
  power_label <- paste0(
    "Power = ", round(power * 100, 1), "% of simulations select K = ", K_best
  )

  # Dummy data frame for the legend entry (invisible point off-plot range)
  legend_df <- data.frame(
    K   = empirical_ic$K[1],
    ic  = empirical_ic[[ic]][1],
    lab = power_label
  )

  p <- ggplot2::ggplot() +
    # Simulated IC profiles — thin, semi-transparent
    ggplot2::geom_line(
      data    = sim_ic,
      mapping = ggplot2::aes(
        x     = K,
        y     = .data[[ic]],
        group = replicate
      ),
      colour  = "#4393c3",
      alpha   = 0.05,
      linewidth = 0.3
    ) +
    # Empirical IC profile — thick, solid
    ggplot2::geom_line(
      data    = empirical_ic,
      mapping = ggplot2::aes(x = K, y = .data[[ic]]),
      colour  = "#d6604d",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data    = empirical_ic,
      mapping = ggplot2::aes(x = K, y = .data[[ic]]),
      colour  = "#d6604d",
      size    = 2.5
    ) +
    # Mark K_best on empirical line
    ggplot2::geom_point(
      data = empirical_ic[empirical_ic$K == K_best, ],
      mapping = ggplot2::aes(x = K, y = .data[[ic]]),
      colour = "#d6604d",
      shape  = 21,
      fill   = "white",
      size   = 4,
      stroke = 1.5
    ) +
    # Power annotation in legend via a transparent dummy layer
    ggplot2::geom_line(
      data    = legend_df,
      mapping = ggplot2::aes(x = K, y = ic, colour = lab),
      alpha   = 0
    ) +
    ggplot2::scale_colour_manual(
      name   = NULL,
      values = c("black"),
      labels = power_label
    ) +
    ggplot2::scale_x_continuous(breaks = empirical_ic$K) +
    ggplot2::labs(
      x = "Number of mixture categories (K)",
      y = ic
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}
