#' Plot the modulus of FFT results
#'
#' Plots the modulus of the FFT results against the frequencies.
#'
#' @param x A `fftab` object.
#' @param ... passed to ggplot.
#' @exportS3Method graphics::plot
plot.fftab <- function(x, ...) {
  to_polr(x) |>
    ggplot2::ggplot(...) +
    ggplot2::aes(x = .data$.dim_1, y = .data$mod) +
    ggplot2::geom_line() +
    ggplot2::ylab("modulus") +
    ggplot2::theme_classic()
}
