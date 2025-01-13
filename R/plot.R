#' Plot the modulus of FFT results
#'
#' Plots the modulus of the FFT results against the frequencies.
#'
#' @param x A `fftab` object. This should contain Fourier-transformed data.
#' @param ... Additional arguments passed to `ggplot2::ggplot`.
#' @return A `ggplot` object representing the modulus of FFT results plotted against the frequencies.
#' The plot shows the modulus (`mod`) on the y-axis and frequency values on the x-axis.
#' @exportS3Method graphics::plot
plot.fftab <- function(x, ...) {
  to_polr(x) |>
    ggplot2::ggplot(...) +
    ggplot2::aes(x = .data$.dim_1, y = .data$mod) +
    ggplot2::geom_line() +
    ggplot2::ylab("modulus") +
    ggplot2::theme_classic()
}
