#' Compute the Cross-Spectrum (Cross FFT)
#'
#' The `cross_spec` function computes the cross-spectrum between two inputs
#' using the Fourier transform. It supports multiple input types including
#' numeric vectors, time series (`ts`), arrays, and `tidy_fft` objects.
#' The function provides options for normalization and controlling whether the
#' conjugate of the second input is used.
#'
#' @param a The first input for the cross FFT. Supported types include numeric
#' vectors, `ts` objects, arrays, and `tidy_fft` objects.
#' @param b The second input for the cross FFT. Must match the dimensions or
#' structure of `a`.
#' @param norm Logical; if `TRUE`, normalizes the Fourier transforms before
#' computation. Default is `FALSE`.
#' @param conj Logical; if `TRUE`, uses the complex conjugate of the Fourier
#' transform of `b`. Default is `TRUE`.
#'
#' @return An object representing the cross-spectrum:
#' \itemize{
#'   \item For `default` and `tidy_fft` methods: A `tidy_fft` object.
#'   \item For `ts` objects: A `tidy_fft` object with `tsp` attributes inherited
#'         from `a`.
#'   \item For arrays: A `tidy_fft` object with `dim` attributes inherited
#'         from `a`.
#' }
#'
#' @details
#' The `cross_spec` function is generic and has specialized methods for each
#' input type. For `ts` and array inputs, additional attributes such as
#' `tsp` or `dim` are preserved. The `%cfft%` operator provides a shorthand
#' for calling this function.
#'
#' @examples
#' # Example with numeric vectors
#' a <- rnorm(8)
#' b <- rnorm(8)
#' result <- cross_spec(a, b)
#'
#' # Example with the infix operator
#' result <- a %cfft% b
#'
#' # Example with time series
#' ts_a <- ts(rnorm(8), frequency = 4)
#' ts_b <- ts(rnorm(8), frequency = 4)
#' result <- cross_spec(ts_a, ts_b)
#'
#' # Example with `tidy_fft` objects
#' tidy_a <- tidy_fft(a)
#' tidy_b <- tidy_fft(b)
#' result <- cross_spec(tidy_a, tidy_b)
#'
#' @seealso [tidy_fft()]
#'
#' @aliases %cfft%
#' @export
cross_spec <- function(a, b, norm = FALSE, conj = TRUE) {
  UseMethod("cross_spec")
}

#' @describeIn cross_spec Default method for computing cross FFT.
#' Converts inputs to `tidy_fft` objects before computation.
#' @export
cross_spec.default <- function(a, b, norm = FALSE, conj = TRUE) {
  cross_spec(tidy_fft(a), tidy_fft(b), norm = norm, conj = conj)
}

#' @describeIn cross_spec Method for time series (`ts`) objects.
#' Ensures the time series frequencies are consistent and preserves the `tsp` attribute.
#' @export
cross_spec.ts <- function(a, b, norm = FALSE, conj = TRUE) {
  stopifnot(frequency(a) == frequency(b))
  cross_spec(tidy_fft(a), tidy_fft(b), norm = norm, conj = conj) |>
    structure(.tsp = attr(a, "tsp"))
}

#' @describeIn cross_spec Method for array inputs.
#' Ensures dimensions are consistent and preserves the `dim` attribute.
#' @export
cross_spec.array <- function(a, b, norm = FALSE, conj = TRUE) {
  stopifnot(dim(a) == dim(b))
  cross_spec(tidy_fft(a), tidy_fft(b), norm = norm, conj = conj) |>
    structure(.dim = dim(a))
}

#' @describeIn cross_spec Method for `tidy_fft` objects.
#' Performs the cross-frequency transform directly using the Fourier transforms of `a` and `b`.
#' @export
cross_spec.tidy_fft <- function(a, b, norm = FALSE, conj = TRUE) {
  stopifnot(nrow(a) == nrow(b))
  fx_a <- get_fx_norm(a, norm)
  fx_b <- get_fx_norm(b, norm)
  if (conj) {
    fx_b <- Conj(fx_b)
  }
  .get_dim_cols(a) |>
    tibble::add_column(fx = fx_a * fx_b) |>
    structure(.is_normalized = norm,
              .is_complex = .is_complex(a) | .is_complex(b))
}

phase_diff <- function(a, b) {
  fa <- tidy_fft(a)
  fb <- tidy_fft(b)
  phase_diff <- cross_spec(fb, fa) |>
    to_polr() |>
    dplyr::mutate(arg = mod * arg / sum(mod)) |>
    get_arg() |>
    sum()
  fb <- shift_phase(fb, phase_diff)
  cor <- cross_spec(fa, fb, norm = TRUE) |>
    get_re() |>
    sum()
  cor <- cor / sqrt(.variance(fa) * .variance(fb))
  c(phase_diff, cor)
}
