#' Compute the Cross-Spectrum (Cross FFT)
#'
#' The `cross_spec` function computes the cross-spectrum between two inputs
#' using the Fourier transform. It supports multiple input types including
#' numeric vectors, time series (`ts`), arrays, and `fftab` objects.
#' The function provides options for normalization and controlling whether the
#' conjugate of the second input is used.
#'
#' @param a The first input for the cross FFT. Supported types include numeric
#' vectors, `ts` objects, arrays, and `fftab` objects.
#' @param b The second input for the cross FFT. Must match the dimensions or
#' structure of `a`.
#' @param norm Logical; if `TRUE`, normalizes the Fourier transforms before
#' computation. Default is `FALSE`.
#' @param conj Logical; if `TRUE`, uses the complex conjugate of the Fourier
#' transform of `b`. Default is `TRUE`.
#'
#' @return An object representing the cross-spectrum:
#' \itemize{
#'   \item For `default` and `fftab` methods: A `fftab` object.
#'   \item For `ts` objects: A `fftab` object with `.tsp` attributes inherited
#'         from `a`.
#'   \item For arrays: A `fftab` object with `.dim` attributes inherited
#'         from `a`.
#' }
#'
#' @examples
#' cross_spec(rnorm(8), rnorm(8), norm = TRUE)
#'
#' cross_spec(
#'   ts(rnorm(8), frequency = 4),
#'   ts(rnorm(8), frequency = 4)
#' )
#'
#' @seealso [fftab()]
#'
#' @export
cross_spec <- function(a, b, norm = FALSE, conj = TRUE) {
  UseMethod("cross_spec")
}

#' @describeIn cross_spec Default method for computing cross FFT.
#' Converts inputs to `fftab` objects before computation.
#' @export
cross_spec.default <- function(a, b, norm = FALSE, conj = TRUE) {
  cross_spec(fftab(a), fftab(b), norm = norm, conj = conj)
}

#' @describeIn cross_spec Method for time series (`ts`) objects.
#' Ensures the time series frequencies are consistent and preserves the `tsp` attribute.
#' @export
cross_spec.ts <- function(a, b, norm = FALSE, conj = TRUE) {
  stopifnot(frequency(a) == frequency(b))
  cross_spec(fftab(a), fftab(b), norm = norm, conj = conj) |>
    structure(.tsp = attr(a, "tsp"))
}

#' @describeIn cross_spec Method for array inputs.
#' Ensures dimensions are consistent and preserves the `dim` attribute.
#' @export
cross_spec.array <- function(a, b, norm = FALSE, conj = TRUE) {
  stopifnot(dim(a) == dim(b))
  cross_spec(fftab(a), fftab(b), norm = norm, conj = conj) |>
    structure(.dim = dim(a))
}

#' @describeIn cross_spec Method for `fftab` objects.
#' Performs the cross-frequency transform directly using the Fourier transforms of `a` and `b`.
#' @export
cross_spec.fftab <- function(a, b, norm = FALSE, conj = TRUE) {
  stopifnot(nrow(a) == nrow(b))
  fx_a <- get_fx_norm(a, norm)
  fx_b <- get_fx_norm(b, norm)
  if (conj) {
    fx_b <- Conj(fx_b)
  }
  .get_dim_cols(a) |>
    tibble::add_column(fx = fx_a * fx_b) |>
    structure(
      .is_normalized = norm,
      .is_complex = .is_complex(a) | .is_complex(b)
    )
}

#' Compute Phase Difference and Maximum Correlation Between Two Signals
#'
#' Computes the phase difference and maximum normalized correlation between two input signals
#' after phase-aligning the second signal (`b`) to the first signal (`a`).
#'
#' `r lifecycle::badge('experimental')`
#'
#' @param a A numeric vector or time series representing the first signal.
#' @param b A numeric vector or time series representing the second signal.
#'
#' @return A numeric vector of length two:
#'   - The first element represents the **phase difference** (in radians) required to maximize alignment between the two signals.
#'   - The second element represents the **maximum normalized correlation** achieved after phase alignment.
#'
#' @details
#' This function performs the following steps:
#' 1. Computes the Fourier Transform of both input signals using `fftab`.
#' 2. Calculates the **cross-spectrum** of the signals.
#' 3. Converts the cross-spectrum to polar form and computes the weighted average phase difference.
#' 4. Adjusts the phase of the second signal (`b`) using `.shift_phase` to maximize alignment with the first signal (`a`).
#' 5. Computes the **normalized correlation** between the phase-aligned signals.
#'
#' The correlation is normalized using the variances of both signals and will generally be **higher** than the correlation
#' between the original signals due to the optimal phase alignment.
#'
#' @seealso
#' - [fftab()]
#' - [cross_spec()]
#'
#' @examples
#' phase_diff(
#'   sin(seq(0, 2 * pi, length.out = 128)),
#'   cos(seq(0, 2 * pi, length.out = 128))
#' )
#'
#' @importFrom lifecycle badge
#' @export
phase_diff <- function(a, b) {
  if (!inherits(a, "fftab")) {
    a <- fftab(a)
  }
  if (!inherits(b, "fftab")) {
    b <- fftab(b)
  }
  phase_diff <- .phase_diff(a, b)
  b <- .shift_phase(b, -phase_diff)
  cor <- .correlation(a, b)
  c(phase_diff = phase_diff, correlation = cor)
}
