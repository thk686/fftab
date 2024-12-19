#' Compute the Cross-Frequency Transform (Cross FFT)
#'
#' The `cross_fft` function computes the cross-frequency transform between two inputs
#' using the Fourier transform. Methods are provided for various input types, including
#' default, `ts`, `array`, and `tidy_fft` objects. The function allows an optional
#' parameter to control whether the conjugate of the second input is used.
#'
#' @param a The first input for the cross FFT. Type depends on the method.
#' @param b The second input for the cross FFT. Must match the dimensions or structure of `a`.
#' @param conj Logical; if `TRUE`, uses the complex conjugate of the Fourier transform of `b`.
#' Default is `TRUE`.
#'
#' @return An object representing the cross-frequency transform. The structure of the output
#' depends on the method used:
#' \itemize{
#'   \item `default`: A `tidy_fft` object.
#'   \item `ts`: A `tidy_fft` object with `tsp` attributes inherited from `a`.
#'   \item `array`: A `tidy_fft` object with `dim` attributes inherited from `a`.
#'   \item `tidy_fft`: A `tidy_fft` object combining the cross-frequency transform of `a` and `b`.
#' }
#'
#' @details
#' The `cross_fft` function is generic and relies on specific methods for different input types.
#' For real-valued time series (`ts`) or arrays, additional attributes are preserved where applicable.
#'
#' @examples
#' # Example with default method
#' a <- rnorm(8)
#' b <- rnorm(8)
#' result <- cross_fft(a, b)
#'
#' # Example with time series
#' ts_a <- ts(rnorm(8), frequency = 4)
#' ts_b <- ts(rnorm(8), frequency = 4)
#' result <- cross_fft(ts_a, ts_b)
#'
#' # Example with tidy_fft objects
#' tidy_a <- tidy_fft(a)
#' tidy_b <- tidy_fft(b)
#' result <- cross_fft(tidy_a, tidy_b)
#'
#' @seealso [tidy_fft()]
#'
#' @export
cross_fft <- function(a, b) {
  UseMethod("cross_fft")
}

#' @describeIn cross_fft Default method for computing cross FFT.
#' Converts inputs to `tidy_fft` objects before computation.
#' @export
cross_fft.default <- function(a, b, conj = TRUE) {
  cross_fft(tidy_fft(a), tidy_fft(b))
}

#' @describeIn cross_fft Method for time series (`ts`) objects.
#' Ensures the time series frequencies are consistent and preserves the `tsp` attribute.
#' @export
cross_fft.ts <- function(a, b, conj = TRUE) {
  stopifnot(frequency(a) == frequency(b))
  cross_fft(tidy_fft(a), tidy_fft(b)) |>
    structure(.tsp = attr(a, "tsp"))
}

#' @describeIn cross_fft Method for array inputs.
#' Ensures dimensions are consistent and preserves the `dim` attribute.
#' @export
cross_fft.array <- function(a, b, conj = TRUE) {
  stopifnot(dim(a) == dim(b))
  cross_fft(tidy_fft(a), tidy_fft(b)) |>
    structure(.dim = .dim(a))
}

#' @describeIn cross_fft Method for `tidy_fft` objects.
#' Performs the cross-frequency transform directly using the Fourier transforms of `a` and `b`.
#' @export
cross_fft.tidy_fft <- function(a, b, conj = TRUE) {
  stopifnot(nrow(a) == nrow(b))
  fx_b <- if (conj) Conj(get_fx(b)) else get_fx(b)
  fourier_frequencies(nrow(a)) |>
    tibble::add_column(fx = get_fx(a) * fx_b) |>
    .as_tidy_fft_obj(.is_complex = .is_complex(a) | .is_complex(b))
}

#'
#' #' Filter for Maximum Magnitude in a `tidy_fft` Object
#' #'
#' #' This function filters a `tidy_fft` object to find the frequency with the maximum
#' #' magnitude (`mod`), while respecting the dimensional and complex/real-valued structure.
#' #'
#' #' @param a A `tidy_fft` object containing FFT results.
#' #' @return A tibble containing a single row corresponding to the frequency with the maximum magnitude.
#' #' @details
#' #' The function first converts the input to the polar representation (`'polr'`) using
#' #' [change_repr()], then filters out rows where all frequency dimensions (`dim_*`) are zero
#' #' (for complex-valued input) or less than or equal to zero (for real-valued input).
#' #'
#' #' @examples
#' #' # Generate a tidy_fft object
#' #' a <- tidy_fft(c(1, 0, -1, 0))
#' #'
#' #' # Compute cross-spectrum (self-correlation for simplicity)
#' #' cross_result <- cross_spec(a, a)
#' #'
#' #' # Filter for maximum magnitude
#' #' max_mod_result <- filter_max_mod(cross_result)
#' #' print(max_mod_result)
#' #'
#' #' @seealso [cross_spec()] for computing the cross-spectrum.
#' #' @export
#' filter_max_mod <- function(a) {
#'   # Define the filter condition based on whether the input is complex
#'   filter_condition <- if (attr(a, "is_complex"))
#'     ~ . != 0
#'   else
#'     ~ . > 0
#'
#'   # Convert to polar representation and filter for maximum magnitude
#'   change_repr(a, "polr") |>
#'     dplyr::filter(if_all(starts_with("dim_"), filter_condition)) |>
#'     dplyr::filter(mod == max(mod))
#' }
#'
#' #' Compute Maximum Correlation and Phase Shift in Original Units
#' #'
#' #' This function computes the maximum correlation and corresponding phase shift
#' #' in the original time or space units from the cross-spectrum of two `tidy_fft` objects.
#' #'
#' #' @param a A `tidy_fft` object representing the first signal.
#' #' @param b A `tidy_fft` object representing the second signal.
#' #' @param freq_scale The scale of the original frequencies (e.g., sampling rate or spatial resolution).
#' #'        Defaults to 1.
#' #' @return A named numeric vector of length 2:
#' #' - `max_correlation`: The maximum normalized correlation.
#' #' - `shift_units`: The phase shift converted to the original time or space units.
#' #' @examples
#' #' # Example signals
#' #' x1 <- c(1, 0, -1, 0)
#' #' x2 <- c(0, 1, 0, -1)
#' #'
#' #' # Compute tidy_fft for the signals
#' #' fft_a <- tidy_fft(x1)
#' #' fft_b <- tidy_fft(x2)
#' #'
#' #' # Compute maximum correlation and shift in original units
#' #' result <- max_correlation_phase_units(fft_a, fft_b, freq_scale = 1)
#' #' print(result) # Returns max correlation and shift in original units
#' #'
#' #' @seealso [cross_spec()], [filter_max_mod()]
#' #' @export
#' max_correlation_phase_units <- function(a, b, freq_scale = 1) {
#'   # Ensure inputs are valid tidy_fft objects
#'   if (!inherits(a, "tidy_fft") || !inherits(b, "tidy_fft")) {
#'     stop("Inputs must be tidy_fft objects.")
#'   }
#'   if (!identical(nrow(a), nrow(b))) {
#'     stop("Inputs must have the same number of rows.")
#'   }
#'
#'   # Compute the cross-spectrum
#'   cross_spectrum <- cross_spec(a, b)
#'
#'   # Compute total variance (energy) for normalization
#'   var_a <- sum(Mod(get_fx(a)) ^ 2)
#'   var_b <- sum(Mod(get_fx(b)) ^ 2)
#'
#'   # Find the maximum magnitude and its corresponding phase
#'   dominant <- filter_max_mod(cross_spectrum)
#'
#'   # Compute normalized correlation
#'   max_correlation <- dominant$mod / sqrt(var_a * var_b)
#'
#'   # Extract the phase shift (radians)
#'   phase_shift <- dominant$arg
#'
#'   # Convert the phase shift to the original units
#'   frequency <- dominant$dim_1 * freq_scale  # Frequency in original units
#'   shift_units <- phase_shift / (2 * pi * frequency)
#'
#'   # Return results as a named vector
#'   c(max_correlation = max_correlation, shift_units = shift_units)
#' }
