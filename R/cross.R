#' Compute the Cross-Spectrum
#'
#' This function computes the cross-spectrum between two `tidy_fft` objects by
#' multiplying the Fourier coefficients (`fx`) of the first object (`a`) with
#' the conjugate of the coefficients from the second object (`b`).
#'
#' @param a A `tidy_fft` object.
#' @param b A `tidy_fft` object with the same number of rows as `a`.
#' @return A `tidy_fft` object with updated `fx` values representing the cross-spectrum.
#' @details
#' The cross-spectrum is calculated as:
#' \deqn{S_{ab}(k) = X_a(k) \cdot \overline{X_b(k)}}
#' where \eqn{X_a(k)} and \eqn{X_b(k)} are the Fourier coefficients of `a` and `b`,
#' respectively, and \eqn{\overline{X_b(k)}} is the complex conjugate of \eqn{X_b(k)}.
#'
#' @examples
#' # Generate example tidy_fft objects
#' a <- tidy_fft(c(1, 0, -1, 0))
#' b <- tidy_fft(c(0, 1, 0, -1))
#'
#' # Compute the cross-spectrum
#' cross_result <- cross_spec(a, b)
#' print(cross_result)
#'
#' @seealso [filter_max_mod()] to filter the resulting cross-spectrum for the maximum magnitude.
#' @export
cross_spec <- function(a, b) {
    stopifnot(inherits(a, "tidy_fft"), inherits(b, "tidy_fft"), identical(nrow(a),
        nrow(b)))
    dplyr::mutate(a, fx = get_fx(a) * Conj(get_fx(b)))
}

#' Filter for Maximum Magnitude in a `tidy_fft` Object
#'
#' This function filters a `tidy_fft` object to find the frequency with the maximum
#' magnitude (`mod`), while respecting the dimensional and complex/real-valued structure.
#'
#' @param a A `tidy_fft` object containing FFT results.
#' @return A tibble containing a single row corresponding to the frequency with the maximum magnitude.
#' @details
#' The function first converts the input to the polar representation (`'polr'`) using
#' [change_repr()], then filters out rows where all frequency dimensions (`dim_*`) are zero
#' (for complex-valued input) or less than or equal to zero (for real-valued input).
#'
#' @examples
#' # Generate a tidy_fft object
#' a <- tidy_fft(c(1, 0, -1, 0))
#'
#' # Compute cross-spectrum (self-correlation for simplicity)
#' cross_result <- cross_spec(a, a)
#'
#' # Filter for maximum magnitude
#' max_mod_result <- filter_max_mod(cross_result)
#' print(max_mod_result)
#'
#' @seealso [cross_spec()] for computing the cross-spectrum.
#' @export
filter_max_mod <- function(a) {
    # Define the filter condition based on whether the input is complex
    filter_condition <- if (attr(a, "is_complex"))
        ~. != 0 else ~. > 0

    # Convert to polar representation and filter for maximum magnitude
    change_repr(a, "polr") |>
        dplyr::filter(if_all(starts_with("dim_"), filter_condition)) |>
        dplyr::filter(mod == max(mod))
}

#' Compute Maximum Correlation and Phase Shift in Original Units
#'
#' This function computes the maximum correlation and corresponding phase shift
#' in the original time or space units from the cross-spectrum of two `tidy_fft` objects.
#'
#' @param a A `tidy_fft` object representing the first signal.
#' @param b A `tidy_fft` object representing the second signal.
#' @param freq_scale The scale of the original frequencies (e.g., sampling rate or spatial resolution).
#'        Defaults to 1.
#' @return A named numeric vector of length 2:
#' - `max_correlation`: The maximum normalized correlation.
#' - `shift_units`: The phase shift converted to the original time or space units.
#' @examples
#' # Example signals
#' x1 <- c(1, 0, -1, 0)
#' x2 <- c(0, 1, 0, -1)
#'
#' # Compute tidy_fft for the signals
#' fft_a <- tidy_fft(x1)
#' fft_b <- tidy_fft(x2)
#'
#' # Compute maximum correlation and shift in original units
#' result <- max_correlation_phase_units(fft_a, fft_b, freq_scale = 1)
#' print(result) # Returns max correlation and shift in original units
#'
#' @seealso [cross_spec()], [filter_max_mod()]
#' @export
max_correlation_phase_units <- function(a, b, freq_scale = 1) {
    # Ensure inputs are valid tidy_fft objects
    if (!inherits(a, "tidy_fft") || !inherits(b, "tidy_fft")) {
        stop("Inputs must be tidy_fft objects.")
    }
    if (!identical(nrow(a), nrow(b))) {
        stop("Inputs must have the same number of rows.")
    }

    # Compute the cross-spectrum
    cross_spectrum <- cross_spec(a, b)

    # Compute total variance (energy) for normalization
    var_a <- sum(Mod(get_fx(a))^2)
    var_b <- sum(Mod(get_fx(b))^2)

    # Find the maximum magnitude and its corresponding phase
    dominant <- filter_max_mod(cross_spectrum)

    # Compute normalized correlation
    max_correlation <- dominant$mod/sqrt(var_a * var_b)

    # Extract the phase shift (radians)
    phase_shift <- dominant$arg

    # Convert the phase shift to the original units
    frequency <- dominant$dim_1 * freq_scale  # Frequency in original units
    shift_units <- phase_shift/(2 * pi * frequency)

    # Return results as a named vector
    c(max_correlation = max_correlation, shift_units = shift_units)
}
