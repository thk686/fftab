# Declare .data as a global variable to avoid R CMD check NOTE
utils::globalVariables(c(
  ".data",
  "arg",
  ".dim_1",
  "fx",
  "im",
  "mod",
  "re",
  "frequency",
  ".correlation"
))

#' Compute the Fast Fourier Transform (FFT) of a Vector
#'
#' Computes the Fast Fourier Transform (FFT) of a numeric vector.
#' Optionally, normalizes the result by dividing it by the length of the input vector.
#'
#' @param x A numeric vector representing the input signal to transform.
#' @param norm A logical value indicating whether to normalize the FFT output
#'             by dividing it by the length of the input vector. Default is `FALSE`.
#'
#' @return A complex vector representing the FFT of the input signal.
#'
#' @details
#' This function wraps around the base R `stats::fft` function and provides an
#' option for normalization.
#'
#' @seealso [stats::fft()]
#' @keywords internal
.fft <- function(x, norm = FALSE) {
  if (norm) {
    stats::fft(x) / length(x)
  } else {
    stats::fft(x)
  }
}

#' Build a FFTAB Object
#'
#' Converts an object into a `fftab` object with additional metadata attributes.
#'
#' @param x The input object to convert, typically a tibble or data frame.
#' @param .is_normalized A logical value indicating if the input is normalized.
#' @param .is_complex A logical value indicating if the original data was complex.
#' @param ... Additional attributes to include in the structured object, such as
#'   metadata or specific attributes required for Fourier Transform analysis.
#'
#' @return The input object `x`, with the `fftab` class and any additional
#'   attributes provided in `...`.
#'
#' @keywords internal
.as_fftab_obj <- function(x, .is_normalized, .is_complex, ...) {
  structure(
    x,
    ...,
    .size = nrow(x),
    .is_angular = FALSE,
    .is_complex = .is_complex,
    .is_normalized = .is_normalized,
    class = c("fftab", class(x))
  )
}

#' Compute the number of samples in an input
#'
#' This helper function determines the number of samples in the input object.
#' For a vector, it returns its length. For a matrix or data frame, it returns
#' the number of rows.
#'
#' @param x An input object (scalar, vector, matrix, or data frame).
#' @return An integer representing the number of samples (rows) in the input
#'   object.
#' @keywords internal
.num_samples <- function(x) {
  input_len <- length(unlist(x))
  if (input_len == 0) {
    stop("Input must not be empty")
  }
  floor(ifelse(input_len == 1, x, NROW(x)))
}

#' Compute Fourier frequencies for default inputs
#'
#' Computes normalized Fourier frequencies for scalar or vector inputs, which
#' are evenly spaced between -0.5 and 0.5.
#'
#' @param x A scalar or vector representing the length of the sequence.
#' @return A numeric vector containing the normalized Fourier frequencies.
#' @keywords internal
.fourier_frequencies <- function(x) {
  n <- .num_samples(x)
  stopifnot(n > 0)
  if (n == 1) {
    return(0)
  }
  k <- 0:(n - 1)
  ifelse(k <= n / 2, k, k - n) / n
}

#' Select Columns Starting with `.dim_`
#'
#' Selects columns from a data frame whose names start with `.dim_`.
#'
#' @param x A data frame or tibble from which to select columns.
#'
#' @return A data frame containing only the columns that start with `.dim_`.
#'
#' @keywords internal
.get_dim_cols <- function(x) {
  dplyr::select(x, dplyr::starts_with(".dim_"))
}

#' Check if an Object is Normalized
#'
#' Checks whether an object has the `.is_normalized` attribute.
#'
#' @param x An R object to check for normalization.
#'
#' @return A logical value indicating whether the `.is_normalized` attribute is `TRUE`.
#'         Returns `NULL` if the attribute does not exist.
#'
#' @keywords internal
.is_normalized <- function(x) {
  attr(x, ".is_normalized")
}

#' Check if an Object is Complex
#'
#' Checks whether an object has the `.is_complex` attribute.
#'
#' @param x An R object to check for complexity.
#'
#' @return A logical value indicating whether the `.is_complex` attribute is `TRUE`.
#' @keywords internal
.is_complex <- function(x) {
  attr(x, ".is_complex")
}

#' Retrieve Object Dimensions
#'
#' Retrieves the `.dim` attribute from an object.
#'
#' @param x An R object to check for dimensions.
#'
#' @return The `.dim` attribute or `NULL` if not present.
#' @keywords internal
.dim <- function(x) {
  attr(x, ".dim")
}

#' Check if an Object is an Array
#'
#' Determines if an object has a `.dim` attribute, indicating it is an array.
#'
#' @param x An R object.
#'
#' @return `TRUE` if the object has a `.dim` attribute, `FALSE` otherwise.
#' @keywords internal
.is_array <- function(x) {
  !is.null(attr(x, ".dim"))
}

#' Retrieve Object Size
#'
#' Retrieves the `.size` attribute from an object.
#'
#' @param x An R object to check for size.
#'
#' @return The `.size` attribute or `NULL` if not present.
#' @keywords internal
.size <- function(x) {
  attr(x, ".size")
}

#' Retrieve Time Series Parameters
#'
#' Retrieves the `.tsp` attribute from an object.
#'
#' @param x An R object.
#'
#' @return The `.tsp` attribute or `NULL` if not present.
#' @keywords internal
.tsp <- function(x) {
  attr(x, ".tsp")
}

#' Check if an Object is a Time Series
#'
#' Determines if an object has a `.tsp` attribute.
#'
#' @param x An R object.
#'
#' @return `TRUE` if the object has a `.tsp` attribute, `FALSE` otherwise.
#' @keywords internal
.is_ts <- function(x) {
  !is.null(.tsp(x))
}

#' Retrieve Frequency
#'
#' Retrieves the frequency of a time series object.
#'
#' @param x A `fftab` object or time series.
#'
#' @return The frequency value or `1` if not a time series.
#' @keywords internal
.frequency <- function(x) {
  if (.is_ts(x)) {
    .tsp(x)[3]
  } else {
    1
  }
}

#' Compute Nyquist Frequency
#'
#' Computes the Nyquist frequency for an object.
#'
#' @param x An R object.
#'
#' @return The Nyquist frequency (`.frequency(x) / 2`).
#' @keywords internal
.nyquist <- function(x) {
  .frequency(x) / 2
}

#' Set Representation Format
#'
#' Converts an object into a specific representation format.
#'
#' @param x An R object.
#' @param repr The representation type (`"cplx"`, `"rect"`, `"polr"`).
#' @param .keep Controls preservation of data.
#'
#' @return The object in the specified representation format.
#' @keywords internal
.set_repr <- function(x, repr, .keep = "unused") {
  switch(repr,
    cplx = to_cplx(x, .keep = .keep),
    rect = to_rect(x, .keep = .keep),
    polr = to_polr(x, .keep = .keep),
    stop("Invalid representation.")
  )
}

#' Shift Phase
#'
#' Adjusts the phase of Fourier coefficients.
#'
#' @param x An R object.
#' @param shift A numeric value indicating the phase shift.
#' @keywords internal
.shift_phase <- function(x, shift) {
  if (.is_array(x)) {
    .NotYetImplemented()
  } else {
    if (.is_complex(x)) {
      to_polr(x) |> dplyr::mutate(arg = arg + shift)
    } else {
      to_polr(x) |>
        dplyr::mutate(
          arg = dplyr::case_when(
            .dim_1 == 0.0 ~ arg,
            .dim_1 == .nyquist(x) ~ arg,
            .dim_1 > 0.0 ~ arg + shift,
            .dim_1 < 0.0 ~ arg - shift,
            TRUE ~ arg
          )
        )
    }
  }
}

#' Compute Variance
#'
#' Computes the variance of Fourier coefficients.
#'
#' @param x An R object.
#' @return A numeric variance value.
#' @keywords internal
.variance <- function(x) {
  if (.is_normalized(x)) {
    to_polr(x) |>
      .remove_dc() |>
      dplyr::mutate(mod = mod^2) |>
      get_mod() |>
      sum()
  } else {
    to_polr(x) |>
      .remove_dc() |>
      dplyr::mutate(mod = mod^2) |>
      get_mod() |>
      sum() / .size(x)^2
  }
}

#' Compute Phase Difference Between Two Signals
#'
#' Computes the phase difference between two signals based on their cross-spectrum,
#' with symmetric and redundant frequency components removed.
#'
#' @param a A `fftab` object or signal representing the first input.
#' @param b A `fftab` object or signal representing the second input.
#'
#' @return A numeric value representing the **phase difference** (in radians) between the two signals.
#'
#' @details
#' This function computes the **cross-spectrum** of two signals, removes symmetric and redundant
#' frequency components, converts the result into polar representation, weights the phase angles
#' by their magnitudes, and calculates the weighted average phase difference.
#'
#' Removing symmetric components ensures accurate phase alignment, avoiding ambiguity caused
#' by redundant negative frequencies.
#'
#' @keywords internal
.phase_diff <- function(a, b) {
  cross_spec(b, a) |>
    to_polr() |>
    .remove_symmetric() |>
    dplyr::mutate(arg = mod * arg / sum(mod)) |>
    get_arg() |>
    sum()
}

#' Compute Normalized Correlation Between Two Signals
#'
#' Computes the normalized correlation between two signals based on their Fourier representations.
#'
#' @param a A `fftab` object or signal representing the first input.
#' @param b A `fftab` object or signal representing the second input.
#'
#' @return A numeric value representing the **normalized correlation** between the two signals.
#'
#' @details
#' This function computes the **cross-spectrum** of two signals, removes the DC component,
#' and calculates the real part of the cross-spectrum sum. The result is normalized using the
#' variances of both signals.
#'
#' Normalization ensures that the correlation value lies between -1 and 1.
#'
#' @keywords internal
.correlation <- function(a, b) {
  cross_spec(a, b, norm = TRUE) |>
    .remove_dc() |>
    get_re() |>
    sum() / sqrt(.variance(a) * .variance(b))
}

#' @keywords internal
.sort_dims <- function(x) {
  dplyr::arrange(x, dplyr::pick(dplyr::starts_with(".dim_")))
}

#' @keywords internal
.lt <- function(a, b) {
  i <- which.max(a != b)
  a[i] < b[i]
}

#' @keywords internal
.gt <- function(a, b) {
  i <- which.max(a != b)
  a[i] > b[i]
}

#' @keywords internal
.le <- function(a, b) {
  i <- which.max(a != b)
  a[i] <= b[i]
}

#' @keywords internal
.ge <- function(a, b) {
  i <- which.max(a != b)
  a[i] >= b[i]
}

#' Binary search for the DC row in lexicographically sorted data
#'
#' @param df A tibble or data.frame with `.dim_*` columns sorted lexicographically.
#' @return Row index of the DC component (0,0,...,0).
#' @keywords internal
.find_dc_row <- function(df) {
  left <- 1
  right <- nrow(df)
  dim_cols <- .get_dim_cols(df)
  zero_row <- rep(0, ncol(dim_cols))
  while (left <= right) {
    mid <- floor((left + right) / 2)
    row <- dim_cols[mid, ]
    if (all(row == zero_row)) {
      return(mid)
    } else if (.lt(row, zero_row)) {
      left <- mid + 1
    } else {
      right <- mid - 1
    }
  }
  stop("Input lacks a DC component.")
}


#' @title Remove DC Component and Symmetric Frequencies
#' @name fft_internal_filters
#' @description
#' Internal functions to manipulate and filter Fourier coefficients in `fftab` objects.
#'
#' @param x A `fftab` object containing Fourier coefficients and associated metadata.
#'
#' @return A `fftab` object with filtered coefficients.
#'
#' @details
#' - **`.remove_dc()`**: Filters out rows where all `.dim_*` columns have a value of `0`.
#' - **`.remove_symmetric()`**:
#'   - For real-valued signals, it filters out redundant, complex-conjugate frequencies.
#'   - For complex-valued signals, no filtering is applied as symmetry isn't relevant.
#' - **`.split_symmetric()`**: Splits the coefficients into symmetric and asymmetric parts.
#'
#' @seealso
#' - [dplyr::filter()]
#' - [fftab()]
#'
#' @keywords internal
NULL

#' @rdname fft_internal_filters
#' @keywords internal
.remove_dc <- function(x) {
  dplyr::filter(x, dplyr::if_any(dplyr::starts_with(".dim_"), ~ . != 0))
}

#' @rdname fft_internal_filters
#' @keywords internal
.get_dc <- function(x) {
  dplyr::filter(x, dplyr::if_all(dplyr::starts_with(".dim_"), ~ . == 0))
}

#' @rdname fft_internal_filters
#' @keywords internal
.remove_symmetric <- function(x) {
  if (.is_complex(x)) {
    return(x)
  }
  x <- .sort_dims(x)
  i <- .find_dc_row(x)
  dplyr::slice(x, i:nrow(x))
}

#' @rdname fft_internal_filters
#' @keywords internal
.split_symmetric <- function(x) {
  if (.is_complex(x)) {
    return(x)
  }
  x <- .sort_dims(x)
  i <- .find_dc_row(x)
  list(
    symmetric = dplyr::slice(x, 1:(i - 1)),
    asymmetric = dplyr::slice(x, i:nrow(x))
  )
}

#' @title Check Whether an Object is Angular
#'
#' @description
#' An internal helper function that checks whether an object `x` is marked as
#' using angular frequency. It looks for the `.is_angular` attribute.
#'
#' @details
#' This function is not exported and is only used internally by `to_angf()`
#' and `to_cycf()`.
#'
#' @keywords internal
.is_angular <- function(x) {
  attr(x, ".is_angular")
}

