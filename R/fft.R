#' Compute FFT and Return a Tidy Result
#'
#' Computes the Fast Fourier Transform (FFT) for various input types, including
#' vectors, time series (`ts`), and arrays. The result is
#' returned as a tidy tibble containing the Fourier frequencies and
#' FFT values.
#'
#' @param x Input object for which to compute the FFT. This can be:
#'   - A numeric vector (default method).
#'   - A time series object (`ts`).
#'   - A multidimensional numeric array.
#'
#' @return A tibble containing:
#'   - **Fourier frequencies**: Represented by columns `dim_1`, `dim_2`, ..., depending on the input dimensions.
#'   - **FFT values**: Stored in the `fx` column as complex values.
#'
#' @details This is a generic function with methods for specific input types:
#'
#' - **Default Input (`tidy_fft.default`)**: Computes the FFT for a numeric vector.
#' - **Time Series Input (`tidy_fft.ts`)**: Computes the FFT for a `ts` object, scaling the frequencies
#' by the time series frequency attribute.
#' - **Array Input (`tidy_fft.array`)**: Computes the FFT for a multidimensional array.
#'
#' Each method returns a tidy tibble where the Fourier frequencies (`dim_1`,
#' `dim_2`, etc.) are paired with their corresponding FFT values.
#'
#' @examples
#' # FFT for a numeric vector
#' tidy_fft(c(1, 0, -1, 0))
#'
#' # FFT for a time series
#' ts_obj <- ts(sin(1:10), frequency = 12)
#' tidy_fft(ts_obj)
#'
#' # FFT for a multidimensional array
#' arr <- array(1:8, dim = c(2, 2, 2))
#' tidy_fft(arr)
#'
#' @seealso [fourier_frequencies()], [stats::fft()]
#' @export
tidy_fft <- function(x, norm = FALSE) {
  UseMethod("tidy_fft")
}

#' @rdname tidy_fft
#' @export
tidy_fft.default <- function(x, norm = FALSE) {
  stopifnot(is.vector(x), is.numeric(x) || is.complex(x))
  fourier_frequencies(x) |>
    tibble::add_column(fx = stats::fft(x)) |>
    .conditional_norm_fx(norm) |>
    .as_tidy_fft_obj(
      .is_normalized = norm,
      .is_complex = is.complex(x))
}

#' @rdname tidy_fft
#' @export
tidy_fft.ts <- function(x, norm = FALSE) {
  tidy_fft(as.vector(x), norm = norm) |>
    dplyr::mutate(dim_1 = dim_1 * frequency(x)) |>
    structure(.tsp = attr(x, "tsp"))
}

#' @rdname tidy_fft
#' @export
tidy_fft.array <- function(x, norm = FALSE) {
  fourier_frequencies(x) |>
    dplyr::mutate(fx = as.vector(stats::fft(x))) |>
    .conditional_norm_fx(norm) |>
    .as_tidy_fft_obj(.is_normalized = norm,
                     .is_complex = is.complex(x),
                     .dim = dim(x))
}

#' Perform Inverse FFT on a Tidy Result
#'
#' Computes the inverse Fast Fourier Transform (IFFT) to reconstruct the
#' original signal from a `tidy_fft` object.
#'
#' @param x A `tidy_fft` object with stored FFT results.
#'
#' @return A vector, array, or time series object representing the reconstructed
#'   signal from the FFT results. If the original signal was real-valued, the
#'   IFFT returns the real part of the reconstruction.
#'
#' @examples
#' # Example with FFT and inverse FFT
#' fft_result <- tidy_fft(c(1, 0, -1, 0))
#' original_signal <- tidy_ifft(fft_result)
#' print(original_signal)
#'
#' @seealso [stats::fft()], [tidy_fft()]
#' @export
tidy_ifft <- function(x) {
  stopifnot(inherits(x, "tidy_fft"))
  if (nrow(x) != .size(x))
    warning("Number of rows does not match the original object size.")
  fx <- get_fx(x)
  if (!.is_normalized(x))
    fx <- fx / length(fx)
  if (.is_array(x))
    dim(fx) <- .dim(x)
  res <- stats::fft(fx, inverse = TRUE)
  if (!.is_complex(x)) {
    res <- Re(res)
  }
  if (.is_ts(x)) {
    attr(res, "tsp") <- .tsp(x)
    class(res) <- "ts"
  }
  res
}
