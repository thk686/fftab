#' Compute FFT and Return a Tidy Result with Normalized Coefficients
#'
#' Computes the Fast Fourier Transform (FFT) for various input types, including
#' vectors, time series (`ts`), and arrays. The FFT values are normalized by the
#' length of the input (or its product of dimensions for arrays). The result is
#' returned as a tidy tibble containing the Fourier frequencies and normalized
#' FFT values.
#'
#' @param x Input object for which to compute the FFT. This can be:
#'   - A numeric vector (default method).
#'   - A time series object (`ts`).
#'   - A multidimensional numeric array.
#'
#' @return A tibble containing:
#'   - **Fourier frequencies**: Represented by columns `dim_1`, `dim_2`, ..., depending on the input dimensions.
#'   - **Normalized FFT values**: Stored in the `fx` column as complex values, normalized by the input size.
#'
#' @details This is a generic function with methods for specific input types:
#'
#' - **Default Input (`tidy_fft.default`)**: Computes the FFT for a numeric vector and normalizes
#' the result by the vector length.
#' - **Time Series Input (`tidy_fft.ts`)**: Computes the FFT for a `ts` object, scaling the frequencies
#' by the time series frequency attribute and normalizing the FFT values.
#' - **Array Input (`tidy_fft.array`)**: Computes the FFT for a multidimensional array and normalizes
#' the FFT values by the product of the array dimensions.
#'
#' Each method returns a tidy tibble where the Fourier frequencies (`dim_1`,
#' `dim_2`, etc.) are paired with their corresponding normalized FFT values.
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
tidy_fft <- function(x) {
  UseMethod("tidy_fft")
}

#' @rdname tidy_fft
#' @export
tidy_fft.default <- function(x) {
  stopifnot(is.vector(x), is.numeric(x) || is.complex(x))
  fourier_frequencies(x) |>
    tibble::add_column(fx = stats::fft(x) / length(x)) |>
    .as_tidy_fft_obj(.is_complex = is.complex(x))
}

#' @rdname tidy_fft
#' @export
tidy_fft.ts <- function(x) {
  tidy_fft(as.vector(x)) |>
    dplyr::mutate(dim_1 = .fourier_frequencies(x) * frequency(x)) |>
    structure(.tsp = attr(x, "tsp"))
}

#' @rdname tidy_fft
#' @export
tidy_fft.array <- function(x) {
  fourier_frequencies(x) |>
    dplyr::mutate(fx = as.vector(stats::fft(x)) / prod(dim(x))) |>
    .as_tidy_fft_obj(.is_complex = is.complex(x), .dim = dim(x))
}

#' Perform Inverse FFT on a Tidy Result
#'
#' Computes the inverse Fast Fourier Transform (IFFT) to reconstruct the
#' original signal from a `tidy_fft` object or other FFT results. This is a
#' generic function with methods for default inputs and `tidy_fft` objects.
#'
#' @param x Input object containing FFT results. This can be:
#'   - A numeric vector of FFT coefficients.
#'   - A `tidy_fft` object with stored FFT results.
#'
#' @return A vector or time series object representing the reconstructed signal
#'   from the FFT results. If the original signal was real-valued, the IFFT
#'   returns the real part of the reconstruction.
#'
#' @details ### Methods:
#' - **Default Input (`tidy_ifft.default`)**: Performs the inverse FFT on a numeric vector of FFT coefficients (with a warning).
#' - **Tidy FFT Input (`tidy_ifft.tidy_fft`)**: Reconstructs the original signal from a `tidy_fft` object,
#' restoring attributes such as dimensions (`.dim`) and time series properties
#' (`.tsp`), if present.
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
  UseMethod("tidy_ifft")
}

#' @rdname tidy_ifft
#' @export
tidy_ifft.default <- function(x) {
  warning("Function tidy_ifft called on a non-tidy-fft object.")
  stats::fft(x, inverse = TRUE)
}

#' @rdname tidy_ifft
#' @export
tidy_ifft.tidy_fft <- function(x) {
  fx <- get_fx(x)
  res <- if (!is.null(attr(x, ".dim"))) {
    dim(fx) <- attr(x, ".dim")
    stats::fft(fx, inverse = TRUE)
  } else {
    stats::fft(fx, inverse = TRUE)
  }
  if (attr(x, ".is_complex") == FALSE) {
    res <- Re(res)
  }
  if (!is.null(attr(x, ".tsp"))) {
    attr(res, "tsp") <- attr(x, ".tsp")
    class(res) <- "ts"
  }
  res
}
