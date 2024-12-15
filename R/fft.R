#' Perform FFT and return a tidy result
#'
#' This is a generic function to compute the FFT and return the result in a tidy format.
#'
#' @param x Input object (vector, time series, or array).
#' @param repr The desired representation of the FFT result: `"polr"`, `"rect"`, or `"cplx"`.
#' @return A tibble containing Fourier frequencies and FFT values in the specified format.
#' @export
tidy_fft <- function(x, repr = c("polr", "rect", "cplx")) {
  UseMethod("tidy_fft")
}

#' Perform FFT for default inputs
#'
#' Computes the FFT for a numeric vector and returns a tidy result.
#'
#' @param x A numeric vector.
#' @param repr The desired representation of the FFT result.
#' @return A tibble containing Fourier frequencies and FFT values in the specified format.
#' @examples
#' tidy_fft(c(1, 0, -1, 0), repr = "rect")
#' @export
tidy_fft.default <- function(x, repr = c("polr", "rect", "cplx")) {
  repr <- match.arg(repr)
  if (!is.vector(x))
    stop('Input must be a vector.')
  fourier_frequencies(x) |>
    tibble::add_column(fx = stats::fft(x)) |>
    structure(
      class = c("tidy_fft", "tbl_df", "tbl", "data.frame"),
      is_complex = is.complex(x)
    ) |>
    change_repr(repr = match.arg(repr))
}

#' Perform FFT for time series inputs
#'
#' Computes the FFT for a time series object and returns a tidy result.
#'
#' @param x A time series (`ts`) object.
#' @param repr The desired representation of the FFT result.
#' @return A tibble containing Fourier frequencies and FFT values in the specified format.
#' @examples
#' ts_obj <- ts(sin(1:10), frequency = 12)
#' tidy_fft(ts_obj, repr = "polr")
#' @export
tidy_fft.ts <- function(x, repr = c("polr", "rect", "cplx")) {
  stride = attr(x, "tsp")[3]
  tidy_fft(as.vector(x), repr = match.arg(repr)) |>
    dplyr::mutate(dim_1 = .fourier_frequencies(x) * stride) |>
    structure(tsp_orig = attr(x, "tsp"))
}

#' Perform FFT for array inputs
#'
#' Computes the FFT for an array and returns a tidy result.
#'
#' @param x A numeric array.
#' @param repr The desired representation of the FFT result.
#' @return A tibble containing Fourier frequencies and FFT values in the specified format.
#' @examples
#' arr <- array(1:8, dim = c(2, 2, 2))
#' tidy_fft(arr, repr = "rect")
#' @export
tidy_fft.array <- function(x, repr = c("polr", "rect", "cplx")) {
  fourier_frequencies(x) |>
    dplyr::mutate(fx = as.vector(stats::fft(x))) |>
    structure(
      class = c("tidy_fft", "tbl_df", "tbl", "data.frame"),
      is_complex = is.complex(x),
      dim_orig = dim(x)
    ) |>
    change_repr(repr = match.arg(repr))
}

#' Perform inverse FFT on a tidy result
#'
#' This is a generic function for performing the inverse FFT on a `tidy_fft` object.
#'
#' @param x A `tidy_fft` object containing FFT results.
#' @return A vector or time series object reconstructed from the FFT results.
#' @export
tidy_ifft <- function(x) {
  UseMethod("tidy_ifft")
}

#' Perform inverse FFT for default inputs
#'
#' Performs the inverse FFT and reconstructs the original signal.
#'
#' @param x A `tidy_fft` object.
#' @return The reconstructed signal as a vector or time series object.
#' @examples
#' fft_res <- tidy_fft(c(1, 0, -1, 0), repr = "cplx")
#' tidy_ifft(fft_res)
#' @export
tidy_ifft.tidy_fft <- function(x) {
  fx <- change_repr(x, "cplx")$fx
  res <- if (!is.null(attr(x, "dim_orig"))) {
    dim(fx) <- attr(x, "dim_orig")
    stats::fft(fx, inverse = TRUE) / prod(dim(fx))
  } else {
    stats::fft(fx, inverse = TRUE) / length(fx)
  }
  if (attr(x, "is_complex") == FALSE)
    res <- Re(res)
  if (is.null(attr(x, "tsp_o$ig")) == FALSE) {
    attr(res, "tsp") <- attr(x, "tsp_orig")
    class(res) <- "ts"
  }
  res
}

