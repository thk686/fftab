#' Perform FFT and return a tidy result
#'
#' This is a generic function to compute the FFT and return the result in a tidy format.
#'
#' @param x Input object (vector, time series, or array).
#' @return A tibble containing Fourier frequencies and FFT values in the specified format.
#' @export
tidy_fft <- function(x) {
  UseMethod("tidy_fft")
}

#' Perform FFT for default inputs
#'
#' Computes the FFT for a numeric vector and returns a tidy result.
#'
#' @param x A numeric vector.
#' @return A tibble containing Fourier frequencies and FFT values in the specified format.
#' @examples
#' tidy_fft(c(1, 0, -1, 0))
#' @export
tidy_fft.default <- function(x) {
  if (!is.vector(x) || !is.numeric(x))
    stop('Input must be a numeric vector.')
  fourier_frequencies(x) |>
    tibble::add_column(fx = stats::fft(x)) |>
    .as_tidy_fft_obj(is_complex = is.complex(x))
}

#' Perform FFT for time series inputs
#'
#' Computes the FFT for a time series object and returns a tidy result.
#'
#' @param x A time series (`ts`) object.
#' @return A tibble containing Fourier frequencies and FFT values in the specified format.
#' @examples
#' ts_obj <- ts(sin(1:10), frequency = 12)
#' tidy_fft(ts_obj)
#' @export
tidy_fft.ts <- function(x) {
  strides = attr(x, "tsp")[3]
  tidy_fft(as.vector(x)) |>
    dplyr::mutate(dim_1 = .fourier_frequencies(x) * strides) |>
    structure(.tsp = attr(x, "tsp"))
}

#' Perform FFT for array inputs
#'
#' Computes the FFT for an array and returns a tidy result.
#'
#' @param x A numeric array.
#' @return A tibble containing Fourier frequencies and FFT values in the specified format.
#' @examples
#' arr <- array(1:8, dim = c(2, 2, 2))
#' tidy_fft(arr)
#' @export
tidy_fft.array <- function(x) {
  fourier_frequencies(x) |>
    dplyr::mutate(fx = as.vector(stats::fft(x))) |>
    .as_tidy_fft_obj(is_complex = is.complex(x), .dim = dim(x))
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

#' @export
tidy_ifft.default <- function(x) {
  stats::fft(x, inverse = TRUE)
}

#' Perform inverse FFT for default inputs
#'
#' Performs the inverse FFT and reconstructs the original signal.
#'
#' @param x A `tidy_fft` object.
#' @return The reconstructed signal as a vector or time series object.
#' @examples
#' fft_res <- tidy_fft(c(1, 0, -1, 0))
#' tidy_ifft(fft_res)
#' @export
tidy_ifft.tidy_fft <- function(x) {
  fx <- change_repr(x, "cplx")$fx
  res <- if (!is.null(attr(x, ".dim"))) {
    dim(fx) <- attr(x, ".dim")
    stats::fft(fx, inverse = TRUE) / prod(dim(fx))
  } else {
    stats::fft(fx, inverse = TRUE) / length(fx)
  }
  if (attr(x, "is_complex") == FALSE)
    res <- Re(res)
  if (is.null(attr(x, ".tsp")) == FALSE) {
    attr(res, "tsp") <- attr(x, ".tsp")
    class(res) <- "ts"
  }
  res
}
