#' Compute the number of samples in an input
#'
#' This helper function determines the number of samples in the input object.
#'
#' @param x An input object (scalar, vector, or matrix).
#' @return The number of samples (rows) in the input object.
#' @keywords internal
.num_samples <- function(x) {
  input_len <- length(unlist(x))
  if (is.null(x) || input_len == 0)
    stop("Input must not be NULL or empty")
  floor(ifelse(input_len == 1, x, NROW(x)))
}

#' Compute Fourier frequencies
#'
#' This is a generic function for computing Fourier frequencies.
#'
#' @param x The input (scalar, vector, or time series object).
#' @return A vector of Fourier frequencies, normalized or scaled based on the input type.
#' @export
fourier_frequencies <- function(x) {
  UseMethod("fourier_frequencies")
}

#' Compute Fourier frequencies for default inputs
#'
#' Computes normalized Fourier frequencies for scalar or vector inputs.
#'
#' @param x A scalar or vector representing the length of the sequence.
#' @return A numeric vector of normalized Fourier frequencies.
#' @examples
#' fourier_frequencies(8)
#' @export
fourier_frequencies.default <- function(x) {
  n <- .num_samples(x)
  if (n < 2)
    stop("Minimum length is 2")
  k <- 0:(n - 1)
  return(ifelse(k <= n / 2, k, k - n) / n)
}

#' Compute Fourier frequencies for time series objects
#'
#' Computes Fourier frequencies scaled by the frequency attribute of a `ts` object.
#'
#' @param x A time series (`ts`) object.
#' @return A numeric vector of scaled Fourier frequencies.
#' @examples
#' ts_obj <- ts(rnorm(12), frequency = 12)
#' fourier_frequencies(ts_obj)
#' @export
fourier_frequencies.ts <- function(x) {
  n <- length(x)
  freq <- attr(x, "tsp")[3]
  return(fourier_frequencies(n) * freq)
}

#' Perform FFT and return a tidy result
#'
#' This is a generic function to compute the FFT and return the result in a tidy format.
#'
#' @param x Input object (vector or `ts` object).
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
    stop('Only vector input is supported at this time.')
  res <- tibble::tibble(frequency = fourier_frequencies(x))
  fx <- stats::fft(x, inverse = FALSE)
  if (repr == "cplx")
    res$fx <- fx
  if (repr == "rect") {
    res$re <- Re(fx)
    res$im <- Im(fx)
  }
  if (repr == "polr") {
    res$mod <- Mod(fx)
    res$arg <- Arg(fx)
  }
  return(structure(
    res,
    class = c("tidy_fft", class(res)),
    is_complex = is.complex(x),
    repr = repr
  ))
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
  res <- tidy_fft(as.vector(x), repr)
  res$frequency <- fourier_frequencies(x)
  attr(res, "tsp_orig") <- attr(x, "tsp")
  return(res)
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
  x <- change_repr(x, "cplx")
  res <- stats::fft(x$fx, inverse = TRUE) / nrow(x)
  if (attr(x, "is_complex") == FALSE)
    res <- Re(res)
  if (!is.null(attr(x, "tsp_orig"))) {
    attr(res, "tsp") <- attr(x, "tsp_orig")
    class(res) <- "ts"
  }
  return(res)
}

#' Retrieve the current representation of a tidy_fft object
#'
#' @param x A `tidy_fft` object.
#' @return The current representation (`"polr"`, `"rect"`, or `"cplx"`).
#' @export
repr <- function(x) {
  return(attr(x, "repr"))
}

# Internal function to convert between FFT representations
# This function is used by `change_repr`.
.convert_fft <- function(x, from, to) {
  with(x, {
    if (from == "cplx" && to == "rect")
      return(tibble::tibble(re = Re(fx), im = Im(fx)))
    if (from == "cplx" && to == "polr")
      return(tibble::tibble(mod = Mod(fx), arg = Arg(fx)))
    if (from == "rect" && to == "cplx")
      return(tibble::tibble(fx = complex(
        real = re, imaginary = im
      )))
    if (from == "rect" && to == "polr")
      return(tibble::tibble(mod = sqrt(re ^ 2 + im ^ 2), arg = atan2(im, re)))
    if (from == "polr" && to == "cplx")
      return(tibble::tibble(fx = complex(
        modulus = mod, argument = arg
      )))
    if (from == "polr" && to == "rect")
      return(tibble::tibble(re = mod * cos(arg), im = mod * sin(arg)))
    stop("Unsupported conversion.")
  })
}

#' Change the representation of FFT results
#'
#' Converts FFT results between different representations.
#'
#' @param x A `tidy_fft` object.
#' @param repr The target representation (`"polr"`, `"rect"`, or `"cplx"`).
#' @return The modified object with the new representation.
#' @examples
#' fft_res <- tidy_fft(c(1, 0, -1, 0), repr = "cplx")
#' change_repr(fft_res, "rect")
#' @export
change_repr <- function(x, repr = c("polr", "rect", "cplx")) {
  if (!inherits(x, "tidy_fft"))
    stop("Input must be of class 'tidy_fft'")
  repr <- match.arg(repr)
  from <- attr(x, "repr")
  if (from == repr)
    return(x)
  res <- .convert_fft(x, from, repr)
  res$frequency <- x$frequency
  return(structure(
    res,
    class = class(x),
    is_complex = attr(x, "is_complex"),
    tsp_orig = attr(x, "tsp_orig"),
    repr = repr
  ))
}

# Declare .data as a global variable to avoid R CMD check NOTE
utils::globalVariables(c(".data"))

#' Plot the modulus of FFT results
#'
#' Plots the modulus of the FFT results against the frequencies.
#'
#' @param x A `tidy_fft` object.
#' @param ... passed to ggplot.
#' @exportS3Method graphics::plot
plot.tidy_fft <- function(x, ...) {
  change_repr(x, "polr") |>
    ggplot2::ggplot(...) +
    ggplot2::aes(x = .data$frequency, y = .data$mod) +
    ggplot2::geom_line() +
    ggplot2::ylab("modulus") +
    ggplot2::theme_classic()
}
