#' Compute the number of samples in an input
#'
#' This helper function determines the number of samples in the input object.
#' For a vector, it returns its length. For a matrix or data frame, it returns
#' the number of rows.
#'
#' @param x An input object (scalar, vector, matrix, or data frame).
#' @return An integer representing the number of samples (rows) in the input object.
#' @keywords internal
.num_samples <- function(x) {
  input_len <- length(unlist(x))
  if (input_len == 0)
    stop("Input must not be empty")
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
  if (n < 2)
    stop("Minimum length is 2")
  k <- 0:(n - 1)
  ifelse(k <= n / 2, k, k - n) / n
}

#' Compute Fourier frequencies
#'
#' This is a generic function for computing Fourier frequencies for various
#' types of inputs, such as scalars, vectors, matrices, time series, or arrays.
#'
#' @param x The input object (scalar, vector, matrix, time series, or array).
#' @return A numeric vector for 1D inputs or a tibble for multidimensional inputs,
#'         containing the Fourier frequencies normalized or scaled based on the input type.
#' @export
fourier_frequencies <- function(x) {
  UseMethod("fourier_frequencies")
}

#' Compute Fourier frequencies for default inputs
#'
#' Computes normalized Fourier frequencies for scalar or vector inputs, which
#' represent evenly spaced values in the Fourier domain.
#'
#' @param x A scalar or vector representing the length of the sequence.
#' @return A tibble with a single column `dim_1` containing the normalized Fourier frequencies.
#' @examples
#' fourier_frequencies(8)
#' @export
fourier_frequencies.default <- function(x) {
  tibble::tibble(dim_1 = .fourier_frequencies(x))
}

#' Compute Fourier frequencies for time series objects
#'
#' Computes Fourier frequencies scaled by the frequency attribute of a `ts` object,
#' making them domain-specific to the time series.
#'
#' @param x A time series (`ts`) object.
#' @return A tibble with a single column `dim_1` containing the normalized Fourier frequencies.
#' @examples
#' ts_obj <- ts(rnorm(12), frequency = 12)
#' fourier_frequencies(ts_obj)
#' @export
fourier_frequencies.ts <- function(x) {
  n <- length(x)
  stride <- attr(x, "tsp")[3]
  fourier_frequencies(n) |>
           dplyr::mutate(dim_1 = stride * dim_1)
}

#' Compute Fourier frequencies for multidimensional arrays
#'
#' Computes Fourier frequencies for each dimension of an array and returns
#' a tidy tibble with frequencies for all dimensions.
#'
#' @param x A numeric array.
#' @return A tibble where each column corresponds to the Fourier frequencies
#'         of a dimension of the input array. The column names are `dim_1`,
#'         `dim_2`, ..., corresponding to each dimension.
#' @examples
#' # Compute Fourier frequencies for a 3D array
#' array_input <- array(1:27, dim = c(3, 3, 3))
#' fourier_frequencies.array(array_input)
#'
#' # Compute Fourier frequencies for a 2D matrix
#' matrix_input <- matrix(1:9, nrow = 3, ncol = 3)
#' fourier_frequencies.array(matrix_input)
#'
#' # Compute Fourier frequencies for a vector (1D case)
#' vector_input <- 1:8
#' fourier_frequencies.array(vector_input)
#' @seealso [fourier_frequencies.default()], [tidyr::expand_grid()]
#' @export
fourier_frequencies.array <- function(x) {
  ff <- list()
  dims <- dim(drop(x))
  for (i in rev(seq_along(dims))) {
    ff[[paste0("dim_", i)]] <- .fourier_frequencies(dims[i])
  }
  tidyr::expand_grid(!!!ff)
}

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
      is_complex = is.complex(x),
      repr = "cplx"
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
      dim_orig = dim(x),
      repr = "cplx"
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
  if (!is.null(attr(x, "tsp_o$ig"))) {
    attr(res, "tsp") <- attr(x, "tsp_orig")
    class(res) <- "ts"
  }
  res
}

#' Retrieve the current representation of a tidy_fft object
#'
#' @param x A `tidy_fft` object.
#' @return The current representation (`"polr"`, `"rect"`, or `"cplx"`).
#' @export
get_repr <- function(x) {
  attr(x, "repr")
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

  to <- match.arg(repr)
  from <- get_repr(x)

  if (from == to) {
    return(x)
  }

  if (from == "cplx" && to == "rect") {
    return(dplyr::mutate(
      structure(x, repr = to),
      re = Re(fx),
      im = Im(fx),
      .keep = "unused"
    ))
  }

  if (from == "cplx" && to == "polr") {
    return(dplyr::mutate(
      structure(x, repr = to),
      mod = Mod(fx),
      arg = Arg(fx),
      .keep = "unused"
    ))
  }

  if (from == "rect" && to == "cplx") {
    return(
      dplyr::mutate(structure(x, repr = to),
      fx = complex(real = re, imaginary = im),
      .keep = "unused")
    )
  }

  if (from == "rect" && to == "polr") {
    return(dplyr::mutate(
      structure(x, repr = to),
      mod = Mod(complex(
        real = re, imaginary = im
      )),
      arg = Arg(complex(
        real = re, imaginary = im
      )),
      .keep = "unused"
    ))
  }

  if (from == "polr" && to == "cplx") {
    return(dplyr::mutate(
      structure(x, repr = to),
      fx = complex(arg = arg, mod = mod),
      .keep = "unused"
    ))
  }

  if (from == "polr" && to == "rect") {
    return(dplyr::mutate(
      structure(x, repr = to),
      re = Re(complex(arg = arg, mod = mod)),
      im = Im(complex(arg = arg, mod = mod)),
      .keep = "unused"
    ))
  }

  stop("Cannot convert to", to)
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
    ggplot2::aes(x = .data$dim_1, y = .data$mod) +
    ggplot2::geom_line() +
    ggplot2::ylab("modulus") +
    ggplot2::theme_classic()
}
