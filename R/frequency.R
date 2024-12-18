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
  if (n == 1) return(0)
  k <- 0:(n - 1)
  ifelse(k <= n / 2, k, k - n) / n
}

#' Compute Fourier Frequencies
#'
#' Computes Fourier frequencies for various types of inputs, such as scalars,
#' vectors, matrices, time series, or arrays. This generic function dispatches
#' appropriate methods based on the input type.
#'
#' @param x The input object. Supported input types:
#'   - **Scalar or vector**: The length of the sequence.
#'   - **Time series (`ts`)**: Frequencies are scaled based on the sampling rate.
#'   - **Multidimensional array or matrix**: Frequencies are computed for each dimension.
#'
#' @return A tibble where:
#'   - `dim_1`, `dim_2`, ..., represent the Fourier frequencies for each dimension.
#'
#' @details This function has the following methods:
#'
#' - **Default Input (`fourier_frequencies.default`)**: Computes normalized Fourier frequencies for scalar or vector inputs.
#' - **Time Series Input (`fourier_frequencies.ts`)**: Computes frequencies scaled by the frequency attribute of a `ts` object.
#' - **Multidimensional Arrays (`fourier_frequencies.array`)**: Computes frequencies for each dimension of a matrix or array.
#'
#' See the examples for details on each case.
#'
#' @examples
#' # Default input (vector)
#' fourier_frequencies(8)
#'
#' # Time series input
#' ts_obj <- ts(rnorm(36), frequency = 12)
#' fourier_frequencies(ts_obj)
#'
#' # Multidimensional array input
#' array_input <- array(1:27, dim = c(3, 3, 3))
#' fourier_frequencies(array_input)
#'
#' # Matrix input
#' matrix_input <- matrix(1:9, nrow = 3, ncol = 3)
#' fourier_frequencies(matrix_input)
#'
#' @seealso [tidyr::expand_grid()], [frequency()]
#'
#' @export
fourier_frequencies <- function(x) {
  UseMethod("fourier_frequencies")
}

# Default method
#' @rdname fourier_frequencies
#' @export
fourier_frequencies.default <- function(x) {
  tibble::tibble(dim_1 = .fourier_frequencies(x))
}

# Time series method
#' @rdname fourier_frequencies
#' @export
fourier_frequencies.ts <- function(x) {
  tibble::tibble(dim_1 = dim_1 * frequency(x))
}

# Array method
#' @rdname fourier_frequencies
#' @export
fourier_frequencies.array <- function(x) {
  ff <- list()
  dims <- dim(x)
  for (i in rev(seq_along(dims))) {
    ff[[paste0("dim_", i)]] <- .fourier_frequencies(dims[i])
  }
  rev(tidyr::expand_grid(!!!ff))
}
