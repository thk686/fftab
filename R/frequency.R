#' Compute Fourier Frequencies
#'
#' Computes Fourier frequencies for various types of inputs, such as scalars, vectors, matrices, time series, or arrays.
#' This generic function dispatches appropriate methods based on the input type.
#'
#' @param x The input object. Supported input types:
#'   - **Scalar or vector**: The length of the sequence.
#'   - **Time series (`ts`)**: Frequencies are scaled based on the sampling rate.
#'   - **Multidimensional array or matrix**: Frequencies are computed for each dimension.
#'
#' @return A numeric vector for 1D inputs or a tibble for multidimensional inputs, where:
#'   - `dim_1`, `dim_2`, ..., represent the Fourier frequencies for each dimension.
#'
#' @details
#' This function has the following methods:
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
  tibble::tibble(dim_1 = .fourier_frequencies(x) * frequency(x))
}

# Array method
#' @rdname fourier_frequencies
#' @export
fourier_frequencies.array <- function(x) {
  ff <- list()
  dims <- dim(drop(x))
  for (i in rev(seq_along(dims))) {
    ff[[paste0("dim_", i)]] <- .fourier_frequencies(dims[i])
  }
  tidyr::expand_grid(!!!ff)
}
