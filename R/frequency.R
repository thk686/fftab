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
#' fourier_frequencies(array_input)
#'
#' # Compute Fourier frequencies for a 2D matrix
#' matrix_input <- matrix(1:9, nrow = 3, ncol = 3)
#' fourier_frequencies(matrix_input)
#'
#' # Compute Fourier frequencies for a vector (1D case)
#' vector_input <- 1:8
#' fourier_frequencies(vector_input)
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

