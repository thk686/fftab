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
#'   - `.dim_1`, `.dim_2`, ..., represent the Fourier frequencies for each dimension.
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
#' ts(rnorm(36), frequency = 12) |> fourier_frequencies()
#'
#' # Multidimensional array input
#' array(1:27, dim = c(3, 3, 3)) |> fourier_frequencies()
#'
#' # Matrix input
#' matrix(1:9, nrow = 3, ncol = 3) |> fourier_frequencies()
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
  tibble::tibble(.dim_1 = .fourier_frequencies(x))
}

# Time series method
#' @rdname fourier_frequencies
#' @export
fourier_frequencies.ts <- function(x) {
  tibble::tibble(.dim_1 = .fourier_frequencies(x) * frequency(x))
}

# Array method
#' @rdname fourier_frequencies
#' @export
fourier_frequencies.array <- function(x) {
  ff <- list()
  dims <- dim(x)
  for (i in rev(seq_along(dims))) {
    ff[[paste0(".dim_", i)]] <- .fourier_frequencies(dims[i])
  }
  rev(tidyr::expand_grid(!!!ff))
}

#' Add L2 Norm and Squared L2 Norm of Dimensions
#'
#' These functions compute and add L2 norms and squared L2 norms of the frequency dimensions
#' (`.dim_*` columns) as additional columns to a `tidy_fft` object.
#'
#' @param x A `tidy_fft` object containing frequency dimensions (`.dim_*`) and associated metadata.
#'
#' @return A `tidy_fft` object with additional columns:
#' - **`l2nm`**: The L2 norm of the frequency dimensions.
#' - **`l2sq`**: The squared L2 norm of the frequency dimensions.
#'
#' @details
#' - **`add_l2nm()`**: Adds a column `l2nm` containing the L2 norm, calculated as the square root
#'   of the sum of squares across `.dim_*` columns.
#' - **`add_l2sq()`**: Adds a column `l2sq` containing the squared L2 norm, calculated as the sum
#'   of squares across `.dim_*` columns.
#'
#' @seealso
#' - [tidy_fft()]
#' - [tibble::add_column()]
#'
#' @examples
#' matrix(1:9, 3) |>
#'   tidy_fft() |>
#'   print(n = 3) |>
#'   add_l2nm() |>
#'   print(n = 3) |>
#'   add_l2sq() |>
#'   print(n = 3)
#'
#' @export
add_l2nm <- function(x) {
  tibble::add_column(x, l2nm = .get_dim_cols(x)^2 |> rowSums() |> sqrt())
}

#' @rdname add_l2nm
#' @export
add_l2sq <- function(x) {
  tibble::add_column(x, l2sq = .get_dim_cols(x)^2 |> rowSums())
}

