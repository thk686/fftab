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

#' Remove DC Component and Symmetric Frequencies
#'
#' These functions operate on `tidy_fft` objects to manipulate and filter Fourier coefficients.
#'
#' @param x A `tidy_fft` object containing Fourier coefficients and associated metadata.
#'
#' @return A `tidy_fft` object with filtered coefficients.
#'
#' @details
#' - **`remove_dc()`**: Filters out rows where all `.dim_*` columns have a value of `0`.
#' - **`remove_symmetric()`**:
#'   - For real-valued signals, it filters out redundant, complex-conjugate frequencies
#'   - For complex-valued signals, no filtering is applied as symmetry isn't relevant.
#'
#' @seealso
#' - [dplyr::filter()]
#' - [tidy_fft()]
#'
#' @examples
#' matrix(rnorm(9), 3) |> tidy_fft() |> remove_dc()
#'
#' matrix(rnorm(9), 3) |> tidy_fft() |> remove_symmetric()
#'
#' @export
remove_dc <- function(x) {
  dplyr::filter(x, dplyr::if_any(dplyr::starts_with(".dim_"), ~ . != 0))
}

#' @rdname remove_dc
#' @export
remove_symmetric <- function(x) {
  if (.is_complex(x)) {
    return(x)
  }
  x <- .sort_dims(x)
  i <- .find_dc_row(x)
  dplyr::slice(x, i:nrow(x))
}

#' @rdname remove_dc
#' @export
split_symmetric <- function(x) {
  if (.is_complex(x)) {
    return(x)
  }
  x <- .sort_dims(x)
  i <- .find_dc_row(x)
  list(symmetric = dplyr::slice(x, 1:(i - 1)),
       asymmetric = dplyr::slice(x, i:nrow(x)))
}

