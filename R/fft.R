#' Perform FFT and IFFT with Tidy Results
#'
#' Provides functions to compute the Fast Fourier Transform (FFT) and its inverse (IFFT)
#' while maintaining results in a tidy format using tibbles. Supports vectors, time series
#' (`ts`), and arrays as inputs.
#'
#' @section FFT:
#' The `tidy_fft` function computes the FFT for different input types:
#'
#' - **Default Input (`tidy_fft.default`)**: Computes FFT for numeric vectors.
#' - **Time Series Input (`tidy_fft.ts`)**: Handles FFT for `ts` objects, scaling frequencies appropriately.
#' - **Array Input (`tidy_fft.array`)**: Processes multidimensional arrays.
#'
#' Results are returned as tidy tibbles containing Fourier frequencies and FFT values.
#'
#' @section IFFT:
#' The `tidy_ifft` function reconstructs the original signal from a `tidy_fft` object.
#' It supports vectors, arrays, and time series inputs. The inverse transform preserves
#' the original structure (e.g., array dimensions or time series attributes).
#'
#' @param x Input object for which to compute the FFT or IFFT. This can be:
#'   - A numeric vector (default method for `tidy_fft`).
#'   - A time series object (`ts`) for `tidy_fft.ts`.
#'   - A multidimensional numeric array for `tidy_fft.array`.
#'   - A `tidy_fft` object for `tidy_ifft`.
#' @param norm Logical. If `TRUE`, computes normalized coefficients for FFT.
#'
#' @return
#' - **`tidy_fft`**: A tibble containing:
#'   - Fourier frequencies (`.dim_1`, `.dim_2`, etc.).
#'   - FFT values stored in the `fx` column as complex values.
#' - **`tidy_ifft`**: A vector, array, or time series object representing the reconstructed signal.
#'
#' @details
#' - `tidy_fft` organizes FFT results into a tidy tibble for downstream analysis.
#' - `tidy_ifft` ensures that reconstructed signals match the input structure (e.g., arrays, `ts`).
#'
#' @examples
#' tidy_fft(c(1, 0, -1, 0))
#'
#' tidy_fft(c(1, 0, -1, 0)) |> tidy_ifft()
#'
#' ts(sin(1:10), frequency = 12) |> tidy_fft()
#'
#' array(1:8, dim = c(2, 2, 2)) |> tidy_fft()
#'
#' @seealso [stats::fft()]
#' @export
tidy_fft <- function(x, norm = FALSE) {
  UseMethod("tidy_fft")
}

#' @rdname tidy_fft
#' @export
tidy_fft.default <- function(x, norm = FALSE) {
  stopifnot(is.numeric(x) || is.complex(x), length(x) > 0)
  fourier_frequencies(x) |>
    tibble::add_column(fx = .fft(x, norm)) |>
    .as_tidy_fft_obj(
      .is_normalized = norm,
      .is_complex = is.complex(x)
    )
}

#' @rdname tidy_fft
#' @export
tidy_fft.ts <- function(x, norm = FALSE) {
  tidy_fft(as.vector(x), norm = norm) |>
    dplyr::mutate(.dim_1 = .dim_1 * frequency(x)) |>
    structure(.tsp = attr(x, "tsp"))
}

#' @rdname tidy_fft
#' @export
tidy_fft.array <- function(x, norm = FALSE) {
  fourier_frequencies(x) |>
    dplyr::mutate(fx = as.vector(.fft(x, norm))) |>
    .as_tidy_fft_obj(
      .is_normalized = norm,
      .is_complex = is.complex(x),
      .dim = dim(x)
    )
}

#' @rdname tidy_fft
#' @export
tidy_ifft <- function(x) {
  stopifnot(inherits(x, "tidy_fft"))
  if (nrow(x) != .size(x)) {
    warning("Number of rows does not match the original object size.")
  }
  fx <- get_fx(x)
  if (!.is_normalized(x)) {
    fx <- fx / length(fx)
  }
  if (.is_array(x)) {
    dim(fx) <- .dim(x)
  }
  res <- stats::fft(fx, inverse = TRUE)
  if (!.is_complex(x)) {
    res <- Re(res)
  }
  if (.is_ts(x)) {
    attr(res, "tsp") <- .tsp(x)
    class(res) <- "ts"
  }
  res
}
