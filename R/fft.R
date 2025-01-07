#' Perform FFT and IFFT with Tidy Results
#'
#' Provides functions to compute the Fast Fourier Transform (FFT) and its inverse (IFFT)
#' while maintaining results in a tabular format. Supports vectors, time series
#' (`ts`), and arrays as inputs.
#'
#' @section FFT:
#' The `fftab` function computes the FFT for different input types:
#'
#' - **Default Input (`fftab.default`)**: Computes FFT for numeric vectors.
#' - **Time Series Input (`fftab.ts`)**: Handles FFT for `ts` objects, scaling frequencies appropriately.
#' - **Array Input (`fftab.array`)**: Processes multidimensional arrays.
#'
#' Results are returned as a tibble containing Fourier frequencies and FFT values.
#'
#' @section IFFT:
#' The `ifftab` function reconstructs the original signal from a `fftab` object.
#' It supports vectors, arrays, and time series inputs. The inverse transform preserves
#' the original structure (e.g., array dimensions or time series attributes).
#'
#' @param x Input object for which to compute the FFT or IFFT. This can be:
#'   - A numeric vector (default method for `fftab`).
#'   - A time series object (`ts`) for `fftab.ts`.
#'   - A multidimensional numeric array for `fftab.array`.
#'   - A `fftab` object for `ifftab`.
#' @param norm Logical. If `TRUE`, computes normalized coefficients for FFT.
#'
#' @return
#' - **`fftab`**: A tibble containing:
#'   - Fourier frequencies (`.dim_1`, `.dim_2`, etc.).
#'   - FFT values stored in the `fx` column as complex values.
#' - **`ifftab`**: A vector, array, or time series object representing the reconstructed signal.
#'
#' @details
#' - `fftab` organizes FFT results into a tibble for downstream analysis.
#' - `ifftab` ensures that reconstructed signals match the input structure (e.g., arrays, `ts`).
#'
#' @examples
#' fftab(c(1, 0, -1, 0))
#'
#' fftab(c(1, 0, -1, 0)) |> ifftab()
#'
#' ts(sin(1:10), frequency = 12) |> fftab()
#'
#' array(1:8, dim = c(2, 2, 2)) |> fftab()
#'
#' @seealso [stats::fft()]
#' @export
fftab <- function(x, norm = FALSE) {
  UseMethod("fftab")
}

#' @rdname fftab
#' @export
fftab.default <- function(x, norm = FALSE) {
  stopifnot(is.numeric(x) || is.complex(x), length(x) > 0)
  fourier_frequencies(x) |>
    tibble::add_column(fx = .fft(x, norm)) |>
    .as_fftab_obj(
      .is_normalized = norm,
      .is_complex = is.complex(x)
    )
}

#' @rdname fftab
#' @export
fftab.ts <- function(x, norm = FALSE) {
  fftab(as.vector(x), norm = norm) |>
    dplyr::mutate(.dim_1 = .dim_1 * frequency(x)) |>
    structure(.tsp = attr(x, "tsp"))
}

#' @rdname fftab
#' @export
fftab.array <- function(x, norm = FALSE) {
  fourier_frequencies(x) |>
    dplyr::mutate(fx = as.vector(.fft(x, norm))) |>
    .as_fftab_obj(
      .is_normalized = norm,
      .is_complex = is.complex(x),
      .dim = dim(x)
    )
}

#' @rdname fftab
#' @export
ifftab <- function(x) {
  stopifnot(inherits(x, "fftab"))
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
