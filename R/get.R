#' Extract Fourier Coefficients and Components
#'
#' These utility functions extract specific components from a `fftab` object.
#' `get_fx` retrieves the raw Fourier coefficients, while `get_fx_norm` ensures the
#' coefficients are either normalized or not normalized based on the `norm` parameter.
#'
#' @param x A `fftab` object containing FFT results.
#' @param norm Logical. If `TRUE`, forces normalized coefficients. If `FALSE`,
#' ensures non-normalized coefficients.
#'
#' @return The requested components:
#' - **`get_fx`**: A complex vector of raw Fourier coefficients (`fx`) as stored in the object.
#' - **`get_fx_norm`**: A complex vector of Fourier coefficients, explicitly normalized
#' or non-normalized based on the `norm` parameter.
#' - **`get_re`**: A numeric vector of real parts (`re`).
#' - **`get_im`**: A numeric vector of imaginary parts (`im`).
#' - **`get_mod`**: A numeric vector of magnitudes (`mod`).
#' - **`get_arg`**: A numeric vector of phase angles (`arg`), in radians.
#'
#' @details
#' - **`get_fx`**: Returns coefficients as they are stored in the `fftab` object.
#' - **`get_fx_norm`**: Adjusts coefficients if they are not in the desired normalization state.
#' - **`get_re`, `get_im`**: Extract real and imaginary components.
#' - **`get_mod`, `get_arg`**: Compute magnitude and phase of coefficients.
#'
#' @examples
#' fftab(c(1, 0, -1, 0)) |> get_fx()
#'
#' fftab(c(1, 0, -1, 0)) |> get_fx_norm(norm = TRUE)
#'
#' fftab(c(1, 0, -1, 0)) |> get_re()
#'
#' fftab(c(1, 0, -1, 0)) |> get_im()
#'
#' fftab(c(1, 0, -1, 0)) |> get_mod()
#'
#' fftab(c(1, 0, -1, 0)) |> get_arg()
#'
#' @seealso [to_cplx()], [to_rect()], [to_polr()]
#'
#' @rdname get_fx
#' @export
get_fx <- function(x) {
  to_cplx(x, .keep = "none")$fx
}

#' @rdname get_fx
#' @export
get_fx_norm <- function(x, norm = FALSE) {
  if (norm) {
    if (.is_normalized(x)) {
      get_fx(x)
    } else {
      get_fx(x) / .size(x)
    }
  } else {
    if (.is_normalized(x)) {
      .size(x) * get_fx(x)
    } else {
      get_fx(x)
    }
  }
}

#' @rdname get_fx
#' @export
get_re <- function(x) {
  to_rect(x, .keep = "none")$re
}

#' @rdname get_fx
#' @export
get_im <- function(x) {
  to_rect(x, .keep = "none")$im
}

#' @rdname get_fx
#' @export
get_mod <- function(x) {
  to_polr(x, .keep = "none")$mod
}

#' @rdname get_fx
#' @export
get_arg <- function(x) {
  to_polr(x, .keep = "none")$arg
}

#' Extract Rectangular or Polar Components
#'
#' The `get_rect` and `get_polr` functions extract specific components from a
#' `fftab` object, representing the Fourier coefficients in either rectangular
#' or polar form.
#'
#' @param x A matrix object containing FFT results.
#'
#' @return
#' - **`get_rect`**: A matrix with two columns:
#'   - `re`: The real part of the coefficients.
#'   - `im`: The imaginary part of the coefficients.
#' - **`get_polr`**: A matrix with two columns:
#'   - `mod`: The magnitude of the coefficients.
#'   - `arg`: The phase angle of the coefficients, in radians.
#'
#' @examples
#' fftab(c(1, 0, -1, 0)) |> get_rect()
#'
#' fftab(c(1, 0, -1, 0)) |> get_polr()
#'
#' @seealso [get_fx()], [get_re()], [get_mod()], [to_rect()], [to_polr()]
#'
#' @rdname get_rect
#' @export
get_rect <- function(x) {
  to_rect(x, .keep = "none") |> as.matrix()
}

#' @rdname get_rect
#' @export
get_polr <- function(x) {
  to_polr(x, .keep = "none") |> as.matrix()
}
