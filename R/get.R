#' Extract Fourier Coefficients and Components
#'
#' These utility functions extract specific components from a `tidy_fft` object.
#'
#' @param x A `tidy_fft` object containing FFT results.
#'
#' @return The requested components:
#' - **`get_fx`**: A complex vector of Fourier coefficients (`fx`).
#' - **`get_re`**: A numeric vector of real parts (`re`).
#' - **`get_im`**: A numeric vector of imaginary parts (`im`).
#' - **`get_mod`**: A numeric vector of magnitudes (`mod`).
#' - **`get_arg`**: A numeric vector of phase angles (`arg`), in radians.
#'
#' @examples
#' fft_result <- tidy_fft(c(1, 0, -1, 0)) |>
#'   to_rect(.keep = "all") |>
#'   to_polr(.keep = "all")
#'
#' # Extract components
#' fx_values <- get_fx(fft_result)
#' re_values <- get_re(fft_result)
#' im_values <- get_im(fft_result)
#' mod_values <- get_mod(fft_result)
#' arg_values <- get_arg(fft_result)
#'
#' print(fx_values)
#' print(re_values)
#' print(im_values)
#' print(mod_values)
#' print(arg_values)
#'
#' @seealso [to_cplx()]
#' @export
get_fx <- function(x) {
  to_cplx(x, .keep = "none")$fx
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

#' @rdname get_fx
#' @export
get_fx_norm <- function(x) {
  if (.is_normalized(x)) {
    get_fx(x)
  } else {
    get_fx(x) / .size(x)
  }
}


#' Extract Rectangular or Polar Components
#'
#' The `get_rect` and `get_polr` functions extract specific components from a
#' `tidy_fft` object, representing the Fourier coefficients in either rectangular
#' or polar form.
#'
#' @param x A `tidy_fft` object containing FFT results.
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
#' fft_result <- tidy_fft(c(1, 0, -1, 0))
#'
#' # Extract rectangular components
#' rect_values <- get_rect(fft_result)
#' print(rect_values)
#'
#' # Extract polar components
#' polr_values <- get_polr(fft_result)
#' print(polr_values)
#'
#' @seealso [get_fx()], [get_re()], [get_mod()], [to_rect()], [to_polr()]
#' @export
get_rect <- function(x) {
  to_rect(x, .keep = "none") |> as.matrix()
}

#' @rdname get_rect
#' @export
get_polr <- function(x) {
  to_polr(x, .keep = "none") |> as.matrix()
}
