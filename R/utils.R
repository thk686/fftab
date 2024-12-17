# Declare .data as a global variable to avoid R CMD check NOTE
utils::globalVariables(c(".data", "arg", "dim_1", "fx", "im", "mod", "re"))

#' Convert an Object to a Tidy FFT Object
#'
#' This internal helper function applies the necessary structure and class
#' attributes to convert a given object into a `tidy_fft` object.
#'
#' A `tidy_fft` object is a tibble-like structure with additional class
#' attributes used to represent Fourier Transform results in a tidy format.
#'
#' @param x The input object to convert, typically a tibble or data frame.
#' @param ... Additional attributes to include in the structured object, such
#'   as metadata or specific attributes required for Fourier Transform analysis.
#' @return The input object `x`, with the `tidy_fft` class and any additional
#'   attributes provided in `...`.
#' @keywords internal
.as_tidy_fft_obj <- function(x, ...) {
  structure(x, ..., class = c("tidy_fft", "tbl_df", "tbl", "data.frame"))
}

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
  if (n < 2) {
    stop("Minimum length is 2")
  }
  k <- 0:(n - 1)
  ifelse(k <= n / 2, k, k - n) / n
}

#' Check if a `tidy_fft` object has a specific representation
#'
#' This function checks if the given `tidy_fft` object contains the specified
#' representation: polar (`'polr'`), rectangular (`'rect'`), or complex (`'cplx'`).
#'
#' @param x A `tidy_fft` object.
#' @param repr The target representation to check (`'polr'`, `'rect'`, or `'cplx'`).
#' @return A logical value (`TRUE` if the object has the specified representation, otherwise `FALSE`).
#' @examples
#' tft <- tidy_fft(c(1, 0, -1, 0))
#' can_repr(tft, "cplx") # Returns TRUE
#' can_repr(tft, "rect") # Returns FALSE
#' @export
can_repr <- function(x, repr = c('cplx', 'rect', 'polr')) {
  switch(match.arg(repr),
         cplx = has_cplx(x),
         rect = has_rect(x),
         polr = has_polr(x))
}

#' Retrieve the current representation of a tidy_fft object
#'
#' @param x A `tidy_fft` object.
#' @return A vector of possible representations
#' @export
get_repr <- function(x) {
  i <- c(has_cplx(x), has_rect(x), has_polr(x))
  c('cplx', 'rect', 'polr')[i]
}

has_cplx <- function(x) {
  "fx" %in% names(x)
}

has_rect <- function(x) {
  all(c("re", "im") %in% names(x))
}

has_polr <- function(x) {
  all(c("mod", "arg") %in% names(x))
}

get_repr <- function(x) {
  c("cplx", "rect", "polr")[c(has_cplx(x), has_rect(x), has_polr(x))]
}

to_cplx <- function(x, .keep = "unused") {
  if (has_cplx(x)) {
    return(x)
  }
  if (has_rect(x)) {
    dplyr::mutate(x, fx = complex(real = re, imaginary = im), .keep = .keep)
  } else {
    dplyr::mutate(x,
                  fx = complex(modulus = mod, argument = arg),
                  .keep = .keep)
  }
}

to_rect <- function(x, .keep = "unused") {
  if (has_rect(x)) {
    return(x)
  }
  if (has_cplx(x)) {
    dplyr::mutate(x,
                  re = Re(fx),
                  im = Im(fx),
                  .keep = .keep)
  } else {
    dplyr::mutate(x,
                  re = mod * cos(arg),
                  im = mod * sin(arg),
                  .keep = .keep)
  }
}

to_polr <- function(x, .keep = "unused") {
  if (has_polr(x)) {
    return(x)
  }
  if (has_cplx(x)) {
    dplyr::mutate(x,
                  mod = Mod(fx),
                  arg = Arg(fx),
                  .keep = .keep)
  } else {
    dplyr::mutate(
      x,
      mod = sqrt(re ^ 2 + im ^ 2),
      arg = atan2(im, re),
      .keep = .keep
    )
  }
}

#' Extract Fourier Coefficients and Derived Components
#'
#' These utility functions convert a `tidy_fft` object to the desired representation
#' (`'cplx'`, `'rect'`, or `'polr'`) and extract specific components.
#'
#' - `get_fx()`: Extracts the complex Fourier coefficients (`fx`) from the `'cplx'` representation.
#' - `get_re()`: Extracts the real part (`re`) of the Fourier coefficients from the `'rect'` representation.
#' - `get_im()`: Extracts the imaginary part (`im`) of the Fourier coefficients from the `'rect'` representation.
#' - `get_mod()`: Extracts the magnitude (`mod`) of the Fourier coefficients from the `'polr'` representation.
#' - `get_arg()`: Extracts the phase angle (`arg`) of the Fourier coefficients from the `'polr'` representation.
#'
#' @param x A `tidy_fft` object containing FFT results in any representation.
#' @return Each function returns the requested component:
#' - `get_fx()`: A complex vector of Fourier coefficients.
#' - `get_re()`: A numeric vector of real parts.
#' - `get_im()`: A numeric vector of imaginary parts.
#' - `get_mod()`: A numeric vector of magnitudes.
#' - `get_arg()`: A numeric vector of phase angles (in radians).
#' @examples
#' # Example usage
#' fft_result <- tidy_fft(c(1, 0, -1, 0))
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
#' @export
get_fx <- function(x) {
  to_cplx(x, .keep = "none")$fx
}

#' @rdname get_fx
#' @export
get_rect <- function(x) {
  to_rect(x, .keep = "none")
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
get_polr <- function(x) {
  to_polr(x, .keep = "none")
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

drop_negf <- function(x) {
  dplyr::filter(x, dplyr::if_all(dplyr::starts_with("dim_"), ~ . >= 0))
}

drop_posf <- function(x) {
  dplyr::filter(x, dplyr::if_any(dplyr::starts_with("dim_"), ~ . < 0))
}

# Generate negative frequency components from positive frequencies
generate_negf <- function(x) {
  # Drop nyquist frequency if odd number of rows
  if (nrow(x) %% 2 == 1) {
    x <- dplyr::filter(x, dplyr::if_any(dplyr::starts_with("dim_"), ~ . != max(dim_1)))
  }
  dplyr::filter(x, dplyr::if_any(dplyr::starts_with("dim_"), ~ . != 0)) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("dim_"), ~ -.), fx = Conj(fx))
}

#' Plot the modulus of FFT results
#'
#' Plots the modulus of the FFT results against the frequencies.
#'
#' @param x A `tidy_fft` object.
#' @param ... passed to ggplot.
#' @exportS3Method graphics::plot
plot.tidy_fft <- function(x, ...) {
  to_polr(x) |>
    ggplot2::ggplot(...) +
    ggplot2::aes(x = .data$dim_1, y = .data$mod) +
    ggplot2::geom_line() +
    ggplot2::ylab("modulus") +
    ggplot2::theme_classic()
}
