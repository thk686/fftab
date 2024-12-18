# Declare .data as a global variable to avoid R CMD check NOTE
utils::globalVariables(c(".data", "arg", "dim_1", "fx", "im", "mod", "re", "frequency"))

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

#' Check and Retrieve Representations of a `tidy_fft` Object
#'
#' These functions check and retrieve specific representations of a `tidy_fft` object.
#' Supported representations include:
#' - **Complex (`"cplx"`)**: Contains a column `fx` with complex Fourier coefficients.
#' - **Rectangular (`"rect"`)**: Contains columns `re` (real) and `im` (imaginary) components.
#' - **Polar (`"polr"`)**: Contains columns `mod` (modulus) and `arg` (argument).
#'
#' @param x A `tidy_fft` object.
#' @param repr For `can_repr()`, the target representation to check. A character string (`"polr"`, `"rect"`, or `"cplx"`).
#'
#' @return
#' - **`can_repr()`**: A logical value (`TRUE` or `FALSE`) indicating if the object has the specified representation.
#' - **`get_repr()`**: A character vector of representations present in the object.
#' - **`has_cplx()`, `has_rect()`, `has_polr()`**: Logical values (`TRUE` or `FALSE`) indicating the presence of specific representations.
#'
#' @examples
#' tft <- tidy_fft(c(1, 0, -1, 0))
#'
#' # Check specific representations
#' can_repr(tft, "cplx") # TRUE
#' can_repr(tft, "rect") # FALSE
#'
#' # Retrieve current representations
#' get_repr(tft) # "cplx"
#'
#' # Check individual representations
#' has_cplx(tft) # TRUE
#' has_rect(tft) # FALSE
#' has_polr(tft) # FALSE
#'
#' @export
can_repr <- function(x, repr) {
  res <- 0
  for (r in repr) {
    res <- res + switch(r,
                        cplx = has_cplx(x),
                        rect = has_rect(x),
                        polr = has_polr(x))
  }
  res > 0
}

#' @rdname can_repr
#' @export
get_repr <- function(x) {
  c("cplx", "rect", "polr")[c(has_cplx(x), has_rect(x), has_polr(x))]
}

#' @rdname can_repr
#' @export
has_cplx <- function(x) {
  "fx" %in% names(x)
}

#' @rdname can_repr
#' @export
has_rect <- function(x) {
  all(c("re", "im") %in% names(x))
}

#' @rdname can_repr
#' @export
has_polr <- function(x) {
  all(c("mod", "arg") %in% names(x))
}

#' Convert a `tidy_fft` Object Between Representations
#'
#' These functions convert a `tidy_fft` object to a specified representation:
#' - **`to_cplx()`**: Converts to complex representation (`fx`).
#' - **`to_rect()`**: Converts to rectangular representation (`re`, `im`).
#' - **`to_polr()`**: Converts to polar representation (`mod`, `arg`).
#'
#' @param x A `tidy_fft` object.
#' @param .keep Specifies which columns to retain. Defaults to `"unused"`.
#'
#' @return A modified `tidy_fft` object containing the specified representation:
#' - **`to_cplx()`**: Adds the `fx` column for complex values.
#' - **`to_rect()`**: Adds the `re` and `im` columns for rectangular components.
#' - **`to_polr()`**: Adds the `mod` and `arg` columns for polar components.
#'
#' @details
#' - **`to_cplx()`**: Converts from rectangular (`re`, `im`) or polar (`mod`, `arg`) components to complex form.
#' - **`to_rect()`**: Converts from complex (`fx`) or polar components to rectangular form.
#' - **`to_polr()`**: Converts from complex (`fx`) or rectangular components to polar form.
#'
#' @examples
#' tft <- tidy_fft(c(1, 0, -1, 0))
#'
#' # Convert to different representations
#' tft_cplx <- to_cplx(tft) # Complex representation
#' tft_rect <- to_rect(tft_cplx) # Rectangular representation
#' tft_polr <- to_polr(tft_cplx) # Polar representation
#'
#' # Print results
#' print(tft_cplx)
#' print(tft_rect)
#' print(tft_polr)
#'
#' @seealso [can_repr()], [get_repr()]
#' @export
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

#' @rdname to_cplx
#' @export
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

#' @rdname to_cplx
#' @export
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
#' These utility functions extract or convert a `tidy_fft` object to the desired representation
#' (`'cplx'`, `'rect'`, or `'polr'`) and extract specific components.
#'
#' @param x A `tidy_fft` object containing FFT results in any representation.
#'
#' @return The requested components or converted representations:
#' - **`get_fx()`**: A complex vector of Fourier coefficients (`fx`).
#' - **`get_rect()`**: A tibble containing rectangular representation (`re`, `im`).
#' - **`get_re()`**: A numeric vector of real parts (`re`).
#' - **`get_im()`**: A numeric vector of imaginary parts (`im`).
#' - **`get_polr()`**: A tibble containing polar representation (`mod`, `arg`).
#' - **`get_mod()`**: A numeric vector of magnitudes (`mod`).
#' - **`get_arg()`**: A numeric vector of phase angles (`arg`), in radians.
#'
#' @details
#' - **`get_fx()`**: Extracts the complex Fourier coefficients from the `'cplx'` representation.
#' - **`get_rect()`**: Converts to the rectangular form and returns a tibble with `re` (real) and `im` (imaginary) components.
#' - **`get_re()`**: Extracts the real part from the rectangular representation.
#' - **`get_im()`**: Extracts the imaginary part from the rectangular representation.
#' - **`get_polr()`**: Converts to the polar form and returns a tibble with `mod` (magnitude) and `arg` (phase angle).
#' - **`get_mod()`**: Extracts the magnitude from the polar representation.
#' - **`get_arg()`**: Extracts the phase angle (in radians) from the polar representation.
#'
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
#' rect_values <- get_rect(fft_result)
#' polr_values <- get_polr(fft_result)
#'
#' print(fx_values)
#' print(re_values)
#' print(im_values)
#' print(mod_values)
#' print(arg_values)
#' print(rect_values)
#' print(polr_values)
#'
#' @seealso [to_cplx()], [to_rect()], [to_polr()]
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
