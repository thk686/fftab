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
#' @examples
#' # Example usage
#' library(tibble)
#' fft_result <- tibble(dim_1 = c(0, 0.25, 0.5, -0.25), fx = complex(real = c(1, 0, -1, 0)))
#' tidy_fft_obj <- .as_tidy_fft_obj(fft_result, is_complex = TRUE)
#' class(tidy_fft_obj) # Returns c("tidy_fft", "tbl_df", "tbl", "data.frame")
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
  if (input_len == 0)
    stop("Input must not be empty")
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
  if (n < 2)
    stop("Minimum length is 2")
  k <- 0:(n - 1)
  ifelse(k <= n / 2, k, k - n) / n
}

#' Check if a `tidy_fft` object has a specific representation
#'
#' This function checks if the given `tidy_fft` object contains the specified
#' representation: polar (`"polr"`), rectangular (`"rect"`), or complex (`"cplx"`).
#'
#' @param x A `tidy_fft` object.
#' @param repr The target representation to check (`"polr"`, `"rect"`, or `"cplx"`).
#' @return A logical value (`TRUE` if the object has the specified representation, otherwise `FALSE`).
#' @examples
#' tft <- tidy_fft(c(1, 0, -1, 0), repr = "cplx")
#' can_repr(tft, "cplx")  # Returns TRUE
#' can_repr(tft, "rect")  # Returns FALSE
#' @export
can_repr <- function(x, repr = c("polr", "rect", "cplx")) {
  switch(
    match.arg(repr),
    polr = all(c("mod", "arg") %in% names(x)),
    rect = all(c("re", "im") %in% names(x)),
    cplx = "fx" %in% names(x)
  )
}

#' Retrieve the current representation of a tidy_fft object
#'
#' @param x A `tidy_fft` object.
#' @return A vector of possible representations
#' @export
get_repr <- function(x) {
  reprs <- c("polr", "rect", "cplx")
  reprs[sapply(reprs, function(repr)
    can_repr(x, repr))]
}

#' Change the representation of FFT results
#'
#' Converts FFT results between different representations.
#'
#' @param x A `tidy_fft` object.
#' @param repr The target representation (`"polr"`, `"rect"`, or `"cplx"`).
#' @param .keep Determines which columns to retain in the output.
#'   Passed directly to \code{dplyr::mutate}. Options include:
#'   \itemize{
#'     \item \code{"all"}: Keep all columns (default).
#'     \item \code{"used"}: Keep only columns used in the computation, plus the new ones.
#'     \item \code{"unused"}: Keep only columns not used in the computation, plus the new ones.
#'     \item \code{"none"}: Keep only the new columns.
#'   }
#' @return The modified object with the new representation.
#' @examples
#' fft_res <- tidy_fft(c(1, 0, -1, 0), repr = "cplx")
#' change_repr(fft_res, "rect", .keep = "none")
#' @export
change_repr <- function(x,
                        repr = c("polr", "rect", "cplx"),
                        .keep = "unused") {
  to <- match.arg(repr)
  if (can_repr(x, "cplx")) {
    switch(to,
           polr = return(dplyr::mutate(
             x,
             mod = Mod(fx),
             arg = Arg(fx),
             .keep = .keep
           )),
           rect = return(dplyr::mutate(
             x,
             re = Re(fx),
             im = Im(fx),
             .keep = .keep
           )),
           cplx = return(x))
  }
  if (can_repr(x, "rect")) {
    switch(to,
           polr = return(dplyr::mutate(
             x,
             mod = Mod(complex(
               real = re, imaginary = im
             )),
             arg = Arg(complex(
               real = re, imaginary = im
             )),
             .keep = .keep
           )),
           rect = return(x),
           cmplx = return(dplyr::mutate(
             x, fx = complex(real = re, imaginary = im), .keep = .keep
           )))
  }
  if (can_repr(x, "polr")) {
    switch(to,
           polr = return(x),
           rect = return(dplyr::mutate(
             x,
             re = Re(complex(
               argument = arg, modulus = mod
             )),
             im = Im(complex(
               argument = arg, modulus = mod
             )),
             .keep = .keep
           )),
           cplx = return(dplyr::mutate(
             x,
             fx = complex(argument = arg, modulus = mod),
             .keep = .keep
           )))
  }
  stop("Cannot convert to ", to, " representation.")
}

#' Plot the modulus of FFT results
#'
#' Plots the modulus of the FFT results against the frequencies.
#'
#' @param x A `tidy_fft` object.
#' @param ... passed to ggplot.
#' @exportS3Method graphics::plot
plot.tidy_fft <- function(x, ...) {
  change_repr(x, "polr") |>
    ggplot2::ggplot(...) +
    ggplot2::aes(x = .data$dim_1, y = .data$mod) +
    ggplot2::geom_line() +
    ggplot2::ylab("modulus") +
    ggplot2::theme_classic()
}
