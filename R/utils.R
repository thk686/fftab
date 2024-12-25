# Declare .data as a global variable to avoid R CMD check NOTE
utils::globalVariables(c(".data", "arg", "dim_1", "fx", "im", "mod", "re", "frequency"))

#' Compute the Fast Fourier Transform (FFT) of a Vector
#'
#' This function computes the Fast Fourier Transform (FFT) of a numeric vector.
#' Optionally, it normalizes the result by dividing it by the length of the input vector.
#'
#' @param x A numeric vector representing the input signal to transform.
#' @param norm A logical value indicating whether to normalize the FFT output
#'             by dividing it by the length of the input vector. Default is `FALSE`.
#'
#' @return A complex vector representing the FFT of the input signal.
#'
#' @details
#' The function wraps around the base R `stats::fft` function and provides an
#' option for normalization.
#'
#' @examples
#' x <- c(1, 2, 3, 4)
#' .fft(x)
#' .fft(x, norm = TRUE)
#'
#' @keywords internal
#' @importFrom stats fft
.fft <- function(x, norm = FALSE) {
  if (norm) {
    stats::fft(x) / length(x)
  } else {
    stats::fft(x)
  }
}

#' @keywords internal
.get_dim_cols <- function(x) {
  dplyr::select(x, dplyr::starts_with("dim_"))
}

#' @keywords internal
.is_normalized <- function(x) {
  attr(x, ".is_normalized")
}

#' Check if an object is complex
#'
#' Retrieves the `.is_complex` attribute of the given object.
#'
#' @param x An object to check for the `.is_complex` attribute.
#' @return The value of the `.is_complex` attribute, or `NULL` if not present.
#' @keywords internal
.is_complex <- function(x) {
  attr(x, ".is_complex")
}

#' Retrieve the dimensions of an object
#'
#' Retrieves the `.dim` attribute of the given object.
#'
#' @param x An object to check for the `.dim` attribute.
#' @return The value of the `.dim` attribute, or `NULL` if not present.
#' @keywords internal
.dim <- function(x) {
  attr(x, ".dim")
}

#' Check if an object is an array
#'
#' Determines if the object has a `.dim` attribute, indicating it is an array.
#'
#' @param x An object to check.
#' @return `TRUE` if the object has a `.dim` attribute, `FALSE` otherwise.
#' @keywords internal
.is_array <- function(x) {
  !is.null(attr(x, ".dim"))
}

#' Retrieve the size of an object
#'
#' Retrieves the `.size` attribute of the given object.
#'
#' @param x An object to check for the `.size` attribute.
#' @return The value of the `.size` attribute, or `NULL` if not present.
#' @keywords internal
.size <- function(x) {
  attr(x, ".size")
}

#' Retrieve the time series parameters of an object
#'
#' Retrieves the `.tsp` attribute of the given object.
#'
#' @param x An object to check for the `.tsp` attribute.
#' @return The value of the `.tsp` attribute, or `NULL` if not present.
#' @keywords internal
.tsp <- function(x) {
  attr(x, ".tsp")
}

#' Check if an object is a time series
#'
#' Determines if the object has a `.tsp` attribute, indicating it is a time series.
#'
#' @param x An object to check.
#' @return `TRUE` if the object has a `.tsp` attribute, `FALSE` otherwise.
#' @keywords internal
.is_ts <- function(x) {
  !is.null(.tsp(x))
}

#' @keywords internal
.frequency <- function(x) {
  if (.is_ts(x)) {
    .tsp(x)[3]
    } else {
      1
    }
}

#' @keywords internal
.nyquist <- function(x) {
  .frequency(x) / 2
}

#' Convert an Object to a Tidy FFT Object
#'
#' This internal helper function applies the necessary structure and class
#' attributes to convert a given object into a `tidy_fft` object.
#'
#' A `tidy_fft` object is a tibble-like structure with additional class
#' attributes used to represent Fourier Transform results in a tidy format.
#'
#' @param x The input object to convert, typically a tibble or data frame.
#' @param .norm Is the input normalized?
#' @param .complex Was the original data complex?
#' @param ... Additional attributes to include in the structured object, such as
#'   metadata or specific attributes required for Fourier Transform analysis.
#' @return The input object `x`, with the `tidy_fft` class and any additional
#'   attributes provided in `...`.
#' @keywords internal
.as_tidy_fft_obj <- function(x, .is_normalized, .is_complex, ...) {
  structure(x,
    ...,
    .size = nrow(x),
    .is_complex = .is_complex,
    .is_normalized = .is_normalized,
    class = c("tidy_fft", class(x))
  )
}

#' @keywords internal
.set_repr <- function(x, repr, .keep = "unused") {
  if (repr == "cplx") {
    return(to_cplx(x, .keep = .keep))
  }
  if (repr == "rect") {
    return(to_rect(x, .keep = .keep))
  }
  if (repr == "polr") {
    return(to_polr(x, .keep = .keep))
  }
  stop("Invalid representation.")
}

#' @keywords internal
.shift_phase <- function(x, shift) {
  if (.is_array(x)) {
    .NotYetImplemented()
  } else {
    if (.is_complex(x)) {
      to_polr(x) |> dplyr::mutate(arg = arg + shift)
    } else {
      to_polr(x) |>
        dplyr::mutate(
          arg = dplyr::case_when(
            dim_1 == 0.0 ~ arg,
            dim_1 == .nyquist(x) ~ arg,
            dim_1 > 0.0 ~ arg + shift,
            dim_1 < 0.0 ~ arg - shift,
            TRUE ~ arg
          )
        )
    }
  }
}

#' @keywords internal
.variance <- function(x) {
  to_polr(x) |>
    remove_dc() |>
    dplyr::mutate(mod = mod^2) |>
    get_mod() |>
    sum() / .size(x)^2
}

# .gen_indices <- function(dims) {
#   ff <- lapply(dims, seq_len) # Generate sequences for each dimension
#   do.call(expand.grid, ff) # Expand the grid in row-major order
# }
#
# # Function to compute conjugate indices with modulo
# find_cc_2d <- function(i, j, ni, nj) {
#   c((ni - i + 1) %% ni + 1, (nj - j + 1) %% nj + 1)
# }
#
# # Function to check if an entry is unique
# is_unique <- function(i, j, ni, nj) {
#   # Row-major index of current entry
#   k <- (i - 1) * nj + j
#   # Row-major index of its conjugate
#   cc <- find_cc_2d(i, j, ni, nj)
#   kc <- (cc[1] - 1) * nj + cc[2]
#   # Keep if kc >= k
#   kc >= k
# }
#
# # Filter unique indices in a matrix of size ni x nj
# filter_unique <- function(ni, nj) {
#   # Generate all indices
#   indices <- .gen_indices(c(4, 4))
#
#   # Apply uniqueness check
#   unique_indices <- indices[apply(indices, 1, function(row) {
#     is_unique(row[1], row[2], ni, nj)
#   }), ]
#
#   unique_indices
# }

