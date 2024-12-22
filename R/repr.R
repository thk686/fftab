#' Check and Retrieve Representations of a `tidy_fft` Object
#'
#' These functions check and retrieve specific representations of a `tidy_fft`
#' object. Supported representations include:
#' - **Complex (`"cplx"`)**: Contains a column `fx` with complex Fourier coefficients.
#' - **Rectangular (`"rect"`)**: Contains columns `re` (real) and `im` (imaginary) components.
#' - **Polar (`"polr"`)**: Contains columns `mod` (modulus) and `arg` (argument).
#'
#' @param x A `tidy_fft` object.
#' @param repr For `can_repr()`, the target representation to check. A character
#'   vector with one or more of (`"polr"`, `"rect"`, or `"cplx"`).
#'
#' @return
#' - **`can_repr()`**: A logical value (`TRUE` or `FALSE`) indicating if the object
#' has any of the specified representations.
#' - **`get_repr()`**: A character vector of representations present in the object.
#' - **`has_cplx()`, `has_rect()`, `has_polr()`**: Logical values (`TRUE` or `FALSE`)
#' indicating the presence of specific representations.
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
      polr = has_polr(x)
    )
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
  !is.null(x[["fx"]])
}

#' @rdname can_repr
#' @export
has_rect <- function(x) {
  !is.null(x[["re"]]) && !is.null(x[["im"]])
}

#' @rdname can_repr
#' @export
has_polr <- function(x) {
  !is.null(x[["mod"]]) && !is.null(x[["arg"]])
}

#' Convert a `tidy_fft` Object Between Representations
#'
#' These functions convert a `tidy_fft` object to a specified representation:
#' - **`to_cplx()`**: Converts to complex representation (`fx`).
#' - **`to_rect()`**: Converts to rectangular representation (`re`, `im`).
#' - **`to_polr()`**: Converts to polar representation (`mod`, `arg`).
#'
#' @param x A `tidy_fft` object.
#' @param .keep Specifies which columns to retain. See [mutate()].
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
    x
  } else {
    if (has_rect(x)) {
      dplyr::mutate(x,
        fx = complex(real = re, imaginary = im),
        .keep = .keep
      )
    } else {
      if (has_polr(x)) {
        dplyr::mutate(x,
          fx = complex(modulus = mod, argument = arg),
          .keep = .keep
        )
      } else {
        stop("No valid fft representation.")
      }
    }
  }
}

#' @rdname to_cplx
#' @export
to_rect <- function(x, .keep = "unused") {
  if (has_cplx(x)) {
    dplyr::mutate(x,
      re = Re(fx),
      im = Im(fx),
      .keep = .keep
    )
  } else {
    if (has_rect(x)) {
      x
    } else {
      if (has_polr(x)) {
        dplyr::mutate(
          x,
          re = mod * cos(arg),
          im = mod * sin(arg),
          .keep = .keep
        )
      } else {
        stop("No valid fft representation.")
      }
    }
  }
}

#' @rdname to_cplx
#' @export
to_polr <- function(x, .keep = "unused") {
  if (has_cplx(x)) {
    dplyr::mutate(x,
      mod = Mod(fx),
      arg = Arg(fx),
      .keep = .keep
    )
  } else {
    if (has_rect(x)) {
      dplyr::mutate(
        x,
        mod = sqrt(re^2 + im^2),
        arg = atan2(im, re),
        .keep = .keep
      )
    } else {
      if (has_polr(x)) {
        x
      } else {
        stop("No valid fft representation.")
      }
    }
  }
}

#' @export
set_repr <- function(x, repr) {
  res <- .get_dim_cols(x) |>
    dplyr::mutate(fx = get_fx(x))
  res <- .set_repr(res, repr[1])
  if (length(repr) > 1) {
    for (i in 2:length(repr)) {
      res <- .set_repr(res, repr[i], .keep = "all")
    }
  }
  res
}
