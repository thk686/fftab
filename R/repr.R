#' Convert a `fftab` Object Between Representations
#'
#' These functions convert a `fftab` object to a specified representation:
#' - **`to_cplx()`**: Converts to complex representation (`fx`).
#' - **`to_rect()`**: Converts to rectangular representation (`re`, `im`).
#' - **`to_polr()`**: Converts to polar representation (`mod`, `arg`).
#'
#' @param x A `fftab` object.
#' @param .keep Specifies which columns to retain. See [dplyr::mutate()].
#'
#' @return A modified `fftab` object containing the specified representation:
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
#' fftab(c(1, 0, -1, 0)) |> to_cplx()
#'
#' fftab(c(1, 0, -1, 0)) |> to_rect()
#'
#' fftab(c(1, 0, -1, 0)) |> to_polr()
#'
#' @seealso [has_cplx()], [get_repr()]
#' @rdname to_cplx
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

#' Check Representations of a `fftab` Object
#'
#' These functions check if specific representations are present in a `fftab` object:
#'
#' - **`has_cplx()`**: Checks if the object has complex representation (`fx` column).
#' - **`has_rect()`**: Checks if the object has rectangular representation (`re`, `im` columns).
#' - **`has_polr()`**: Checks if the object has polar representation (`mod`, `arg` columns).
#'
#' @param x A `fftab` object.
#'
#' @return Logical value (`TRUE` or `FALSE`) indicating whether the specified representation exists.
#'
#' @examples
#' fftab(c(1, 0, -1, 0)) |> has_cplx()
#'
#' fftab(c(1, 0, -1, 0)) |> has_rect()
#'
#' fftab(c(1, 0, -1, 0)) |> has_polr()
#'
#' @seealso [to_cplx()], [get_repr()]
#' @rdname has_cplx
#' @export
has_cplx <- function(x) {
  !is.null(x[["fx"]])
}

#' @rdname has_cplx
#' @export
has_rect <- function(x) {
  !is.null(x[["re"]]) && !is.null(x[["im"]])
}

#' @rdname has_cplx
#' @export
has_polr <- function(x) {
  !is.null(x[["mod"]]) && !is.null(x[["arg"]])
}

#' Manage Representations of a `fftab` Object
#'
#' These functions handle representation management for a `fftab` object:
#'
#' - **`get_repr()`**: Retrieve current representations.
#' - **`can_repr()`**: Check if the object supports specific representations.
#' - **`set_repr()`**: Convert the object to one or more specified representations.
#'
#' @param x A `fftab` object.
#' @param repr For `can_repr()`, a character vector specifying representations (`"polr"`, `"rect"`, `"cplx"`).
#'
#' @return
#' - **`can_repr()`**: Logical value (`TRUE` or `FALSE`) indicating if the object supports the specified representations.
#' - **`get_repr()`**: A character vector of current representations.
#' - **`set_repr()`**: A modified `fftab` object with the specified representation(s).
#'
#' @examples
#' fftab(c(1, 0, -1, 0)) |> can_repr("cplx")
#'
#' fftab(c(1, 0, -1, 0)) |> get_repr()
#'
#' fftab(c(1, 0, -1, 0)) |> set_repr(c("polr", "rect"))
#'
#' @seealso [to_cplx()], [has_cplx()]
#' @rdname can_repr
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
set_repr <- function(x, repr) {
  res <- .get_dim_cols(x) |> dplyr::mutate(fx = get_fx(x))
  res <- .set_repr(res, repr[1])
  if (length(repr) > 1) {
    for (i in 2:length(repr)) {
      res <- .set_repr(res, repr[i], .keep = "all")
    }
  }
  res
}

#' Add Additional Representations to Fourier Coefficients
#'
#' These functions add additional representations to a `fftab` object without removing or modifying existing representations.
#'
#' @param x A `fftab` object containing Fourier coefficients and associated metadata.
#'
#' @return A `fftab` object with the additional representation included.
#'
#' @details
#' - **`add_cplx()`**: Adds a **complex** (`"cplx"`) representation to the Fourier coefficients.
#' - **`add_rect()`**: Adds a **rectangular** (`"rect"`) representation to the Fourier coefficients.
#' - **`add_polr()`**: Adds a **polar** (`"polr"`) representation to the Fourier coefficients.
#'
#' These functions are useful for working with multiple representations simultaneously without overwriting existing data.
#'
#' @seealso
#' - [fftab()]
#'
#' @examples
#' matrix(1:9, 3) |>
#'   fftab() |>
#'   print(n = 3) |>
#'   add_polr() |>
#'   print(n = 3) |>
#'   add_rect() |>
#'   print(n = 3) |>
#'   add_cplx() |>
#'   print(n = 3)
#'
#' @export
add_cplx <- function(x) {
  .set_repr(x, "cplx", .keep = "all")
}

#' @rdname add_cplx
#' @export
add_rect <- function(x) {
  .set_repr(x, "rect", .keep = "all")
}

#' @rdname add_cplx
#' @export
add_polr <- function(x) {
  .set_repr(x, "polr", .keep = "all")
}

#' Convert between cyclic and angular frequencies
#'
#' @description
#' These functions convert dimensions of frequency from cyclic (measured in cycles)
#' to angular (measured in radians) or vice versa. This transformation scales
#' dimensions by a factor of 2 * pi.
#'
#' - `to_angf()`: Converts from cyclic to angular frequency.
#' - `to_cycf()`: Converts from angular to cyclic frequency.
#'
#' @param x An `fftab` object containing frequency dimensions. Must include
#'   columns prefixed with `.dim_` and have an attribute `.is_angular` indicating
#'   the frequency type.
#'
#' @returns
#' An `fftab` object with dimensions scaled appropriately and the `.is_angular`
#' attribute updated.
#'
#' @examples
#' # Convert cyclic to angular frequencies
#' rnorm(64) |>
#'   fftab() |>
#'   to_angf() |>
#'   to_rect() |>
#'   dplyr::slice_max(abs(.dim_1), n = 5)
#'
#' # Convert angular back to cyclic frequencies
#' rnorm(64) |>
#'   fftab() |>
#'   to_angf() |>
#'   to_cycf() |>
#'   to_rect() |>
#'   dplyr::slice_max(abs(.dim_1), n = 5)
#'
#' @export
to_angf <- function(x) {
  if (.is_angular(x)) {
    x
  } else {
    x |>
      dplyr::mutate(
        dplyr::across(
          dplyr::starts_with(".dim_"), ~ . * 2 * pi)) |>
      structure(.is_angular = TRUE)
  }
}

#' @rdname to_angf
#' @export
to_cycf <- function(x) {
  if (.is_angular(x)) {
    x |>
      dplyr::mutate(
        dplyr::across(
          dplyr::starts_with(".dim_"), ~ . / 2 / pi)) |>
      structure(.is_angular = FALSE)
  } else {
    x
  }
}

