#' Convert a `tidy_fft` Object Between Representations
#'
#' These functions convert a `tidy_fft` object to a specified representation:
#' - **`to_cplx()`**: Converts to complex representation (`fx`).
#' - **`to_rect()`**: Converts to rectangular representation (`re`, `im`).
#' - **`to_polr()`**: Converts to polar representation (`mod`, `arg`).
#'
#' @param x A `tidy_fft` object.
#' @param .keep Specifies which columns to retain. See [dplyr::mutate()].
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
#' tidy_fft(c(1, 0, -1, 0)) |> to_cplx()
#'
#' tidy_fft(c(1, 0, -1, 0)) |> to_rect()
#'
#' tidy_fft(c(1, 0, -1, 0)) |> to_polr()
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

#' Check Representations of a `tidy_fft` Object
#'
#' These functions check if specific representations are present in a `tidy_fft` object:
#'
#' - **`has_cplx()`**: Checks if the object has complex representation (`fx` column).
#' - **`has_rect()`**: Checks if the object has rectangular representation (`re`, `im` columns).
#' - **`has_polr()`**: Checks if the object has polar representation (`mod`, `arg` columns).
#'
#' @param x A `tidy_fft` object.
#'
#' @return Logical value (`TRUE` or `FALSE`) indicating whether the specified representation exists.
#'
#' @examples
#' tidy_fft(c(1, 0, -1, 0)) |> has_cplx()
#'
#' tidy_fft(c(1, 0, -1, 0)) |> has_rect()
#'
#' tidy_fft(c(1, 0, -1, 0)) |> has_polr()
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

#' Manage Representations of a `tidy_fft` Object
#'
#' These functions handle representation management for a `tidy_fft` object:
#'
#' - **`get_repr()`**: Retrieve current representations.
#' - **`can_repr()`**: Check if the object supports specific representations.
#' - **`set_repr()`**: Convert the object to one or more specified representations.
#'
#' @param x A `tidy_fft` object.
#' @param repr For `can_repr()`, a character vector specifying representations (`"polr"`, `"rect"`, `"cplx"`).
#'
#' @return
#' - **`can_repr()`**: Logical value (`TRUE` or `FALSE`) indicating if the object supports the specified representations.
#' - **`get_repr()`**: A character vector of current representations.
#' - **`set_repr()`**: A modified `tidy_fft` object with the specified representation(s).
#'
#' @examples
#' tidy_fft(c(1, 0, -1, 0)) |> can_repr("cplx")
#'
#' tidy_fft(c(1, 0, -1, 0)) |> get_repr()
#'
#' tidy_fft(c(1, 0, -1, 0)) |> set_repr(c("polr", "rect"))
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
