#' @export
gauss_filter <- function(x, sd) {
  kern <- exp(-pi^2 * sd^2 * get_l2sq(x))
  to_polr(x) |>
    dplyr::mutate(mod = kern * mod)
}
