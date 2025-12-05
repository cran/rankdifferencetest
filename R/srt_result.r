#' @title
#' Prepare the return object for `srt()` and `rdt()`.
#'
#' @description
#' Prepares the return object for [srt()] and [rdt()].
#'
#' @param x
#' (named list)\cr
#' The object returned by [srt_pmedian()] in an srt pipeline.
#'
#' @returns
#' Named `list`
#'
#' @keywords internal
#' @export
srt_result <- function(x) {
  info <- c(
    "p_value_method",
    "pseudomedian_method",
    "conf_method",
    "conf_level_achieved",
    "n_sample",
    "n_analytic",
    "n_zeros",
    "n_signed",
    "n_ties",
    "data_type",
    "focal_name",
    "reference_name"
  )
  res <- list(
    p_value = x$p_value,
    statistic = x$statistic,
    pseudomedian = x$pseudomedian,
    lower = x$lower,
    upper = x$upper,
    method = x$method,
    info = x[info],
    call = x$call
  )
  class(res) <- c(x$class, "list")

  res
}
