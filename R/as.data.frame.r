#' Coerce to a data frame
#'
#' Coerce objects to a `data.frame` in the style of `broom::tidy()`.
#'
#' @param x
#' Objects returned by:
#' - [rdt()] or [rdt2()]
#' - [srt()] or [srt2()]
#' - [rdpmedian()] or [rdpmedian2()]
#' - [pmedian()] or [pmedian2()]
#'
#' @param ...
#' Unused arguments.
#'
#' @returns
#' `data.frame`
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # as.data.frame() examples
#' #----------------------------------------------------------------------------
#' library(rankdifferencetest)
#'
#' # Use example data from Kornbrot (1990)
#' data <- kornbrot_table1
#'
#' rdt(
#'   data = data,
#'   formula = placebo ~ drug,
#'   conf_level = 0.95,
#'   alternative = "two.sided",
#'   distribution = "asymptotic",
#'   zero_method = "wilcoxon",
#'   correct = FALSE
#' ) |>
#'   as.data.frame()
#'
#' rdpmedian(
#'   data = data,
#'   formula = placebo ~ drug,
#'   conf_level = 0.95
#' ) |>
#'   as.data.frame()
#'
#' @rdname as.data.frame
#' @export
as.data.frame.rdt <- function(x, ...) {
  .as.data.frame.srt(x)
}

#' @rdname as.data.frame
#' @export
as.data.frame.srt <- function(x, ...) {
  .as.data.frame.srt(x)
}

#' @rdname as.data.frame
#' @export
as.data.frame.rdpmedian <- function(x, ...) {
  .as.data.frame.pmedian(x)
}

#' @rdname as.data.frame
#' @export
as.data.frame.pmedian <- function(x, ...) {
  .as.data.frame.pmedian(x)
}

#-------------------------------------------------------------------------------
# Generic helpers
#-------------------------------------------------------------------------------
.as.data.frame.srt <- function(x) {
  res <- if (is.null(x$lower)) {
    list(
      estimate = x$pseudomedian,
      statistic = x$statistic,
      p.value = x$p_value,
      method = x$method,
      alternative = x$call$alternative
    )
  } else {
    list(
      estimate = x$pseudomedian,
      statistic = x$statistic,
      p.value = x$p_value,
      conf.low = x$lower,
      conf.high = x$upper,
      conf.level = x$call$conf_level,
      method = x$method,
      alternative = x$call$alternative
    )
  }

  attr(res, "names") <- names(res)
  attr(res, "row.names") <- .set_row_names(1L)
  attr(res, "class") <- "data.frame"

  # Return
  res
}

.as.data.frame.pmedian <- function(x) {
  res <- if (is.null(x$lower)) {
    list(
      estimate = x$pseudomedian,
      method = x$method
    )
  } else {
    list(
      estimate = x$pseudomedian,
      conf.low = x$lower,
      conf.high = x$upper,
      conf.level = x$call$conf_level,
      method = x$method
    )
  }

  attr(res, "names") <- names(res)
  attr(res, "row.names") <- .set_row_names(1L)
  attr(res, "class") <- "data.frame"

  # Return
  res
}
