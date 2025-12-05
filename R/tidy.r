#' Coerce to a data frame
#'
#' Coerce objects to a 'tibble' data frame in the style of `broom::tidy()`.
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
#' @return `tbl_df`
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # tidy() examples
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
#'   tidy()
#'
#' rdpmedian(
#'   data = data,
#'   formula = placebo ~ drug,
#'   conf_level = 0.95
#' ) |>
#'   tidy()
#'
#' @rdname tidy
#' @export
tidy.rdt <- function(x, ...) {
  .tidy.srt(x)
}

#' @rdname tidy
#' @export
tidy.srt <- function(x, ...) {
  .tidy.srt(x)
}

#' @rdname tidy
#' @export
tidy.rdpmedian <- function(x, ...) {
  .tidy.pmedian(x)
}

#' @rdname tidy
#' @export
tidy.pmedian <- function(x, ...) {
  .tidy.pmedian(x)
}

#' @importFrom generics tidy
#' @export
generics::tidy

#-------------------------------------------------------------------------------
# Generic helpers
#-------------------------------------------------------------------------------
.tidy.srt <- function(x) {
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
  attr(res, "class") <- c("tbl_df", "tbl", "data.frame")

  # Return
  res
}

.tidy.pmedian <- function(x) {
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
  attr(res, "class") <- c("tbl_df", "tbl", "data.frame")

  # Return
  res
}
