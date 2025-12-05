#' Print results
#'
#' Prints results in the style of `htest`.
#'
#' @param x
#' Objects returned by:
#' - [rdt()] or [rdt2()]
#' - [srt()] or [srt2()]
#' - [rdpmedian()] or [rdpmedian2()]
#' - [pmedian()] or [pmedian2()]
#'
#' @param digits
#' (Scalar integer: `3`; `[1, Inf)`)\cr
#' The number of significant digits to use.
#'
#' @param ...
#' Unused additional arguments.
#'
#' @returns
#' `invisible(x)`
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # print() examples
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
#'   print()
#'
#' srt(
#'   data = data,
#'   formula = placebo ~ drug,
#'   conf_level = 0.95,
#'   alternative = "two.sided",
#'   distribution = "asymptotic",
#'   zero_method = "wilcoxon",
#'   correct = FALSE
#' ) |>
#'   print()
#'
#' rdpmedian(
#'   data = data,
#'   formula = placebo ~ drug,
#'   conf_level = 0.95
#' ) |>
#'   print()
#'
#' pmedian(
#'   data = data,
#'   formula = placebo ~ drug,
#'   conf_level = 0.95
#' ) |>
#'   print()
#'
#' @rdname print
#' @export
print.rdt <- function(x, digits = 3, ...) {
  .print.generic(x, digits)
  invisible(x)
}

#' @rdname print
#' @export
print.srt <- function(x, digits = 3, ...) {
  .print.generic(x, digits)
  invisible(x)
}

#' @rdname print
#' @export
print.rdpmedian <- function(x, digits = 3, ...) {
  .print.generic(x, digits)
  invisible(x)
}

#' @rdname print
#' @export
print.pmedian <- function(x, digits = 3, ...) {
  .print.generic(x, digits)
  invisible(x)
}

#-------------------------------------------------------------------------------
# Generic helper
#-------------------------------------------------------------------------------
.print.generic <- function(x, digits) {
  if (
    !(is.finite(digits) &&
      is.numeric(digits) &&
      length(digits) == 1L &&
      is.finite(digits) &&
      digits >= 1)
  ) {
    stop("Argument 'digits' must be a positive scalar integer.")
  }

  method <- x$method

  focal_name <- x$info$focal_name
  reference_name <- x$info$reference_name
  data_names <- c(focal_name, reference_name)
  data_type <- x$info$data_type
  data_diffs <- paste0(data_names, collapse = " - ")
  data <- paste0(data_type, ": ", data_diffs)

  tstat <- x$statistic
  if (!is.null(tstat)) {
    tstat <- paste0(names(tstat), " = ", format(tstat, digits = digits))
  }

  p <- x$p_value
  if (!is.null(p)) {
    p <- format.pval(x$p_value, digits = digits)
    p <- paste("p", if (startsWith(p, "<")) p else paste("=", p))
  }

  alternative <- x$call$alternative
  if (!is.null(alternative)) {
    alternative <- switch(
      alternative,
      two.sided = "not equal to",
      less = "less than",
      greater = "greater than"
    )
    alternative <- paste0(
      "Alternative hypothesis: True location shift of ",
      tolower(data_type),
      " is ",
      alternative,
      " ",
      x$call$mu
    )
  }

  estimate <- paste0("Pseudomedian: ", format(x$pseudomedian, digits = digits))

  ci <- if (!is.null(x$lower)) {
    conf_level <- if (
      !is.null(x$info$conf_level_achieved) &&
        !is.na(x$info$conf_level_achieved) &&
        x$info$conf_level_achieved < x$call$conf_level
    ) {
      x$info$conf_level_achieved
    } else {
      x$call$conf_level
    }
    paste0(
      format(100 * conf_level, digits = digits),
      "% CI: ",
      paste0(
        format(x$lower, digits = digits),
        ", ",
        format(x$upper, digits = digits)
      )
    )
  }

  cat("\n")
  cat(strwrap(method, width = 62), sep = "\n")
  if (!is.null(p)) {
    cat("\n")
    cat(p)
  }
  if (!is.null(tstat)) {
    cat("\n")
    cat(tstat)
  }
  cat("\n")
  cat(strwrap(data, width = 62), sep = "\n")
  if (!is.null(alternative)) {
    cat(strwrap(alternative, width = 62), sep = "\n")
  }
  cat("\n")
  cat(estimate)
  if (!is.null(ci)) {
    cat("\n")
    cat(ci)
  }
}
