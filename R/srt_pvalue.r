#' @title
#' Calculate the p-value for a signed-rank test.
#'
#' @description
#' Calculates the p-value for a signed-rank test.
#'
#' @param x
#' (named list)\cr
#' The object returned by [rankdifferencetest::srt_statistic()] in a 'srt pipeline'.
#'
#' @param ...
#' dots arguments are not used.
#'
#' @returns
#' list
#'
#' @keywords internal
#' @export
srt_pvalue <- function(x, ...) {
  if (x$invalid_data) {
    return(x)
  }
  UseMethod("srt_pvalue")
}

#' @importFrom stats psignrank
#' @keywords internal
#' @export
srt_pvalue.wilcoxon <- function(x, ...) {
  n <- x$n_signed
  statistic <- x$statistic

  p <- switch(
    x$call$alternative,
    two.sided = {
      p <- if (statistic > (n * (n + 1) / 4)) {
        psignrank(statistic - 1, n, lower.tail = FALSE)
      } else {
        psignrank(statistic, n)
      }
      min(2 * p, 1)
    },
    greater = psignrank(statistic - 1, n, lower.tail = FALSE),
    less = psignrank(statistic, n)
  )

  # Return
  x$p_value <- p
  x$p_value_method <- "Wilcoxon"
  x
}

#' @importFrom exactRankTests pperm
#' @keywords internal
#' @export
srt_pvalue.shift <- function(x, ...) {
  p <- pperm(
    q = x$statistic,
    scores = x$ranks,
    m = length(x$ranks),
    alternative = x$call$alternative,
    pprob = FALSE
  )

  # Return
  x$p_value <- p
  x$p_value_method <- "Shift-algorithm"
  x
}

#' @importFrom stats pnorm
#' @keywords internal
#' @export
srt_pvalue.asymptotic <- function(x, ...) {
  statistic <- x$statistic

  p <- switch(
    x$call$alternative,
    two.sided = 2 * min(pnorm(statistic), pnorm(statistic, lower.tail = FALSE)),
    greater = pnorm(statistic, lower.tail = FALSE),
    less = pnorm(statistic)
  )

  # Return
  x$p_value <- p
  x$p_value_method <- "Asymptotic"
  x
}
