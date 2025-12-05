#' @title
#' Calculate the test statistic for a signed-rank test.
#'
#' @description
#' Calculates the test statistic for a signed-rank test.
#'
#' @param x
#' (named list)\cr
#' The object returned by [rankdifferencetest::srt_method()] in a 'srt pipeline'.
#'
#' @param warn_zd
#' (Scalar logical: `c(FALSE, TRUE)`)\cr
#' Used for root-finding algorithm 'Z(d)' in confidence interval test inversion.
#' If `TRUE` return a warning and exit with value `statistic = 0`.
#'
#' @param ...
#' dots arguments are not used.
#'
#' @returns
#' list
#'
#' It carries forward the class selected by [rankdifferencetest::srt_method()].
#'
#' @importFrom stats setNames
#' @keywords internal
#' @rdname srt_statistic
#' @export
srt_statistic <- function(x, ...) {
  if (x$invalid_data) {
    return(x)
  }
  UseMethod("srt_statistic")
}

#' @keywords internal
#' @rdname srt_statistic
#' @export
srt_statistic.wilcoxon <- function(x, ...) {
  x$statistic <- setNames(x$wplus, "W+")
  x
}

#' @keywords internal
#' @rdname srt_statistic
#' @export
srt_statistic.shift <- function(x, ...) {
  x$statistic <- setNames(x$wplus, "W+")
  x
}

#' @keywords internal
#' @rdname srt_statistic
#' @export
srt_statistic.asymptotic <- function(x, warn_zd = FALSE, ...) {
  # Expected values under null
  e_wplus <- 0.5 * sum(x$ranks)
  var_wplus <- 0.25 * sum(x$ranks^2)

  if (warn_zd) {
    if (!is.finite(var_wplus) || var_wplus <= 0) {
      msg <- "When evaluating Z(d), asymptotic variance was zero or nonfinite.\nReturning Z(d)=0"
      warning(msg)
      z <- 0
      z <- setNames(z, "Z")
      x$statistic <- z
      return(x)
    }
  }

  sd_wplus <- sqrt(var_wplus)

  # Continuity correction
  cc <- if (x$call$correct) {
    switch(
      x$call$alternative,
      "two.sided" = 0.5 * sign(x$wplus - e_wplus),
      "greater" = 0.5,
      "less" = -0.5
    )
  } else {
    0
  }

  # Test statistic
  z <- (x$wplus - e_wplus - cc) / sd_wplus
  z <- setNames(z, "Z")

  # Return
  x$statistic <- z
  x
}
