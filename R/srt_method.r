#' @title
#' Methods selector for the signed-rank test.
#'
#' @description
#' Performs four tasks:
#' 1. Selects the appropriate methods for a signed-rank test.
#' 2. Selects the appropriate methods for a confidence interval.
#' 3. Checks for invalid data before continuing with calculating test results.
#' 4. Calculates the number of tied values among the ranked absolute differences.
#'
#' @details
#' The signed-rank method is chosen by:
#'
#' \itemize{
#'   \item If ties are absent and zeros are absent, the exact null distribution of \eqn{W^+} depends only on \eqn{n_{\text{signed}} = |\mathcal{I}_{\neq 0}|} and can be obtained from [stats::psignrank()] or [stats::qsignrank()].
#'   \item If ties are present or when zeros are present, the exact null distribution is the distribution of a random subset-sum of the `ranks` weights.
#'         Methods such as [exactRankTests::pperm()] or [exactRankTests::qperm()] are appropriate.
#'   \item For normal approximations, the mean and variance of \eqn{W^+} under the symmetry null are \eqn{\mathrm{E}(W^+) = \tfrac{1}{2}\sum r_i} and \eqn{\mathrm{Var}(W^+) = \tfrac{1}{4}\sum r_i^2} computed over `ranks`.
#' }
#'
#' The `srt` pipeline uses two intermediate classes for internal routing of calculations.
#' 1. The test class (`"wilcoxon"`, `"shift"`, or `"asymptotic"`) is determined by properties of the data and `x$call$distribution`.
#' 2. The confidence interval class (`"inversion"`, `"bootstrap"`, `"none"`) is determined by `x$call$conf_level` and `x$call$conf_method`.
#' The final user facing result class is either `srt` or `rdt`, determined by the calling function.
#'
#' @param x
#' (named list)\cr
#' The object returned by [rankdifferencetest::srt_ranks()].
#'
#' @param test_class
#' (Scalar character)\cr
#' This argument should only be used inside `.zd()` (asymptotic test inversion for CI), for speed, and should be left to the default value in an 'srt pipeline'.
#' Must be one of `"auto"`, `"wilcoxon"`, `"shift"`, or `"asymptotic"`.
#' If `"auto"` (default), the test class class will automatically be chosen based on properties of the data and `x$call$distribution`.
#' Otherwise, the function will return early with class set to `test_class`.
#'
#' @returns
#' list
#'
#' @keywords internal
srt_method <- function(x, test_class = "auto") {
  #-----------------------------------------------------------------------------
  # Short circuit
  # Used for .zd() asymptotic test inversion for CI. Not used in full srt
  # pipeline, so setting `x$n_ties` and `x$class` not required.
  #-----------------------------------------------------------------------------
  if (test_class != "auto") {
    x$call$distribution <- switch(
      test_class,
      "wilcoxon" = "exact",
      "shift" = "exact",
      "asymptotic" = "asymptotic"
    )
    x$invalid_data <- FALSE
    class(x) <- c(test_class, "list")
    return(x)
  }

  #-----------------------------------------------------------------------------
  # Set distribution if needed
  #-----------------------------------------------------------------------------
  if (x$call$distribution == "auto") {
    if (x$n_signed < 50L) {
      x$call$distribution <- "exact"
    } else {
      x$call$distribution <- "asymptotic"
    }
  }

  #-----------------------------------------------------------------------------
  # Ties (among ranks with nonzero differences)
  #-----------------------------------------------------------------------------
  n_ties <- if (anyDuplicated(x$ranks) > 0L) {
    sum(ftabulate(x$ranks) - 1L)
  } else {
    0L
  }

  #-----------------------------------------------------------------------------
  # Test class
  #-----------------------------------------------------------------------------
  plan <- paste(x$call$distribution, x$call$zero_method, sep = "_")

  test_class <- switch(
    plan,
    "exact_wilcoxon" = if (n_ties == 0L && x$n_zeros == 0L) {
      "wilcoxon"
    } else {
      "shift"
    },
    "exact_pratt" = if (n_ties == 0L && x$n_zeros == 0L) {
      "wilcoxon"
    } else {
      "shift"
    },
    "asymptotic_wilcoxon" = "asymptotic",
    "asymptotic_pratt" = "asymptotic"
  )

  #-----------------------------------------------------------------------------
  # Confidence interval class
  #-----------------------------------------------------------------------------
  conf_class <- if (x$call$conf_level == 0) {
    "none"
  } else if (x$call$conf_method == "inversion") {
    test_class
  } else {
    "bootstrap"
  }

  #-----------------------------------------------------------------------------
  # Result class
  #-----------------------------------------------------------------------------
  result_class <- switch(
    as.character(x$call[[1L]]),
    "rdt" = "rdt",
    "rdt2" = "rdt",
    "srt" = "srt",
    "srt2" = "srt"
  )

  #-----------------------------------------------------------------------------
  # Method description
  #-----------------------------------------------------------------------------
  pratt <- if (x$call$zero_method == "pratt") "-Pratt"

  switch(
    x$method,
    "Kornbrot-Wilcoxon" = {
      test_type <- "rank difference test"
      data_type <- "Paired differences of ranks"
    },
    "Wilcoxon" = {
      test_type <- "signed-rank test"
      data_type <- if (is.null(x$reference_name)) {
        "Observations"
      } else {
        "Paired differences"
      }
    },
    test_type <- data_type <- NULL
  )

  method <- paste0(
    x$call$distribution,
    " ",
    x$method,
    pratt,
    " ",
    test_type
  )

  #-----------------------------------------------------------------------------
  # Invalid data checks
  #-----------------------------------------------------------------------------
  # skip test if data is bad
  invalid_data <- x$n_analytic < 3L ||
    x$n_signed == 0L ||
    all(x$diffs == x$diffs[1L])

  if (invalid_data) {
    msg <- if (x$n_analytic < 3L) {
      valid_data <- switch(
        x$method,
        "Kornbrot-Wilcoxon" = "(non-missing)",
        "Wilcoxon" = "(non-missing and finite)"
      )
      paste0(
        "The number of valid observations ",
        valid_data,
        " must be at least 3.\nReturned NA for results."
      )
    } else if (x$n_signed == 0L) {
      paste0(
        "No ",
        tolower(data_type),
        " left after removing zeros.\nReturned NA for results."
      )
    } else if (all(x$diffs == x$diffs[1L])) {
      paste0(
        data_type,
        " are constant/zero variance.\nReturned NA for results."
      )
    }
    warning(msg)

    statistic <- NA_real_
    statistic_name <- switch(
      test_class,
      "wilcoxon" = "W+",
      "shift" = "W+",
      "asymptotic" = "Z"
    )

    pseudomedian <- if (x$n_analytic > 0L) {
      Hmisc::pMedian(x$diffs)
    } else {
      NA_real_
    }
    conf <- if (x$call$conf_level > 0) NA_real_ else NULL

    out <- list(
      p_value = NA_real_,
      statistic = setNames(statistic, statistic_name),
      pseudomedian = pseudomedian,
      lower = conf,
      upper = conf,
      conf_method = conf,
      conf_level_achieved = conf,
      p_value_method = NA_real_,
      pseudomedian_method = NA_real_
    )
    x <- c(x, out)
  }

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  x$method <- method
  out <- list(
    n_ties = n_ties,
    class = result_class,
    test_type = test_type,
    data_type = data_type,
    invalid_data = invalid_data
  )
  x <- c(x, out)
  class(x) <- c(conf_class, test_class, "list")

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  x
}
