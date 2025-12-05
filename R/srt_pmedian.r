# Code and methods below are based on and modified from:
# - stats::wilcox.test()
# - exactRankTests::wilcox.exact()
# The stats package is Copyright R Core Team, GPL-2 | GPL-3
# The exactRankTests package is Copyright Torsten Hothorn, GPL-2 | GPL-3

#' @title
#' Calculate the pseudomedian and confidence interval for a signed-rank test.
#'
#' @description
#' Calculates the pseudomedian and confidence interval for a signed-rank test.
#'
#' @references
#'
#' \insertRef{exactRankTests}{rankdifferencetest}
#'
#' \insertRef{R}{rankdifferencetest}
#'
#' @param x
#' (named list)\cr
#' The object returned by [rankdifferencetest::srt_pvalue()] in a 'srt pipeline'.
#'
#' @param ... dots arguments are not used.
#'
#' @returns
#' list
#'
#' @keywords internal
#' @export
srt_pmedian <- function(x, ...) {
  if (x$invalid_data) {
    return(x)
  }
  UseMethod("srt_pmedian")
}

#' @keywords internal
#' @export
srt_pmedian.none <- function(x, ...) {
  #-----------------------------------------------------------------------------
  # Only pseudomedian; no CI
  # Use original data with possible zeros that would have otherwise been removed
  # by 'wilcoxon' zero handling.
  #-----------------------------------------------------------------------------
  out <- list(
    pseudomedian = Hmisc::pMedian(x$diffs_original),
    pseudomedian_method = "Hodges-Lehmann",
    conf_method = NULL,
    conf_level_achieved = NULL
  )
  return(c(x, out))
}

#' @importFrom stats psignrank qsignrank
#' @keywords internal
#' @export
srt_pmedian.wilcoxon <- function(x, ...) {
  #-----------------------------------------------------------------------------
  # The srt_pmedian.wilcoxon() function is only used when possible mu-shifted
  # data results in zero-free data and tie-free ranks.
  # `x$diffs + x$call$mu` provides mu-shift-free differences. However, this may
  # introduce zeros (that 'wilcoxon' would have otherwise removed) and/or tied
  # ranks. I'll follow wilcox.test approach and ignore this.
  #-----------------------------------------------------------------------------
  diffs <- if (x$call$mu != 0) x$diffs + x$call$mu else x$diffs

  #-----------------------------------------------------------------------------
  # Confidence interval
  #-----------------------------------------------------------------------------
  n <- length(diffs)
  alpha <- 1 - x$call$conf_level
  walsh <- walsh(x = diffs, sort = TRUE)
  pseudomedian <- fmedian(walsh, is_sorted = TRUE)

  switch(
    x$call$alternative,
    two.sided = {
      # Lower critical value at tail alpha/2
      cv_lower <- max(1, qsignrank(alpha / 2, n))
      # Upper critical value at tail 1-alpha/2 (equal to qsignrank(1 - alpha / 2, n))
      cv_upper <- n * (n + 1) / 2 - cv_lower
      # psignrank(q, n, lower.tail = TRUE) returns P(W+ <= q).
      # Substituting q = (cv_lower - 1) makes it P(W+ < cv_lower).
      # In this exact case, both tails are symmetric, so multiply by 2 to get alpha.
      alpha_achieved <- 2 * psignrank(cv_lower - 1, n)
      lower <- walsh[cv_lower]
      upper <- walsh[cv_upper + 1]
    },
    greater = {
      # Lower critical value at tail alpha
      cv_lower <- max(1, qsignrank(alpha, n))
      # P(W+ < cv_lower)
      alpha_achieved <- psignrank(cv_lower - 1, n)
      lower <- walsh[cv_lower]
      upper <- Inf
    },
    less = {
      # Upper critical value at tail 1-alpha
      cv_upper <- max(1, qsignrank(alpha, n, lower.tail = FALSE))
      # P(W+ > cv_lower)
      alpha_achieved <- psignrank(cv_upper, n, lower.tail = FALSE)
      lower <- -Inf
      upper <- walsh[cv_upper + 1]
    }
  )

  if (alpha_achieved - alpha > alpha / 2) {
    warning("Requested 'conf_level' not achievable.")
  }

  conf_level_achieved <- 1 - signif(alpha_achieved, 2)
  pseudomedian_method <- "Hodges-Lehmann"
  conf_method <- "Inversion of exact Wilcoxon signed-rank test"

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  out <- list(
    pseudomedian = pseudomedian,
    lower = lower,
    upper = upper,
    conf_level_achieved = conf_level_achieved,
    pseudomedian_method = pseudomedian_method,
    conf_method = conf_method
  )

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  c(x, out)
}

#' @importFrom exactRankTests qperm pperm
#' @keywords internal
#' @export
srt_pmedian.shift <- function(x, ...) {
  #-----------------------------------------------------------------------------
  # `x$diffs + x$call$mu` provides mu-shift-free differences. However, this may
  # introduce zeros (that 'wilcoxon' would have otherwise removed) and/or tied
  # ranks. I'll follow wilcox.exact approach and ignore this.
  #-----------------------------------------------------------------------------
  diffs <- if (x$call$mu != 0) x$diffs + x$call$mu else x$diffs
  ranks <- x$ranks
  n <- x$n_signed
  digits_rank <- x$call$digits_rank
  zero_method <- x$call$zero_method
  alpha <- 1 - x$call$conf_level

  #-----------------------------------------------------------------------------
  # .wd(): W+(d), sum of positive ranks at shift d
  # wd_vec: Vector of W+(d) across all Walsh candidates
  #
  # Note that wilcox.exact() always uses the 'pratt' method of zero handling here.
  # This may be because it preserves the rank structure across candidate shifts,
  # maintains consistency with the test statistic definition, and removing zeros
  # dynamically during CI construction could lead to inconsistencies between the
  # test statistic and the CI inversion logic...
  #
  # Small candidate shifts d imply many positive differences leading to large .wd(d)
  # Large candidate shifts d imply few positive differences leading to small .wd(d)
  #-----------------------------------------------------------------------------
  .wd <- function(d) {
    srt_ranks(
      x = list(diffs = diffs - d),
      call = list(
        mu = 0,
        zero_method = "pratt",
        digits_rank = digits_rank
      ),
      warn_wd = TRUE
    ) |>
      getElement("wplus")
  }

  walsh <- walsh(x = diffs, sort = TRUE)

  wd_vec <- if (x$n_ties == 0L && zero_method == "wilcoxon") {
    # Short circuit
    ranks_wd <- rank(abs(diffs))
    if (length(walsh) != sum(ranks_wd)) {
      stop("length(walsh) != sum(ranks_wd)") # Validate
    }
    sum(ranks_wd):1L
  } else {
    # Are duplicate walsh redundant for inversion? probably not, so keep them...
    vapply(walsh, .wd, numeric(1L))
  }

  #-----------------------------------------------------------------------------
  # Invert the exact distribution with qperm to get CI
  #
  # The step function .wd(d) decreases with d.
  # To find the lower bound, you look for the largest d where the statistic is
  # still above the upper critical value.
  # To find the upper bound, you look for the smallest d where the statistic is
  # still below the lower critical value.
  #-----------------------------------------------------------------------------
  # Helpers to get P(W <= w) and P(W >= w) for the achieved confidence level
  .alpha_lower <- function(w) {
    exactRankTests::pperm(
      q = w,
      scores = ranks,
      m = length(ranks),
      paired = TRUE,
      alternative = "less"
    )
  }
  .alpha_upper <- function(w) {
    exactRankTests::pperm(
      q = w,
      scores = ranks,
      m = length(ranks),
      paired = TRUE,
      alternative = "greater"
    )
  }

  if (x$call$alternative == "two.sided") {
    # Use these to define acceptance region.
    # lower critical value for signed rank sum at tail probability alpha/2
    cv_lower <- exactRankTests::qperm(
      p = alpha / 2,
      scores = ranks,
      m = length(ranks),
      paired = TRUE
    )
    # upper critical value for signed rank sum at tail probability 1-alpha/2
    cv_upper <- exactRankTests::qperm(
      p = 1 - alpha / 2,
      scores = ranks,
      m = length(ranks),
      paired = TRUE
    )

    # CI inversion finds shift value d where w(d) crosses the acceptance region.
    # lower bound
    if (cv_upper >= max(wd_vec, na.rm = TRUE)) {
      # edge case, set lower as smallest shift available.
      lower <- min(walsh, na.rm = TRUE)
    } else {
      # largest candidate d where the statistic is at or above the upper critical value
      lower <- max(walsh[wd_vec > cv_upper], na.rm = TRUE)
    }
    # upper bound
    if (cv_lower <= min(wd_vec, na.rm = TRUE)) {
      # edge case, set upper as largest shift available.
      upper <- max(walsh, na.rm = TRUE)
    } else {
      # smallest candidate d where the statistic is at or below the lower critical value
      upper <- min(walsh[wd_vec <= cv_lower], na.rm = TRUE)
    }
    alpha_achieved <- .alpha_lower(cv_lower - 1) + .alpha_upper(cv_upper)
  } else if (x$call$alternative == "greater") {
    # (lower, +Inf)
    cv_upper <- exactRankTests::qperm(
      p = 1 - alpha,
      scores = ranks,
      m = length(ranks),
      paired = TRUE
    )
    if (cv_upper >= max(wd_vec, na.rm = TRUE)) {
      lower <- min(walsh, na.rm = TRUE)
    } else {
      lower <- max(walsh[wd_vec > cv_upper], na.rm = TRUE)
    }
    upper <- Inf
    alpha_achieved <- .alpha_upper(cv_upper)
  } else if (x$call$alternative == "less") {
    # (-Inf, upper)
    cv_lower <- exactRankTests::qperm(
      p = alpha,
      scores = ranks,
      m = length(ranks),
      paired = TRUE
    )
    if (cv_lower <= min(wd_vec, na.rm = TRUE)) {
      upper <- max(walsh, na.rm = TRUE)
    } else {
      upper <- min(walsh[wd_vec <= cv_lower], na.rm = TRUE)
    }
    lower <- -Inf
    alpha_achieved <- .alpha_lower(cv_lower - 1)
  } else {
    msg <- "Argument 'alternative' must be one of 'two.sided', 'greater', or 'less'."
    stop(msg)
  }

  #-----------------------------------------------------------------------------
  # Hodges-Lehmann estimator via mid-quantile under exact distribution
  #-----------------------------------------------------------------------------
  min_diffs <- min(diffs)
  max_diffs <- max(diffs)
  # Zero variance data has already been handled. So can immediately calculate.
  wmean <- sum(ranks) / 2
  pseudomedian <- fmean(c(
    min(walsh[wd_vec <= ceiling(wmean)], na.rm = TRUE),
    max(walsh[wd_vec > wmean], na.rm = TRUE)
  ))
  pseudomedian_method <- "Midpoint around expected signed-rank statistic"
  conf_method <- "Inversion of exact Wilcoxon signed-rank test"

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  out <- list(
    pseudomedian = pseudomedian,
    lower = lower,
    upper = upper,
    conf_level_achieved = 1 - alpha_achieved,
    pseudomedian_method = pseudomedian_method,
    conf_method = conf_method
  )

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  c(x, out)
}

#' @importFrom stats uniroot qnorm
#' @keywords internal
#' @export
srt_pmedian.asymptotic <- function(x, ...) {
  #-----------------------------------------------------------------------------
  # `x$diffs + x$call$mu` provides mu-shift-free differences. However, this may
  # introduce zeros (that 'wilcoxon' would have otherwise removed) and/or tied
  # ranks. I'll follow wilcox.exact approach and ignore this.
  #-----------------------------------------------------------------------------
  diffs <- if (x$call$mu != 0) x$diffs + x$call$mu else x$diffs

  #-----------------------------------------------------------------------------
  # Main
  #-----------------------------------------------------------------------------
  tol_root <- x$call$tol_root
  correct <- x$call$correct
  alternative <- x$call$alternative
  zero_method <- x$call$zero_method
  digits_rank <- x$call$digits_rank

  mumin <- min(diffs)
  mumax <- max(diffs)
  # Zero variance data has already been handled. So can immediately calculate.
  pseudomedian <- uniroot(
    f = .zd,
    lower = mumin,
    upper = mumax,
    tol = tol_root,
    x = diffs,
    correct = correct,
    alternative = alternative,
    zero_method = zero_method,
    digits_rank = digits_rank
  )$root

  # The result of evaluating the signed-rank z-statistic at the lower bound of
  # the search interval (mumin).
  Zmumin <- .zd(
    d = mumin,
    x = diffs,
    correct = correct,
    alternative = alternative,
    zero_method = zero_method,
    digits_rank = digits_rank
  )

  # The result of evaluating the signed-rank z-statistic at the upper bound of
  # the search interval (mumax).
  Zmumax <- if (!is.finite(Zmumin)) {
    NA
  } else {
    .zd(
      d = mumax,
      x = diffs,
      correct = correct,
      alternative = alternative,
      zero_method = zero_method,
      digits_rank = digits_rank
    )
  }

  # Confidence interval
  conf_level <- x$call$conf_level
  alpha <- 1 - conf_level
  conf_level_achieved <- 1 - alpha

  if (!is.finite(Zmumax)) {
    cint <- c(
      if (alternative == "less") -Inf else NaN,
      if (alternative == "greater") +Inf else NaN
    )
    conf_level_achieved <- 0
    pseudomedian <- (mumin + mumax) / 2
    warning("Requested 'conf_level' not achievable.")
  } else {
    cint <- switch(
      alternative,
      two.sided = {
        repeat {
          mindiff <- Zmumin - qnorm(alpha / 2, lower.tail = FALSE)
          maxdiff <- Zmumax - qnorm(alpha / 2)
          if (mindiff < 0 || maxdiff > 0) {
            alpha <- alpha * 2
          } else {
            break
          }
        }
        if (alpha >= 1 || 1 - conf_level < alpha * 0.75) {
          conf_level_achieved <- 1 - pmin(1, alpha)
          warning("Requested 'conf_level' not achievable.")
        }
        if (alpha < 1) {
          l <- .root(
            zq = qnorm(alpha / 2, lower.tail = FALSE),
            x = diffs,
            mumin = mumin,
            mumax = mumax,
            Zmumin = Zmumin,
            Zmumax = Zmumax,
            correct = correct,
            alternative = alternative,
            zero_method = zero_method,
            digits_rank = digits_rank,
            tol_root = tol_root
          )
          u <- .root(
            zq = qnorm(alpha / 2),
            x = diffs,
            mumin = mumin,
            mumax = mumax,
            Zmumin = Zmumin,
            Zmumax = Zmumax,
            correct = correct,
            alternative = alternative,
            zero_method = zero_method,
            digits_rank = digits_rank,
            tol_root = tol_root
          )
          c(l, u)
        } else {
          rep(fmedian(diffs), 2)
        }
      },
      greater = {
        repeat {
          mindiff <- Zmumin - qnorm(alpha, lower.tail = FALSE)
          if (mindiff < 0) {
            alpha <- alpha * 2
          } else {
            break
          }
        }
        if (alpha >= 1 || 1 - conf_level < alpha * 0.75) {
          conf_level_achieved <- 1 - pmin(1, alpha)
          warning("Requested 'conf_level' not achievable.")
        }
        l <- if (alpha < 1) {
          .root(
            zq = qnorm(alpha, lower.tail = FALSE),
            x = diffs,
            mumin = mumin,
            mumax = mumax,
            Zmumin = Zmumin,
            Zmumax = Zmumax,
            correct = correct,
            alternative = alternative,
            zero_method = zero_method,
            digits_rank = digits_rank,
            tol_root = tol_root
          )
        } else {
          fmedian(diffs)
        }
        c(l, Inf)
      },
      less = {
        repeat {
          maxdiff <- Zmumax - qnorm(alpha / 2)
          if (maxdiff > 0) {
            alpha <- alpha * 2
          } else {
            break
          }
        }
        if (alpha >= 1 || 1 - conf_level < alpha * 0.75) {
          conf_level_achieved <- 1 - pmin(1, alpha)
          warning("Requested 'conf_level' not achievable.")
        }
        u <- if (alpha < 1) {
          .root(
            zq = qnorm(alpha),
            x = diffs,
            mumin = mumin,
            mumax = mumax,
            Zmumin = Zmumin,
            Zmumax = Zmumax,
            correct = correct,
            alternative = alternative,
            zero_method = zero_method,
            digits_rank = digits_rank,
            tol_root = tol_root
          )
        } else {
          fmedian(diffs)
        }
        c(-Inf, u)
      }
    )
  }

  pseudomedian_method <- "Root of standardized signed-rank statistic"
  conf_method <- "Inversion of asymptotic Wilcoxon signed-rank test"

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  out <- list(
    pseudomedian = pseudomedian,
    lower = cint[1L],
    upper = cint[2L],
    conf_level_achieved = conf_level_achieved,
    pseudomedian_method = pseudomedian_method,
    conf_method = conf_method
  )

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  c(x, out)
}

#' @keywords internal
#' @export
srt_pmedian.bootstrap <- function(x, ...) {
  #-----------------------------------------------------------------------------
  # Percentile bootstrap
  # Use original data with possible zeros that would have otherwise been removed
  # by 'wilcoxon' zero handling.
  #-----------------------------------------------------------------------------
  res <- Hmisc::pMedian(
    x = x$diffs_original,
    na.rm = TRUE,
    conf.int = x$call$conf_level,
    B = x$call$n_resamples,
    type = x$call$conf_method
  )
  conf_method <- switch(
    x$call$conf_method,
    "percentile" = "Percentile bootstrap",
    "bca" = "BCa bootstrap"
  )
  out <- list(
    pseudomedian = res[["estimate"]],
    lower = res[["lower"]],
    upper = res[["upper"]],
    pseudomedian_method = "Hodges-Lehmann",
    conf_method = conf_method,
    conf_level_achieved = NULL
  )
  return(c(x, out))
}

#-----------------------------------------------------------------------------
# Helper functions for asymptotic confidence intervals
#-----------------------------------------------------------------------------
#' @title
#' Standardized signed-rank statistic for pseudomedian CI
#'
#' @description
#' The `.zd()` function computes the standardized signed-rank statistic used to
#' estimate the pseudomedian and construct its confidence interval in a
#' one-sample Wilcoxon-type test.
#'
#' @param d
#' (Scalar numeric)\cr
#' Candidate shift at which to evaluate the statistic (i.e., test `median(x) = d`).
#'
#' @param x
#' (numeric)\cr
#' numeric observations.
#'
#' @param correct
#' (Scalar logical)\cr
#' If `TRUE`, apply a continuity correction of \eqn{\pm 0.5} on the \eqn{T^+} scale before standardization.
#'
#' @param alternative
#' (Scalar character)\cr
#' One of `"two.sided"`, `"greater"`, or `"less"`.
#' Controls the sign of the continuity correction when `correct = TRUE`.
#'
#' @param zero_method
#' (Scalar character)\cr
#' One of `"wilcoxon"` or `"pratt"`.
#' Specifies how zero differences are handled.
#'
#' @param digits_rank
#' (Scalar numeric)\cr
#' If finite, differences are rounded with `signif(., digits_rank)` to stabilize ties.
#' If `Inf`, no rounding is applied.
#'
#' @return
#' Scalar numeric z-score.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(20, mean = 0.2)
#'
#' # Evaluate z at d = 0 under Wilcoxon zero handling
#' rankdifferencetest:::.zd(
#'   d = 0,
#'   x = x,
#'   correct = TRUE,
#'   alternative = "two.sided",
#'   zero_method = "wilcoxon",
#'   digits_rank = Inf
#' )
#'
#' # Evaluate z at d = median(x) under Pratt handling with rounding
#' rankdifferencetest:::.zd(
#'   d = median(x),
#'   x = x,
#'   correct = FALSE,
#'   alternative = "two.sided",
#'   zero_method = "pratt",
#'   digits_rank = 6
#' )
#' @noRd
.zd <- function(d, x, correct, alternative, zero_method, digits_rank) {
  xd <- x - d

  list(diffs = xd) |>
    srt_ranks(
      call = list(
        mu = 0,
        correct = correct,
        alternative = alternative,
        zero_method = zero_method,
        digits_rank = digits_rank
      ),
      warn_wd = TRUE
    ) |>
    srt_method(test_class = "asymptotic") |>
    srt_statistic(warn_zd = TRUE) |>
    getElement("statistic")
}

# The .zdzq() function computes the difference between the standardized
# signed-rank statistic (from the .zd() function) and a target z-score zq.
# It is used as the objective function in root-finding (uniroot) to determine
# the bounds of a confidence interval for the pseudomedian.
.zdzq <- function(d, zq, x, correct, alternative, zero_method, digits_rank) {
  .zd(
    d = d,
    x = x,
    correct = correct,
    alternative = alternative,
    zero_method = zero_method,
    digits_rank = digits_rank
  ) -
    zq
}

# Solves for a confidence bound by finding the value of d where the
# signed-rank z-statistic equals a target quantile zq. This is done by by
# calling uniroot() on the function .zdzq(d) = .zd(d, ...) - zq, over the
# interval [mumin, mumax], using the precomputed endpoint values Zmumin and
# Zmumax to avoid redundant recomputation and to set f.lower/f.upper
# efficiently. i.e. a bracketed root-finder specialized for locating CI
# endpoints for the pseudomedian.
.root <- function(
  zq,
  mumin,
  mumax,
  Zmumin,
  Zmumax,
  x,
  correct,
  alternative,
  zero_method,
  digits_rank,
  tol_root
) {
  uniroot(
    f = .zdzq,
    lower = mumin,
    upper = mumax,
    f.lower = Zmumin - zq,
    f.upper = Zmumax - zq,
    tol = tol_root,
    zq = zq,
    x = x,
    correct = correct,
    alternative = alternative,
    zero_method = zero_method,
    digits_rank = digits_rank
  )$root
}
