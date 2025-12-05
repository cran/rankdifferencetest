#' @title
#' Ranks for the signed-rank test.
#'
#' @description
#' Computes ranks of the absolute differences and the signed-rank test statistic \eqn{W^+}.
#' The returned list is designed to be reused by higher-level signed-rank functions.
#'
#' @details
#' Consider a numeric vector of paired differences \eqn{x = (x_1, \dots, x_n)}.
#' After removing `NA` values, let \eqn{\mathcal{I}_0 = \{i : x_i = 0\}} and \eqn{\mathcal{I}_{\neq 0} = \{i : x_i \neq 0\}}.
#'
#' ## Zero handling
#'
#' For the Wilcoxon method (`zero_method = "wilcoxon"`), zeros are removed prior to ranking, so the ranking is performed on \eqn{\{x_i : i \in \mathcal{I}_{\neq 0}\}} only.
#'
#' For the Pratt method (`zero_method = "pratt"`), zeros are retained when computing ranks, but their corresponding ranks do not contribute to the signed-rank sum \eqn{W^+}.
#'
#' ## Ranking
#'
#' Ranks are assigned to \eqn{|x_i|} using average-tie handling.
#' The argument `digits_rank` controls rounding for ranking only.
#' Ranks are computed from \eqn{|\mathrm{signif}(x_i, \text{digits\_rank})|} when `digits_rank` is finite, and from \eqn{|x_i|} otherwise.
#' This rounding may induce ties, which changes both the values of the averaged ranks and the variance of the statistic in asymptotic procedures.
#'
#' The next function in the 'srt pipeline', [srt_method()], calculates the number of ties among the ranks (`n_ties`).
#'
#' ## Wilcoxon signed-Rank statistic
#'
#' Let \eqn{r_i} denote the absolute-value rank of the \eqn{i}-th observation after applying the chosen zero handling and ranking precision.
#' The sum of positive ranks is
#'
#' \deqn{W^+ = \sum_{i : x_i > 0} r_i,}
#'
#' which is the canonical Wilcoxon signed-rank statistic used for both exact and asymptotic inference.
#'
#' @param x
#' (named list)\cr
#' A named list which contains a vector of paired differences named `"diffs"`.
#' Typically the object returned by [srt_data()].
#'
#' @param call
#' (call or named list)\cr
#' A call or named list which contains the arguments from the parent function [srt()] or [rdt()].
#' The minimum list of arguments required for `srt_ranks()` are:
#' - `mu`
#'   - Scalar numeric `(-Inf, Inf)`.
#'     Under the null hypothesis, `diffs` is assumed to be symmetric around `mu`.
#' - `zero_method`
#'   - Scalar character `c("wilcoxon", "pratt")`.
#'     String for zero handling.
#' - `digits_rank`
#'   - Scalar numeric `(0, Inf]`.
#'     Controls ranking precision.
#'
#' Next in line, [srt_method()] needs:
#' - `[[1L]]`
#'   - `as.name()` of calling function
#'   - One of `quote(srt)`, `quote(srt2)`, `quote(rdt)`, or `quote(rdt2)`
#' - `distribution`
#'   - Scalar character: `c("auto", "exact", "asymptotic")`
#'   - The method used to calculate the p-value.
#'
#' @param warn_wd
#' (Scalar logical: `c(FALSE, TRUE)`)\cr
#' Used for
#' 1. Midpoint algorithm 'W+(d)' in confidence interval test inversion.
#' 2. Recalculating ranks after removing mu-shift.
#' If `TRUE` return a warning and exit with value `wplus = 0`.
#'
#' @returns
#' list
#'
#' @examples
#' library(rankdifferencetest)
#'
#' # Synthetic paired differences with zeros and ties
#' set.seed(1)
#' diffs <- c(rnorm(8, mean = 0.3), 0, 0, 0, round(rnorm(8, mean = -0.2), 1))
#' x <- list(diffs = diffs)
#'
#' # Wilcoxon zero method: zeros dropped before ranking
#' call <- list(mu = 0, zero_method = "wilcoxon", digits_rank = Inf)
#' cw <- rankdifferencetest:::srt_ranks(x, call)
#' cw$ranks
#' cw$wplus
#'
#' # Pratt zero method: zeros retained for ranking
#' call <- list(mu = 0, zero_method = "pratt", digits_rank = Inf)
#' cp <- rankdifferencetest:::srt_ranks(x, call)
#' cp$wplus
#' cp$n_signed
#'
#' # Induce ties via ranking precision
#' call <- list(mu = 0, zero_method = "wilcoxon", digits_rank = 1)
#' ctied <- rankdifferencetest:::srt_ranks(x, call)
#' ctied$ranks
#'
#' @keywords internal
srt_ranks <- function(
  x,
  call,
  warn_wd = FALSE
) {
  if (call$mu != 0) {
    x$diffs <- x$diffs - call$mu
  }
  diffs <- x$diffs

  n_zeros <- sum(diffs == 0)
  n_signed <- length(diffs) - n_zeros

  if (warn_wd) {
    if (n_signed == 0L) {
      msg <- "When evaluating W+(d), all differences are zero at chosen d.\nReturning W+(d)=0"
      warning(msg)
      return(c(x, list(call = call), wplus = 0))
    }
  }

  #-----------------------------------------------------------------------------
  # Ranks
  #-----------------------------------------------------------------------------
  are_zeros <- diffs == 0
  any_zeros <- any(are_zeros)

  if (any_zeros && call$zero_method == "wilcoxon") {
    diffs <- diffs[!are_zeros]
  }

  # Rounding only affects ranking
  abs_diffs <- if (is.finite(call$digits_rank)) {
    abs(signif(diffs, call$digits_rank))
  } else {
    abs(diffs)
  }

  # Ranks with average ties
  ranks <- rank(abs_diffs, na.last = NA, ties.method = "average")

  #-----------------------------------------------------------------------------
  # Signed-rank statistic W+
  #-----------------------------------------------------------------------------
  wplus <- sum(ranks[diffs > 0])

  # Reduce size after wplus (ensures length(ranks) == length(diffs))
  if (call$zero_method == "pratt") {
    ranks <- ranks[!are_zeros]
  }

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  # Update x
  if (any_zeros && call$zero_method == "wilcoxon") {
    x$diffs <- diffs
  }

  out <- list(
    n_zeros = n_zeros,
    n_signed = n_signed,
    ranks = ranks,
    wplus = wplus
  )

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  c(x, list(call = call), out)
}
