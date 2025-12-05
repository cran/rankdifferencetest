#' Rank difference test
#'
#' Performs Kornbrot's rank difference test.
#' The rank difference test is a modified Wilcoxon signed-rank test that produces consistent and meaningful results for ordinal or monotonically transformed data.
#'
#' For paired data, the Wilcoxon signed-rank test results in subtraction of the paired values.
#' However, this subtraction is not meaningful for ordinal scale variables.
#' In addition, any monotone transformation of the data will result in different signed ranks, thus different p-values.
#' However, ranking the original data allows for meaningful addition and subtraction of ranks and preserves ranks over monotonic transformation.
#' Kornbrot developed the rank difference test for these reasons.
#'
#' Kornbrot recommends that the rank difference test be used in preference to the Wilcoxon signed-rank test in all paired comparison designs where the data are not both of interval scale and of known distribution.
#' The rank difference test preserves good power compared to Wilcoxon's signed-rank test, is more powerful than the sign test, and has the benefit of being a true distribution-free test.
#'
#' The procedure for Kornbrot's rank difference test is as follows:
#'
#' 1. Combine all \eqn{2n} paired observations.
#' 2. Order the values from smallest to largest.
#' 3. Assign ranks \eqn{1, 2, \dots, 2n} with average rank for ties.
#' 4. Perform the Wilcoxon signed-rank test using the paired ranks.
#'
#' The test statistic for the rank difference test \eqn{(D)} is not exactly equal to the test statistic of the naive rank-transformed Wilcoxon signed-rank test \eqn{(W^+)}.
#' However, using \eqn{W^+} should result in a conservative estimate for \eqn{D}, and they approach in distribution as the sample size increases.
#' \insertCite{kornbrot1990;textual}{rankdifferencetest} discusses methods for calculating \eqn{D} when \eqn{n<7} and \eqn{8 < n \leq 20}.
#' [rankdifferencetest::rdt()] uses \eqn{W^+} instead of \eqn{D}.
#'
#' See [rankdifferencetest::srt()] for additional details about implementation of Wilcoxon's signed-rank test.
#'
#' @references
#' \insertRef{kornbrot1990}{rankdifferencetest}
#'
#' @param x
#' (numeric)\cr
#' Numeric vector of data.
#' Differences of ranks correspond with `x - y`.
#' Pairs with missing values are silently dropped.
#'
#' @param y
#' (numeric)\cr
#' Numeric vector of data.
#' Differences of ranks correspond with `x - y`.
#' Pairs with missing values are silently dropped.
#'
#' @param data
#' (data.frame)\cr
#' The data frame of interest.
#'
#' @param formula
#' (formula)\cr
#' A formula of form:
#' \describe{
#'   \item{y ~ group | block}{
#'     Use when `data` is in tall format.
#'     `y` is the numeric outcome, `group` is the binary grouping variable, and `block` is the subject/item-level variable indicating pairs of observations.
#'     `group` will be converted to a factor and the first level will be the reference value.
#'     For example, when `levels(data$group) <- c("pre", "post")`, the focal level is 'post', so differences are `post - pre`.
#'     Pairs with missing values are silently dropped.
#'     See `agg_fun` for handling duplicate cases of grouping/blocking combinations.
#'   }
#'   \item{y ~ x}{
#'     Use when `data` is in wide format.
#'     `y` and `x` must be numeric vectors.
#'     Differences of ranks correspond with `data$y - data$x`.
#'     Pairs with missing values are silently dropped.
#'   }
#' }
#'
#' @param conf_level
#' (Scalar numeric: `[0, 1)`)\cr
#' The confidence level.
#' Typically `0.95`.
#' If `0` (default), no confidence interval is calculated.
#'
#' @param conf_method
#' (Scalar character: `c("inversion", "percentile", "bca")`)\cr
#' The type of confidence interval.
#' If `"inversion"` (default), the bounds are computed by inverting the hypothesis test.
#' If `"percentile"`, the bounds are computed using a percentile bootstrap.
#' If `"bca"`, the bounds are computed using a bias-corrected and accelerated (BCa) bootstrap.
#'
#' @param n_resamples
#' (Scalar integer: `1000L`; `[10, Inf)`)\cr
#' The number of bootstrap resamples.
#' Only used if `"percentile"` or `"bca"` confidence intervals are chosen.
#'
#' @param alternative
#' (Scalar character: `c("two.sided", "greater", "less")`)\cr
#' The alternative hypothesis.
#' Must be one of `"two.sided"` (default), `"greater"`, or `"less"`.
#'
#' @param mu
#' (Scalar numeric: `0`; `(-Inf, Inf)`)\cr
#' Under the null hypothesis, differences of ranks are assumed to be symmetric around `mu`.
#'
#' @param distribution
#' (Scalar character: `c("auto", "exact", "asymptotic")`)\cr
#' The method used to calculate the p-value.
#' If `"auto"` (default), an appropriate method will automatically be chosen (`distribution = "exact"` when n < 50 or `distribution = "asymptotic"` otherwise).
#' If `"exact"`, the exact Wilcoxon signed-rank distribution is used.
#' If `"asymptotic"`, the asymptotic normal approximation is used.
#'
#' @param correct
#' (Scalar logical: `c(TRUE, FALSE)`)\cr
#' Whether or not to apply a continuity correction to the Z-statistic for the asymptotic approximation of the p-value.
#'
#' @param zero_method
#' (Scalar character: `c("wilcoxon", "pratt")`)\cr
#' The method used to handle differences of ranks equal to zero.
#' If `"wilcoxon"` (default), zeros are removed prior to ranking (classic Wilcoxon convention).
#' If `"pratt"`, zeros are retained for ranking but excluded from the signed-rank sum.
#'
#' @param agg_fun
#' (Scalar character or function: `"error"`)\cr
#' Used for aggregating duplicate cases of grouping/blocking combinations when data is in tall format and `formula` has structure `y ~ group | block`.
#' `"error"` (default) will return an error if duplicate grouping/blocking combinations are encountered.
#' Select one of `"first"`, `"last"`, `"sum"`, `"mean"`, `"median"`, `"min"`, or `"max"` for built in aggregation handling (each applies `na.rm = TRUE`).
#' Or define your own function.
#' For example, `myfun <- function(x) {as.numeric(quantile(x, 0.75, na.rm = TRUE))}`.
#'
#' @param digits_rank
#' (Scalar integer: `Inf`; `(0, Inf]`)\cr
#' Controls ranking precision.
#' If finite, ranks are computed from [base::signif]`(abs(diffs), digits_rank)`.
#' If `Inf` (default), ranks are computed from `abs(diffs)`.
#' Smaller values may introduce ties (because they no longer depend on extremely small numeric differences) and thus change averaged ranks and tie counts.
#'
#' @param tol_root
#' (Scalar numeric: `1e-4`; `(0, Inf)`)\cr
#' For [stats::uniroot]`(tol=tol_root)` calls when `conf_level > 0` and `distribution = "asymptotic"`.
#'
#' @seealso
#' [rankdifferencetest::srt()]
#'
#' @returns
#' A list with the following elements:
#' \tabular{llll}{
#'   Slot \tab Subslot \tab Name \tab Description \cr
#'   1 \tab \tab `p_value`       \tab p-value. \cr
#'   2 \tab \tab `statistic`     \tab Test statistic. \eqn{W^+} for the exact Wilcoxon signed-rank distribution. \eqn{Z} for the asymptotic normal approximation. \cr
#'   3 \tab \tab `pseudomedian`  \tab Measure of centrality. \cr
#'   4 \tab \tab `lower`         \tab Lower bound of confidence interval for the pseudomedian. `NULL` if no CI requested. \cr
#'   5 \tab \tab `upper`         \tab Upper bound of confidence interval for the pseudomedian. `NULL` if no CI requested. \cr
#'   6 \tab \tab `method`        \tab Test method. \cr
#'   7 \tab   \tab `info`        \tab Additional test information. \cr
#'   7 \tab 1 \tab `p_value_method` \tab Method used to calculate the p-value. \cr
#'   7 \tab 2 \tab `pseudomedian_method` \tab Method used to calculate the pseudomedian. \cr
#'   7 \tab 3 \tab `conf_method` \tab Method used to calculate the confidence interval. \cr
#'   7 \tab 4 \tab `conf_level_achieved` \tab Achieved confidence level. \cr
#'   7 \tab 5 \tab `n_sample`    \tab Number of observations in the original data. \cr
#'   7 \tab 6 \tab `n_analytic`  \tab Number of observations after removing missing values from the original data. \cr
#'   7 \tab 7 \tab `n_zeros`     \tab Number of zeros among differences of ranks in the analytic data set. \cr
#'   7 \tab 8 \tab `n_signed`    \tab Number of nonzero differences of ranks in the analytic data set. \cr
#'   7 \tab 9 \tab `n_ties`      \tab Number of tied ranks after ranking the absolute differences of ranks. \cr
#'   7 \tab 10 \tab `data_type`  \tab Data type. \cr
#'   7 \tab 11 \tab `focal_name` \tab Name of the focal variable (differences are focal - reference). \cr
#'   7 \tab 12 \tab `reference_name` \tab Name of the reference variable (differences are focal - reference). \cr
#'   8 \tab   \tab `call`       \tab A named list of the function's arguments (use `as.call()` to convert to a call; `call$distribution` may be updated from `"exact"` to `"asymptotic"`).
#' }
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # rdt() example
#' #----------------------------------------------------------------------------
#' library(rankdifferencetest)
#'
#' # Use example data from Kornbrot (1990)
#' data <- kornbrot_table1
#'
#' # Create tall-format data for demonstration purposes
#' data_tall <- reshape(
#'   data = kornbrot_table1,
#'   direction = "long",
#'   varying = c("placebo", "drug"),
#'   v.names = c("time"),
#'   idvar = "subject",
#'   times = c("placebo", "drug"),
#'   timevar = "treatment",
#'   new.row.names = seq_len(prod(length(c("placebo", "drug")), nrow(kornbrot_table1)))
#' )
#'
#' # Subject and treatment should be factors. The ordering of the treatment factor
#' # will determine the difference (placebo - drug).
#' data_tall$subject <- factor(data_tall$subject)
#' data_tall$treatment <- factor(data_tall$treatment, levels = c("drug", "placebo"))
#'
#' # Recreate analysis and results from table 3 (page 248) in Kornbrot (1990)
#' ## Divide p-value by 2 for one-tailed probability.
#' rdt(
#'   data = data,
#'   formula = placebo ~ drug,
#'   alternative = "two.sided",
#'   distribution = "asymptotic",
#'   zero_method = "wilcoxon",
#'   correct = TRUE,
#'   conf_level = 0.95
#' )
#'
#' rdt2(
#'   x = data$placebo,
#'   y = data$drug,
#'   alternative = "two.sided",
#'   distribution = "asymptotic",
#'   zero_method = "wilcoxon",
#'   correct = TRUE,
#'   conf_level = 0.95
#' )
#'
#' # The same outcome is seen after transforming time to rate.
#' ## Rate transformation inverts the rank ordering.
#' data$placebo_rate <- 60 / data$placebo
#' data$drug_rate <- 60 / data$drug
#' data_tall$rate <- 60 / data_tall$time
#'
#' rdt(
#'   data = data_tall,
#'   formula = rate ~ treatment | subject,
#'   alternative = "two.sided",
#'   distribution = "asymptotic",
#'   zero_method = "wilcoxon",
#'   correct = TRUE,
#'   conf_level = 0.95
#' )
#'
#' # In contrast to the rank difference test, the Wilcoxon signed-rank test
#' # produces differing results. See table 1 and table 2 (page 245) in
#' # Kornbrot (1990).
#' ## Divide p-value by 2 for one-tailed probability.
#' srt(
#'   data = data,
#'   formula = placebo ~ drug,
#'   alternative = "two.sided",
#'   distribution = "asymptotic",
#'   zero_method = "wilcoxon",
#'   correct = TRUE,
#'   conf_level = 0.95
#' )
#'
#' srt(
#'   data = data_tall,
#'   formula = rate ~ treatment | subject,
#'   alternative = "two.sided",
#'   distribution = "asymptotic",
#'   zero_method = "wilcoxon",
#'   correct = TRUE,
#'   conf_level = 0.95
#' )
#'
#' @rdname rdt
#' @export
rdt <- function(
  data,
  formula,
  conf_level = 0,
  conf_method = "inversion",
  n_resamples = 1000L,
  alternative = "two.sided",
  mu = 0,
  distribution = "auto",
  correct = TRUE,
  zero_method = "wilcoxon",
  agg_fun = "error",
  digits_rank = Inf,
  tol_root = 1e-4
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be an object of class 'data.frame'.")
  }
  if (!is_formula(formula)) {
    stop("Argument 'formula' must be an object of class 'formula'.")
  }

  if (
    !(is.finite(conf_level) &&
      is.numeric(conf_level) &&
      length(conf_level) == 1L &&
      conf_level >= 0 &&
      conf_level < 1)
  ) {
    stop("Argument 'conf_level' must be a number between 0 and 1.")
  }
  if (
    !(is.character(conf_method) &&
      length(conf_method) == 1L &&
      conf_method %in% c("inversion", "percentile", "bca"))
  ) {
    msg <- "Argument 'conf_method' must be one of 'inversion', 'percentile', or 'bca'."
    stop(msg)
  }
  if (
    !(is.finite(n_resamples) &&
      is.numeric(n_resamples) &&
      length(n_resamples) == 1L &&
      n_resamples >= 10L)
  ) {
    stop(
      "Argument 'n_resamples' must be a scalar integer greater than or equal to 10."
    )
  }
  n_resamples <- as.integer(n_resamples)

  if (
    !(is.character(alternative) &&
      length(alternative) == 1L &&
      alternative %in% c("two.sided", "greater", "less"))
  ) {
    msg <- "Argument 'alternative' must be one of 'two.sided', 'greater', or 'less'."
    stop(msg)
  }
  if (!(is.finite(mu) && is.numeric(mu) && length(mu) == 1L)) {
    stop("Argument 'mu' must be a scalar vector of class 'numeric'.")
  }
  if (
    !(is.character(distribution) &&
      length(distribution) == 1L &&
      distribution %in% c("auto", "exact", "asymptotic"))
  ) {
    msg <- "Argument 'distribution' must be one of 'auto', 'exact', or 'asymptotic'."
    stop(msg)
  }
  if (!(!is.na(correct) && is.logical(correct) && length(correct) == 1L)) {
    stop("Argument 'correct' must be a scalar vector of class 'logical'.")
  }
  if (
    !(is.character(zero_method) &&
      length(zero_method) == 1L &&
      zero_method %in% c("pratt", "wilcoxon"))
  ) {
    stop("Argument 'zero_method' must be one of 'wilcoxon' or 'pratt'.")
  }

  if (is.character(agg_fun)) {
    if (
      !(length(agg_fun) == 1L &&
        agg_fun %in%
          c("error", "first", "last", "sum", "mean", "median", "min", "max"))
    ) {
      msg <- "Argument 'agg_fun' must be one of 'error', 'first', 'last', 'sum', 'mean', 'median', 'min', or 'max'."
      stop(msg)
    }
  } else {
    if (!is.function(agg_fun)) {
      stop("Argument 'agg_fun' must be a function or scalar character.")
    }
  }

  if (
    !(!is.na(digits_rank) &&
      is.numeric(digits_rank) &&
      length(digits_rank) == 1L &&
      digits_rank > 0)
  ) {
    stop("Argument 'digits_rank' must be a scalar integer > 0 (may be Inf).")
  }
  if (
    !(is.finite(tol_root) &&
      is.numeric(tol_root) &&
      length(tol_root) == 1L &&
      tol_root > 0)
  ) {
    stop("Argument 'tol_root' must be a scalar numeric > 0.")
  }

  # match.call() isn't compatible with passing values through pipeline method.
  # So build manually, can as.call() later.
  call <- list(
    quote(rdt),
    data = substitute(data),
    formula = formula,
    conf_level = conf_level,
    conf_method = conf_method,
    n_resamples = n_resamples,
    alternative = alternative,
    mu = mu,
    distribution = distribution,
    correct = correct,
    zero_method = zero_method,
    agg_fun = agg_fun,
    digits_rank = digits_rank,
    tol_root = tol_root
  )

  #-----------------------------------------------------------------------------
  # Prepare results
  #-----------------------------------------------------------------------------
  res <- rdt_data(x = data, formula = formula, agg_fun = agg_fun) |>
    srt_ranks(call = call) |>
    srt_method() |>
    srt_statistic() |>
    srt_pvalue() |>
    srt_pmedian() |>
    srt_result()

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

#' @rdname rdt
#' @export
rdt2 <- function(
  x,
  y,
  conf_level = 0,
  conf_method = "inversion",
  n_resamples = 1000L,
  alternative = "two.sided",
  mu = 0,
  distribution = "auto",
  correct = TRUE,
  zero_method = "wilcoxon",
  digits_rank = Inf,
  tol_root = 1e-4
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.numeric(x)) {
    stop("Argument 'x' must be an object of class 'numeric'.")
  }
  if (!is.numeric(y)) {
    stop("Argument 'y' must be an object of class 'numeric'.")
  }
  if (length(x) != length(y)) {
    stop("Arguments 'x' and 'y' must be the same length.")
  }
  focal_name <- deparse1(substitute(x))
  reference_name <- deparse1(substitute(y))

  if (
    !(is.finite(conf_level) &&
      is.numeric(conf_level) &&
      length(conf_level) == 1L &&
      conf_level >= 0 &&
      conf_level < 1)
  ) {
    stop("Argument 'conf_level' must be a number between 0 and 1.")
  }
  if (
    !(is.character(conf_method) &&
      length(conf_method) == 1L &&
      conf_method %in% c("inversion", "percentile", "bca"))
  ) {
    msg <- "Argument 'conf_method' must be one of 'inversion', 'percentile', or 'bca'."
    stop(msg)
  }
  if (
    !(is.finite(n_resamples) &&
      is.numeric(n_resamples) &&
      length(n_resamples) == 1L &&
      n_resamples >= 10L)
  ) {
    stop(
      "Argument 'n_resamples' must be a scalar integer greater than or equal to 10."
    )
  }
  n_resamples <- as.integer(n_resamples)

  if (
    !(is.character(alternative) &&
      length(alternative) == 1L &&
      alternative %in% c("two.sided", "greater", "less"))
  ) {
    msg <- "Argument 'alternative' must be one of 'two.sided', 'greater', or 'less'."
    stop(msg)
  }
  if (!(is.finite(mu) && is.numeric(mu) && length(mu) == 1L)) {
    stop("Argument 'mu' must be a scalar vector of class 'numeric'.")
  }
  if (
    !(is.character(distribution) &&
      length(distribution) == 1L &&
      distribution %in% c("auto", "exact", "asymptotic"))
  ) {
    msg <- "Argument 'distribution' must be one of 'auto', 'exact', or 'asymptotic'."
    stop(msg)
  }
  if (!(!is.na(correct) && is.logical(correct) && length(correct) == 1L)) {
    stop("Argument 'correct' must be a scalar vector of class 'logical'.")
  }
  if (
    !(is.character(zero_method) &&
      length(zero_method) == 1L &&
      zero_method %in% c("pratt", "wilcoxon"))
  ) {
    stop("Argument 'zero_method' must be one of 'wilcoxon' or 'pratt'.")
  }

  if (
    !(!is.na(digits_rank) &&
      is.numeric(digits_rank) &&
      length(digits_rank) == 1L &&
      digits_rank > 0)
  ) {
    stop("Argument 'digits_rank' must be a scalar integer > 0 (may be Inf).")
  }
  if (
    !(is.finite(tol_root) &&
      is.numeric(tol_root) &&
      length(tol_root) == 1L &&
      tol_root > 0)
  ) {
    stop("Argument 'tol_root' must be a scalar numeric > 0.")
  }

  # match.call() isn't compatible with passing values through pipeline method.
  # So build manually, can as.call() later.
  call <- list(
    quote(rdt2),
    x = substitute(x),
    y = substitute(y),
    conf_level = conf_level,
    conf_method = conf_method,
    n_resamples = n_resamples,
    alternative = alternative,
    mu = mu,
    distribution = distribution,
    correct = correct,
    zero_method = zero_method,
    digits_rank = digits_rank,
    tol_root = tol_root
  )

  #-----------------------------------------------------------------------------
  # Prepare data
  #-----------------------------------------------------------------------------
  res <- rdt_data(x = x, y = y, x_name = focal_name, y_name = reference_name) |>
    srt_ranks(call = call) |>
    srt_method() |>
    srt_statistic() |>
    srt_pvalue() |>
    srt_pmedian() |>
    srt_result()

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}
