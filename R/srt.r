#' Signed-rank test
#'
#' Performs Wilcoxon's signed-rank test.
#'
#' The procedure for Wilcoxon's signed-rank test is as follows:
#'
#' 1. For one-sample data `x` or paired samples `x` and `y`, where `mu` is the measure of center under the null hypothesis, define the values used for analysis as `(x - mu)` or `(x - y - mu)`.
#' 2. Define 'zero' values as `(x - mu == 0)` or `(x - y - mu == 0)`.
#'    - `zero_method = "wilcoxon"`: Remove values equal to zero.
#'    - `zero_method = "pratt"`: Keep values equal to zero.
#' 3. Order the absolute values from smallest to largest.
#' 4. Assign ranks \eqn{1, 2, \dots, n} to the ordered absolute values, using mean rank for ties.
#'    - `zero_method = "pratt"`: remove values equal to zero and their corresponding ranks.
#' 5. Calculate the test statistic.
#'
#'    Calculate \eqn{W^+} as the sum of the ranks for positive values.
#'    Let \eqn{r_i} denote the absolute-value rank of the \eqn{i}-th observation after applying the chosen zero handling and ranking precision and \eqn{v_i} denote the values used for analysis.
#'    Then
#'    \deqn{W^+ = \sum_{i : v_i > 0} r_i.}
#'    \eqn{W^+ + W^- = n(n+1)/2}, thus either can be calculated from the other.
#'    If the null hypothesis is true, \eqn{W^+} and \eqn{W^-} are expected to be equal.
#'
#'    - `distribution = "exact"`: Use \eqn{W^+} which takes values between \eqn{0} and \eqn{n(n+1)/2}.
#'    - `distribution = "asymptotic"`: Use the standardized test statistic \eqn{Z=\frac{W^+ - E_0(W^+) - cc}{Var_0(W^+)^{1/2}}} where \eqn{Z \sim N(0, 1)} asymptotically.
#'      - `zero_method = "wilcoxon"`: Under the null hypothesis, the expected mean and variance are
#'        \deqn{
#'        \begin{aligned}
#'        E_0(W^+) &= \frac{n(n+1)}{4} \\
#'        Var_0(W^+) &= \frac{n(n+1)(2n+1)}{24} - \frac{\sum(t^3-t)}{48},
#'        \end{aligned}
#'        }
#'        where \eqn{t} is the number of ties for each unique ranked absolute value.
#'      - `zero_method = "pratt"`: Under the null hypothesis, the expected mean and variance are
#'        \deqn{
#'        \begin{aligned}
#'        E_0(W^+) &= \frac{n(n+1)}{4} - \frac{n_{zeros}(n_{zeros}+1)}{4} \\
#'        Var_0(W^+) &= \frac{n(n+1)(2n+1) - n_{zeros}(n_{zeros}+1)(2n_{zeros}+1)}{24} - \frac{\sum(t^3-t)}{48},
#'        \end{aligned}
#'        }
#'        where \eqn{t} is the number of ties for each unique ranked absolute value.
#'      - `correct = TRUE`: The continuity correction is defined by:
#'        \deqn{
#'          cc = \begin{cases}
#'            0.5 \times \text{sign}(W^+ - E_0(W^+)) & \text{two-sided alternative} \\
#'            0.5 & \text{greater than alternative} \\
#'            -0.5 & \text{less than alternative.}
#'            \end{cases}
#'        }
#'      - `correct = FALSE`: Set \eqn{cc = 0}.
#'
#'    Alternatively, \eqn{E_0(W^+)} and \eqn{Var_0(W^+)} can be calculated without need for closed form expressions that include zero correction and/or tie correction.
#'    Consider each rank \eqn{r_i} (with averaged rank for ties) as a random variable \eqn{X_i} that contributes to \eqn{W^+}.
#'    Under the null hypothesis, each rank has 50% chance of being positive or negative. So \eqn{X_i} can take values
#'    \deqn{
#'      X_i = \begin{cases}
#'            r_i & \text{with probability } p = 0.5 \\
#'            0 & \text{with probability } 1 - p = 0.5.
#'            \end{cases}
#'    }
#'    Using \eqn{Var(X_i) = E(X_i^2) - E(X_i)^2} where
#'    \deqn{
#'    \begin{aligned}
#'    E(X_i) &= p \cdot r_i + (1 - p) \cdot 0 \\
#'           &= 0.5r_i \\
#'    E(X_i^2) &= p \cdot r_i^2 + (1 - p) \cdot 0^2 \\
#'             &= 0.5r_i^2 \\
#'    E(X_i)^2 &= (0.5r_i)^2,
#'    \end{aligned}
#'    }
#'    it follows that \eqn{Var(X_i) = 0.5r_i^2 - (0.5r_i)^2 = 0.25r_i^2.}
#'    Hence, \eqn{E_0(W^+) = \frac{\sum r_i}{2}} and \eqn{Var_0(W^+) = \frac{\sum r_i^2}{4}.}
#' 6. Calculate the p-value.
#'    - `distribution = "exact"`
#'      - No zeros among the differences and ranks are tie free
#'        - Use the Wilcoxon signed-rank distribution as implemented in [stats::psignrank()] to calculate the probability of being as or more extreme than \eqn{W^+.}
#'      - Zeros present and/or ties present
#'        - Use the Shift-Algorithm from \insertCite{streitberg1984;textual}{rankdifferencetest} as implemented in [exactRankTests::pperm()].
#'    - `distribution = "asymptotic"`
#'      - Use the standard normal distribution as implemented in [stats::pnorm()] to calculate the probability of being as or more extreme than \eqn{Z.}
#'
#' ## Hypotheses and test assumptions
#'
#' The signed-rank test hypotheses are stated as:
#' - Null: `(x)` or `(x - y)` is centered at `mu`.
#' - Two-sided alternative: `(x)` or `(x - y)` is not centered at `mu`.
#' - Greater than alternative: `(x)` or `(x - y)` is centered at a value greater than `mu`.
#' - Less than alternative: `(x)` or `(x - y)` is centered at a value less than `mu`.
#'
#' The signed-rank test traditionally assumes the differences are independent with identical, continuous, and symmetric distribution.
#' However, not all of these assumptions may be required \insertCite{pratt1981}{rankdifferencetest}.
#' The 'identically distributed' assumption is not required, keeping the level of test as expected for the hypotheses as stated above.
#' The symmetry assumption is not required when using one-sided alternative hypotheses.
#' For example:
#' - Null: `(x)` or `(x - y)` is symmetric and centered at `mu`.
#' - Greater than alternative: `(x)` or `(x - y)` is stochastically larger than `mu`.
#' - Less than alternative: `(x)` or `(x - y)` is stochastically smaller than `mu`.
#'
#' ## Zero handling
#'
#' `zero_method = "pratt"` uses the method by \insertCite{pratt1959;textual}{rankdifferencetest}, which first rank-transforms the absolute values, including zeros, and then removes the ranks corresponding to the zeros.
#' `zero_method = "wilcoxon"` uses the method by \insertCite{wilcoxon1950;textual}{rankdifferencetest}, which first removes the zeros and then rank-transforms the remaining absolute values.
#' \insertCite{conover1973;textual}{rankdifferencetest} found that when comparing a discrete uniform distribution to a distribution where probabilities linearly increase from left to right, Pratt's method outperforms Wilcoxon's.
#' When testing a binomial distribution centered at zero to see whether the parameter of each Bernoulli trial is \eqn{\frac{1}{2}}, Wilcoxon's method outperforms Pratt's.
#'
#' ## Pseudomedians
#'
#' ### Hodges-Lehmann estimator
#'
#' The Hodges-Lehmann estimator is the median of all pairwise averages of the sample values.
#' \deqn{\mathrm{HL} = \mathrm{median} \left\{ \frac{x_i + x_j}{2} \right\}_{i \le j}}
#' This pseudomedian is a robust, distribution-free estimate of central tendency for a single sample, or a location-shift estimator for paired data.
#' It's resistant to outliers and compatible with rank-based inference.
#' This statistic is returned when `conf_level = 0` (for all test methods) or when confidence intervals are requested for an exact Wilcoxon test with zero-free data and tie-free ranks.
#'
#' ### Midpoint around expected signed-rank statistic (exact CI-centered estimator)
#'
#' The exact test for data which contains zeros or whose ranks contain ties uses the Streitberg-Rohmel Shift-algorithm.
#' When a confidence interval is requested under this scenario, the estimated pseudomedian is a midpoint around the expected rank sum.
#' This midpoint is the average of the largest shift value \eqn{d} whose signed-rank statistic does not exceed the null expectation and the smallest \eqn{d} whose statistic exceeds it.
#'
#' In detail, let \eqn{W^+(d)} be the Wilcoxon signed-rank sum at shift \eqn{d}, and let \eqn{E_0} denote the null expectation (e.g., \eqn{\sum r_i / 2} when ranks are \eqn{r_i}).
#' Then \deqn{\hat{d}_{\mathrm{mid}} = \frac{1}{2} \left( \min\{ d : W^+(d) \le \lceil E_0 \rceil \} + \max\{ d : W^+(d) > E_0 \} \right)}
#'
#' This pseudomedian is a discrete-compatible point estimate that centers the exact confidence interval obtained by inverting the exact signed-rank distribution.
#' It may differ from the Hodges-Lehmann estimator when data are tied or rounded.
#'
#' ### Root of standardized signed-rank statistic (asymptotic CI-centered estimator)
#'
#' A similar algorithm is used to estimate the pseudomedian when a confidence interval is requested under the asymptotic test scenario.
#' This pseudomedian is the value of the shift \eqn{d} for which the standardized signed-rank statistic equals zero under the asymptotic normal approximation.
#'
#' In detail, let \eqn{W^+(d)} be the signed-rank sum, with null mean \eqn{E_0(d)} and null variance \eqn{\mathrm{Var}_0(d)} (with possible tie and continuity corrections).
#' Define \deqn{Z(d) = \frac{W^+(d) - E_0(d)}{\sqrt{\mathrm{Var}_0(d)}}.}
#' The pseudomedian is the root \deqn{\hat{d}_{\mathrm{root}} = \{ d : Z(d) = 0 \}.}
#'
#' This pseudomedian is a continuous-compatible point estimate that centers the asymptotic confidence interval.
#' It's the solution to the test-inversion equation under a normal approximation.
#' It's sensitive to tie/zero patterns through \eqn{\mathrm{Var}_0(d)}, may include a continuity correction, and is not guaranteed to equal the Hodges-Lehmann estimator or the exact midpoint.
#'
#' ## Confidence intervals
#'
#' ### Exact Wilcoxon confidence interval
#'
#' The exact Wilcoxon confidence interval is obtained by inverting the exact distribution of the signed-rank statistic.
#' It uses the permutation distribution of the Wilcoxon statistic and finds bounds where cumulative probabilities cross \eqn{\alpha/2} and \eqn{1-\alpha/2}.
#' Endpoints correspond to quantiles from [stats::qsignrank()].
#' This interval guarantees nominal coverage under the null hypothesis without relying on asymptotic approximations.
#' It respects discreteness of the data and may produce conservative intervals when the requested confidence level is not achievable (with warnings).
#'
#' ### Exact Wilcoxon confidence interval using the Shift-algorithm
#'
#' The exact Wilcoxon confidence interval using the Shift-algorithm is obtained by enumerating all candidate shifts and inverting the exact signed-rank distribution.
#' In detail, it generates all pairwise averages \eqn{\frac{x_i + x_j}{2}}, evaluates the signed-rank statistic for each candidate shift, and determines bounds using [exactRankTests::pperm()] and [exactRankTests::qperm()].
#'
#' ### Asymptotic Wilcoxon confidence interval
#'
#' The asymptotic Wilcoxon confidence interval is obtained by inverting the asymptotic normal approximation of the signed-rank statistic.
#' In detail, Define a standardized statistic:
#' \deqn{Z(d) = \frac{W^+(d) - E_0(d) - cc}{\sqrt{\mathrm{Var}_0(d)}}}
#' where \eqn{W^+(d)} is the signed-rank sum at shift \eqn{d}.
#' Then solve for \eqn{d} such that \eqn{Z(d)} equals the normal quantiles using [stats::uniroot()].
#' This interval may not be reliable for small samples or heavily tied data.
#'
#' ### Bootstrap confidence intervals
#'
#' The percentile and BCa bootstrap confidence interval methods are described in chapter 5.3 of \insertCite{davison1997;textual}{rankdifferencetest}.
#'
#' ## History
#'
#' [stats::wilcox.test()] is the canonical function for the Wilcoxon signed-rank test.
#' [exactRankTests::wilcox.exact()] implemented the Streitberg-Rohmel Shift-algorithm for exact inference when zeros and/or ties are present.
#' [coin::coin] superseded \pkg{exactRankTests} and includes additional methods.
#' [srt()] reimplements these functions so the best features of each is available in a fast and easy to use format.
#'
#' @references
#' \insertRef{wilcoxon1945}{rankdifferencetest}
#'
#' \insertRef{wilcoxon1950}{rankdifferencetest}
#'
#' \insertRef{pratt1981}{rankdifferencetest}
#'
#' \insertRef{pratt1959}{rankdifferencetest}
#'
#' \insertRef{cureton1967}{rankdifferencetest}
#'
#' \insertRef{conover1973}{rankdifferencetest}
#'
#' \insertRef{hollander2014}{rankdifferencetest}
#'
#' \insertRef{bauer1972}{rankdifferencetest}
#'
#' \insertRef{streitberg1984}{rankdifferencetest}
#'
#' \insertRef{hothorn2001}{rankdifferencetest}
#'
#' \insertRef{hothorn2002}{rankdifferencetest}
#'
#' \insertRef{hothorn2008}{rankdifferencetest}
#'
#' \insertRef{davison1997}{rankdifferencetest}
#'
#' \insertRef{exactRankTests}{rankdifferencetest}
#'
#' \insertRef{R}{rankdifferencetest}
#'
#' @param x
#' (numeric)\cr
#' Numeric vector of data.
#' Values with non-finite values (infinite or missing) are silently dropped.
#'
#' @param y
#' (numeric: `NULL`)\cr
#' Numeric vector of data or `NULL`.
#' If `NULL` (default), a one-sample test is performed using `x`.
#' If numeric, differences are calculated as `x - y`.
#' Pairs with non-finite values (infinite or missing) are silently dropped.
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
#'     Pairs with non-finite values (infinite or missing) are silently dropped.
#'     See `agg_fun` for handling duplicate cases of grouping/blocking combinations.
#'   }
#'   \item{y ~ x}{
#'     Use when `data` is in wide format.
#'     `y` and `x` must be numeric vectors.
#'     Differences are calculated as `data$y - data$x`.
#'     Pairs with non-finite values (infinite or missing) are silently dropped.
#'   }
#'   \item{ ~ x}{
#'     Use when `data$x` represents pre-calculated differences or for the one-sample case.
#'     Values with non-finite values (infinite or missing) are silently dropped.
#'   }
#' }
#'
#' @param mu
#' (Scalar numeric: `0`; `(-Inf, Inf)`)\cr
#' Under the null hypothesis, `x` or `x - y` is assumed to be symmetric around `mu`.
#'
#' @param zero_method
#' (Scalar character: `c("wilcoxon", "pratt")`)\cr
#' The method used to handle differences equal to zero.
#' If `"wilcoxon"` (default), zeros are removed prior to ranking (classic Wilcoxon convention).
#' If `"pratt"`, zeros are retained for ranking but excluded from the signed-rank sum.
#'
#' @inheritParams rdt
#'
#' @seealso
#' [rankdifferencetest::rdt()],
#' [stats::wilcox.test()],
#' [coin::wilcoxsign_test()]
#'
#' @returns
#' A list with the following elements:
#' \tabular{llll}{
#'   Slot \tab Subslot \tab Name \tab Description \cr
#'   1 \tab \tab `p_value`            \tab p-value. \cr
#'   2 \tab \tab `statistic`    \tab Test statistic. \eqn{W^+} for the exact Wilcoxon signed-rank distribution. \eqn{Z} for the asymptotic normal approximation. \cr
#'   3 \tab \tab `pseudomedian` \tab Measure of centrality. \cr
#'   4 \tab \tab `lower`        \tab Lower bound of confidence interval for the pseudomedian. `NULL` if no CI requested. \cr
#'   5 \tab \tab `upper`        \tab Upper bound of confidence interval for the pseudomedian. `NULL` if no CI requested. \cr
#'   6 \tab \tab `method`       \tab Test method. \cr
#'   7 \tab   \tab `info` \tab Additional test information. \cr
#'   7 \tab 1 \tab `p_value_method` \tab Method used to calculate the p-value. \cr
#'   7 \tab 2 \tab `pseudomedian_method` \tab Method used to calculate the pseudomedian. \cr
#'   7 \tab 3 \tab `conf_method` \tab Method used to calculate the confidence interval. \cr
#'   7 \tab 4 \tab `conf_level_achieved` \tab Achieved confidence level. \cr
#'   7 \tab 5 \tab `n_sample` \tab Number of observations in the original data. \cr
#'   7 \tab 6 \tab `n_analytic` \tab Number of observations after removing non-finite values from the original data. \cr
#'   7 \tab 7 \tab `n_zeros` \tab Number of zeros among differences in the analytic data set. \cr
#'   7 \tab 8 \tab `n_signed` \tab Number of nonzero differences in the analytic data set. \cr
#'   7 \tab 9 \tab `n_ties` \tab Number of tied ranks after ranking the absolute differences. \cr
#'   7 \tab 10 \tab `data_type` \tab Data type. \cr
#'   7 \tab 11 \tab `focal_name` \tab Name of the focal variable (differences are focal - reference). \cr
#'   7 \tab 12 \tab `reference_name` \tab Name of the reference variable (differences are focal - reference). \cr
#'   8 \tab   \tab `call` \tab A named list of the function's arguments (use `as.call()` to convert to a call; `call$distribution` may be updated from `"exact"` to `"asymptotic"`).
#' }
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # srt() example
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
#' # Rate transformation inverts the rank ordering.
#' data$placebo_rate <- 60 / data$placebo
#' data$drug_rate <- 60 / data$drug
#' data_tall$rate <- 60 / data_tall$time
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
#' @rdname srt
#' @export
srt <- function(
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
    quote(srt),
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
  res <- srt_data(x = data, formula = formula, agg_fun = agg_fun) |>
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

#' @rdname srt
#' @export
srt2 <- function(
  x,
  y = NULL,
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
  if (!is.null(y)) {
    if (!is.numeric(y)) {
      stop("Argument 'y' must be an object of class 'numeric'.")
    }
    if (length(x) != length(y)) {
      stop("Arguments 'x' and 'y' must be the same length.")
    }
    reference_name <- deparse1(substitute(y))
  } else {
    reference_name <- NULL
  }
  focal_name <- deparse1(substitute(x))

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
    quote(srt2),
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
  res <- srt_data(x = x, y = y, x_name = focal_name, y_name = reference_name) |>
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
