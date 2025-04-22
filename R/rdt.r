#' Rank difference test
#'
#' Performs Kornbrot's rank difference test. It's a modified Wilcoxon signed-rank
#' test which produces consistent and meaningful results for ordinal or monotonically
#' transformed data.
#'
#' For ordinal scale data, the Wilcoxon signed-rank test results in subtraction
#' of those ordinal scale values. However, this subtraction is not meaningful on
#' the ordinal scale. In addition, any monotone transformation of the data will
#' result in different signed ranks, thus different p-values. However, ranking
#' the original data allows for meaningful addition and subtraction of ranks and
#' preserves ranks over monotonic transformation. Kornbrot described the rank
#' difference test for this reason.
#'
#' Kornbrot recommends that the rank difference procedure be used in preference
#' to the Wilcoxon signed-rank test in all paired comparison designs where the data
#' are not both of interval scale type and of known distribution type. The rank
#' difference test preserves good power compared to Wilcoxon's signed-rank test,
#' is more powerful than the sign test, and has the benefit of being a true
#' distribution-free test.
#'
#' The procedure for Wilcoxon's signed-rank test is as follows:
#'
#' 1. Calculate differences for each paired observation.
#' 2. Remove differences equal to zero.
#' 3. Order the absolute differences from smallest to largest.
#' 4. Assign ranks \eqn{1, \dots, n} with average rank for ties.
#' 5. Calculate W+ = sum of the ranks for positive differences. The sum of W+ and
#'    W- is \eqn{n(n+1)/2}, so either can be calculated from the other.
#' 6. Choose the smaller of W+ and W- as the test statistic W.
#' 7. Since the test statistic is the smaller of W+ and W-, the critical region
#'    is the left tail of the distribution. W is distributed approximately normal
#'    with mean \eqn{mu = (n(n+1))/4} and variance
#'    \eqn{sigma^2 = (Tn(n+1)(2n+1))/24}, where T is a correction for ties and
#'    \eqn{T = 1-(sum(t^3-t)/(N^3-N))}, summed over all ties, where t is the
#'    length of a tie. The continuity corrected mean \eqn{mu = ((n(n+1))/4)+0.5}.
#'
#' The procedure for Kornbrot's rank difference test is as follows:
#'
#' 1. Combine all 2n observations.
#' 2. Assign ranks \eqn{1, \dots, 2n} with average rank for ties.
#' 3. Perform the Wilcoxon signed-rank test using the paired ranks.
#'
#' The test statistic for the rank difference test (D) is not exactly equal to the
#' test statistic (W) of the naive rank-transformed Wilcoxon signed-rank test
#' (the latter being implemented in `rdt()`). Using W should result in a
#' conservative estimate for D, and they approach in distribution as the sample
#' size increases. \insertCite{kornbrot1990;textual}{rankdifferencetest}
#' discusses methods for calculating D when n<7 and 8<n<=20.
#'
#' `zero.method = "Pratt"` uses the method by Pratt (1959), which first
#' rank-transforms the absolute differences (including zeros) and then removes
#' the ranks corresponding to zero-differences. `zero.method = "Wilcoxon"` uses
#' the method by Wilcoxon (1950), which first removes the zero-differences and
#' then rank-transforms the remaining absolute differences.
#'
#' @param data A data.frame
#' @param formula  A formula of either form:
#'        \describe{
#'          \item{y ~ x | block}{
#'            For use when `data` is in long format. `y` is the numeric outcome,
#'            `x` is the binary grouping factor, and `block` is the
#'            observation-level grouping factor. `x` and `block` must be factors.
#'          }
#'          \item{y ~ x}{
#'            For use when `data` is in wide format. `y` is the first measurement
#'            and `x` is the second measurement on the same observation.
#'          }
#'        }
#'        The differences in ranks are calculated as `y - x` or, when `x` is a
#'        factor, the first factor level minus the second.
#' @param zero.method A string for the method used to handle differences equal to
#'        zero: `"Pratt"` (default) or `"Wilcoxon"`.
#' @param distribution A string for the method used to calculate the conditional
#'        null distribution of the test statistic: Asymptotic distribution
#'        `"asymptotic"` (default), Monte Carlo resampling `"approximate"`,
#'        or the exact distribution `"exact"`.
#' @param alternative A string for the alternative hypothesis: `"two.sided"`
#'        (default), `"greater"`, or `"less"`.
#' @param return A string for the return object: `"data.frame"` (default) or
#'        `"coin"`.
#' @param ... Further arguments passed to \code{coin::\link[coin]{wilcoxsign_test}}.
#'
#' @references
#' \insertRef{kornbrot1990}{rankdifferencetest}
#'
#' \insertRef{pratt1959}{rankdifferencetest}
#'
#' \insertRef{wilcoxon1950}{rankdifferencetest}
#'
#' @return
#' If `return = "data.frame"`, a data.frame with columns `p.value`, `z.statistic`,
#' `formula`, `alternative`, and `method`.
#'
#' If `return = "coin"`, an object inheriting from class \link[coin]{IndependenceTest}.
#'
#' @importFrom coin wilcoxsign_test pvalue
#' @importFrom stats as.formula
#'
#' @export
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
#' # Create long-format data for demonstration purposes
#' data_long <- reshape(
#'   data = kornbrot_table1,
#'   direction = "long",
#'   varying = c("placebo", "drug"),
#'   v.names = c("time"),
#'   idvar = "subject",
#'   times = c("placebo", "drug"),
#'   timevar = "treatment",
#'   new.row.names = seq_len(prod(length(c("placebo", "drug")), nrow(kornbrot_table1)))
#' )
#' # Subject and treatment must be factors. The ordering of the treatment factor
#' # will determine the difference (placebo - drug).
#' data_long$subject <- factor(data_long$subject)
#' data_long$treatment <- factor(data_long$treatment, levels = c("placebo", "drug"))
#'
#' # Recreate analysis and results from section 7.1 in Kornbrot (1990)
#' ## The p-value shown in Kornbrot (1990) was continuity corrected. rdt() does
#' ## not apply a continuity correction, so the p-value here will be slightly
#' ## lower. It does match the uncorrected p-value shown in footnote on page 246.
#' rdt(
#'   data = data,
#'   formula = placebo ~ drug,
#'   alternative = "two.sided",
#'   distribution = "asymptotic"
#' )$p.value/2
#' rdt(
#'   data = data_long,
#'   formula = time ~ treatment | subject,
#'   alternative = "two.sided",
#'   distribution = "asymptotic"
#' )$p.value/2
#'
#' # The same outcome is seen after transforming time to rate.
#' ## The rate transformation inverts the rank ordering.
#' data$placebo_rate <- 60 / data$placebo
#' data$drug_rate <- 60 / data$drug
#' data_long$rate <- 60 / data_long$time
#'
#' rdt(
#'   data = data,
#'   formula = placebo_rate ~ drug_rate,
#'   alternative = "two.sided",
#'   distribution = "asymptotic"
#' )$p.value/2
#' rdt(
#'   data = data_long,
#'   formula = rate ~ treatment | subject,
#'   alternative = "two.sided",
#'   distribution = "asymptotic"
#' )$p.value/2
#'
#' # In contrast to the rank difference test, the Wilcoxon signed-rank test
#' # produces differing results. See table 1 and table 2 in Kornbrot (1990).
#' wilcox.test(
#'   x = data$placebo,
#'   y = data$drug,
#'   paired = TRUE,
#'   exact = TRUE,
#'   alternative = "two.sided"
#' )$p.value/2
#' wilcox.test(
#'   x = data$placebo_rate,
#'   y = data$drug_rate,
#'   paired = TRUE,
#'   exact = TRUE,
#'   alternative = "two.sided"
#' )$p.value/2
#'
rdt <- function(
  data,
  formula,
  zero.method = c("Pratt", "Wilcoxon"),
  distribution = c("asymptotic", "approximate", "exact"),
  alternative = c("two.sided", "greater", "less"),
  return = c("data.frame", "coin"),
  ...
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if(!is.data.frame(data)) {
    stop("Check the 'data' argument in rdt(). The object must be a data.frame.")
  }
  if(!inherits(formula, "formula")) {
    stop("Check the 'formula' argument in rdt(). The object must be a formula.")
  }
  formula_vars <- all.vars(formula)
  missing_vars <- formula_vars[!formula_vars %in% names(data)]
  if(length(missing_vars) > 0) {
    stop(paste0(
      "Check the formula argument in rdt(). Variables not found in the data: ",
      paste0(shQuote(missing_vars), collapse = ", ")
    ))
  }
  if(!length(formula_vars) %in% 2:3) {
    stop("Check the formula argument in rdt(). It must be of form y ~ x or y ~ x | block.")
  }
  zero.method <- match.arg(zero.method)
  distribution <- match.arg(distribution)
  alternative <- match.arg(alternative)
  return <- match.arg(return)

  #-----------------------------------------------------------------------------
  # Perform test
  #-----------------------------------------------------------------------------
  ## Step 1: Generate model data
  md <- model2data(formula = formula, data = data)

  ## Step 2: Rank order data and prepare for coin::wilcoxsign_test()
  if(is.null(md$block)) {
    # Intermediate argument check: if no block, must only have two variables.
    if(!all(lengths(md) == c(1, 1, 0))) {
      stop("Check the formula used in rdt(). It must correspond with a pairwise comparison as y ~ x.")
    }
    # Note: Kornbrot orders opposite direction: rank(-xtfrm(data)).
    ranks <- rank(c(md$y[[1]], md$x[[1]]), na.last = "keep")
    coin_data <- data.frame(
      y = ranks[seq_along(md$y[[1]])],
      x = ranks[seq_along(md$x[[1]]) + length(md$y[[1]])]
    )
    coin_formula <- as.formula("y ~ x")
  } else {
    # Intermediate argument check: must be a balanced pairwise comparison.
    if(length(unique(md$x[[1]])) != 2) {
      stop(paste0(
        "Check the data and/or formula used in rdt(). The variable ",
        shQuote(names(md$x)),
        " must be a balanced 2 level factor."
      ))
    }
    coin_data <- data.frame(
      y = rank(md$y[[1]], na.last = "keep"),
      x = as.factor(md$x[[1]]),
      block = as.factor(md$block[[1]])
    )
    coin_formula <- as.formula("y ~ x | block")
  }

  ## Step 3: Run Wilcoxon signed-rank test
  coin <- wilcoxsign_test(
    formula = coin_formula,
    data = coin_data,
    zero.method = zero.method,
    distribution = distribution,
    alternative = alternative,
    ...
  )

  #-----------------------------------------------------------------------------
  # Prepare results
  #-----------------------------------------------------------------------------
  diff <- if(is.null(md$block)) {
    paste0(names(md$y), " - ", names(md$x))
  } else {
    paste0(levels(coin_data$x)[1], " - ", levels(coin_data$x)[2])
  }

  alternative <- paste0(
    "True location shift of ranks (",
    diff,
    ") is ",
    ifelse(
      test = identical(coin@statistic@alternative, "two.sided"),
      yes = "not equal to ",
      no = paste0(coin@statistic@alternative, " than ")
    ),
    coin@nullvalue
  )

  method <- paste(
    "Kornbrot's Rank Difference Test using the",
    paste0(
      toupper(substring(distribution, 1,1)),
      substring(distribution, 2)
    ),
    coin@method,
    sep = " "
  )

  df <- data.frame(
    p.value = as.numeric(pvalue(coin)),
    z.statistic = coin@statistic@teststatistic,
    formula = deparse(formula),
    alternative = alternative,
    method = method
  )

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  if(identical(x = return, y = "data.frame")) {
    df
  } else {
    coin
  }
}
