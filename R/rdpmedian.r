#' Rank difference pseudomedian
#'
#' Computes the Hodges-Lehmann pseudomedian and bootstrap confidence interval for the differences of ranks.
#'
#' This function generates a confidence interval for the pseudomedian based on the observed differences of ranks, not based on an inversion of the rank difference test [rdt()].
#'
#' The Hodges-Lehmann estimator is the median of all pairwise averages of the sample values.
#' \deqn{\mathrm{HL} = \mathrm{median} \left\{ \frac{x_i + x_j}{2} \right\}_{i \le j}}
#' This pseudomedian is a robust, distribution-free estimate of central tendency for a single sample, or a location-shift estimator for paired data.
#' It's resistant to outliers and compatible with rank-based inference.
#'
#' The percentile and BCa bootstrap confidence interval methods are described in chapter 5.3 of \insertCite{davison1997;textual}{rankdifferencetest}.
#'
#' This function is mainly a wrapper for the function [Hmisc::pMedian()].
#'
#' @references
#' \insertRef{davison1997}{rankdifferencetest}
#'
#' \insertRef{Hmisc}{rankdifferencetest}
#'
#' @inheritParams rdt
#'
#' @param conf_level
#' (Scalar numeric: `0.95`; `[0, 1)`)\cr
#' The confidence level.
#' If `conf_level = 0`, no confidence interval is calculated.
#'
#' @param conf_method
#' (Scalar character: `c("percentile", "bca")`)\cr
#' The type of bootstrap confidence interval.
#'
#' @param n_resamples
#' (Scalar integer: `[10L, Inf)`\cr
#' The number of bootstrap resamples.
#' If `conf_level = 0`, no resampling is performed.
#'
#' @seealso
#' [pmedian()]
#'
#' @return
#' A list with the following elements:
#' \tabular{llll}{
#'   Slot \tab Subslot \tab Name \tab Description \cr
#'   1 \tab \tab `pseudomedian` \tab Measure of centrality. \cr
#'   2 \tab \tab `lower`        \tab Lower bound of confidence interval for the pseudomedian. \cr
#'   3 \tab \tab `upper`        \tab Upper bound of confidence interval for the pseudomedian. \cr
#'   4 \tab \tab `method`       \tab Estimate method. \cr
#'   5 \tab   \tab `info` \tab Additional information. \cr
#'   5 \tab 1 \tab `n_sample` \tab Number of observations in the original data. \cr
#'   5 \tab 2 \tab `n_analytic` \tab Number of observations after removing missing values from the original data. \cr
#'   5 \tab 3 \tab `data_type` \tab Data type. \cr
#'   5 \tab 4 \tab `focal_name` \tab Name of the focal variable (differences are focal - reference). \cr
#'   5 \tab 5 \tab `reference_name` \tab Name of the reference variable (differences are focal - reference). \cr
#'   6 \tab   \tab `call` \tab A named list of the function's arguments (use `as.call()` to convert to a call).
#' }
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # rdpmedian() example
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
#' # Estimates
#' rdpmedian(
#'   data = data,
#'   formula = placebo ~ drug
#' )
#'
#' rdpmedian(
#'   data = data_tall,
#'   formula = time ~ treatment | subject
#' )
#'
#' rdpmedian2(
#'   x = data$placebo_rate,
#'   y = data$drug_rate
#' )
#'
#' rdpmedian(
#'   data = data_tall,
#'   formula = rate ~ treatment | subject
#' )
#'
#' @rdname rdpmedian
#' @export
rdpmedian <- function(
  data,
  formula,
  conf_level = 0.95,
  conf_method = "percentile",
  n_resamples = 1000L,
  agg_fun = "error"
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
    !(!is.na(conf_method) &&
      is.character(conf_method) &&
      length(conf_method) == 1L &&
      conf_method %in% c("percentile", "bca"))
  ) {
    msg <- "Argument 'conf_method' must be one of 'percentile' or 'bca'."
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

  # Could use match.call(), but will stay consistent with srt() for now.
  call <- list(
    quote(rdpmedian),
    data = substitute(data),
    formula = formula,
    conf_level = conf_level,
    conf_method = conf_method,
    n_resamples = n_resamples,
    agg_fun = agg_fun
  )

  #-----------------------------------------------------------------------------
  # Prepare results
  #-----------------------------------------------------------------------------
  lst <- rdt_data(x = data, formula = formula, agg_fun = agg_fun)

  res <- Hmisc::pMedian(
    x = lst$diffs,
    na.rm = TRUE,
    conf.int = conf_level,
    B = n_resamples,
    type = conf_method
  )

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  data_type <- "Paired differences of ranks"

  res <- if (conf_level == 0) {
    list(
      pseudomedian = as.numeric(res),
      lower = NULL,
      upper = NULL,
      method = "Hodges-Lehmann estimator"
    )
  } else {
    method <- switch(
      conf_method,
      "percentile" = "Hodges-Lehmann estimator and percentile bootstrap confidence interval",
      "bca" = "Hodges-Lehmann estimator and BCa bootstrap confidence interval",
      stop("Unknown bootstrap CI")
    )
    list(
      pseudomedian = res[["estimate"]],
      lower = res[["lower"]],
      upper = res[["upper"]],
      method = method
    )
  }

  info <- list(
    n_sample = lst$n_sample,
    n_analytic = lst$n_analytic,
    data_type = data_type,
    focal_name = lst$focal_name,
    reference_name = lst$reference_name
  )

  res <- c(res, list(info = info), list(call = call))
  class(res) <- c("rdpmedian", "list")

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}

#' @rdname rdpmedian
#' @export
rdpmedian2 <- function(
  x,
  y,
  conf_level = 0.95,
  conf_method = "percentile",
  n_resamples = 1000L
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
    !(!is.na(conf_method) &&
      is.character(conf_method) &&
      length(conf_method) == 1L &&
      conf_method %in% c("percentile", "bca"))
  ) {
    msg <- "Argument 'conf_method' must be one of 'percentile' or 'bca'."
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

  # Could use match.call(), but will stay consistent with srt() for now.
  call <- list(
    quote(rdpmedian2),
    x = substitute(x),
    y = substitute(y),
    conf_level = conf_level,
    conf_method = conf_method,
    n_resamples = n_resamples
  )

  #-----------------------------------------------------------------------------
  # Prepare results
  #-----------------------------------------------------------------------------
  lst <- rdt_data(x = x, y = y, x_name = focal_name, y_name = reference_name)

  res <- Hmisc::pMedian(
    x = lst$diffs,
    na.rm = TRUE,
    conf.int = conf_level,
    B = n_resamples,
    type = conf_method
  )

  #-----------------------------------------------------------------------------
  # Prepare return
  #-----------------------------------------------------------------------------
  data_type <- "Paired differences of ranks"

  res <- if (conf_level == 0) {
    list(
      pseudomedian = as.numeric(res),
      lower = NULL,
      upper = NULL,
      method = "Hodges-Lehmann estimator"
    )
  } else {
    method <- switch(
      conf_method,
      "percentile" = "Hodges-Lehmann estimator and percentile bootstrap confidence interval",
      "bca" = "Hodges-Lehmann estimator and BCa bootstrap confidence interval",
      stop("Unknown bootstrap CI")
    )
    list(
      pseudomedian = res[["estimate"]],
      lower = res[["lower"]],
      upper = res[["upper"]],
      method = method
    )
  }

  info <- list(
    n_sample = lst$n_sample,
    n_analytic = lst$n_analytic,
    data_type = data_type,
    focal_name = lst$focal_name,
    reference_name = lst$reference_name
  )

  res <- c(res, list(info = info), list(call = call))
  class(res) <- c("rdpmedian", "list")

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  res
}
