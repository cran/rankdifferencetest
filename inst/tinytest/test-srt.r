library(rankdifferencetest)
library(tinytest)

# Data (paired differences are zero-free and tie-free)
data <- kornbrot_table1
data_tall <- reshape(
  data = kornbrot_table1,
  direction = "long",
  varying = c("placebo", "drug"),
  v.names = c("time"),
  idvar = "subject",
  times = c("placebo", "drug"),
  timevar = "treatment",
  new.row.names = seq_len(prod(
    length(c("placebo", "drug")),
    nrow(kornbrot_table1)
  ))
)
data_tall$subject <- factor(data_tall$subject)
data_tall$treatment <- factor(
  data_tall$treatment,
  levels = c("drug", "placebo")
)
data_tall$xtra <- seq_len(nrow(data_tall))

#-------------------------------------------------------------------------------
# Successful run
#-------------------------------------------------------------------------------
lst_names <- c(
  "p_value",
  "statistic",
  "pseudomedian",
  "lower",
  "upper",
  "method",
  "info",
  "call"
)
info_names <- c(
  "p_value_method",
  "pseudomedian_method",
  "conf_method",
  "conf_level_achieved",
  "n_sample",
  "n_analytic",
  "n_zeros",
  "n_signed",
  "n_ties",
  "data_type",
  "focal_name",
  "reference_name"
)
call_names <- c(
  "",
  "data",
  "formula",
  "conf_level",
  "conf_method",
  "n_resamples",
  "alternative",
  "mu",
  "distribution",
  "correct",
  "zero_method",
  "agg_fun",
  "digits_rank",
  "tol_root"
)
call2_names <- c(
  "",
  "x",
  "y",
  "conf_level",
  "conf_method",
  "n_resamples",
  "alternative",
  "mu",
  "distribution",
  "correct",
  "zero_method",
  "digits_rank",
  "tol_root"
)

# wide format
expect_silent(res <- srt(data = data, formula = placebo ~ drug))
expect_true(inherits(res, "srt"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call_names)

expect_silent(res <- srt2(x = data$placebo, y = data$drug))
expect_true(inherits(res, "srt"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call2_names)

# one sample
expect_silent(res <- srt2(x = data$placebo))
expect_true(inherits(res, "srt"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call2_names)

# tall format
expect_silent(
  res <- srt(data = data_tall, formula = time ~ treatment | subject)
)
expect_true(inherits(res, "srt"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call_names)

# with inversion CI
expect_silent(res <- srt2(x = data$placebo, y = data$drug, conf_level = 0.95))
expect_true(inherits(res, "srt"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call2_names)

# with bootstrap CI
expect_silent(
  res <- srt2(
    x = data$placebo,
    y = data$drug,
    conf_level = 0.95,
    conf_method = "percentile"
  )
)
expect_true(inherits(res, "srt"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call2_names)

#-------------------------------------------------------------------------------
# srt is same as srt2
#-------------------------------------------------------------------------------
# Exact test
res1 <- srt(data = data, formula = placebo ~ drug)
res2 <- srt(data = data_tall, formula = time ~ treatment | subject)
res3 <- srt2(x = data$placebo, y = data$drug)

res1$call$data <- res1$call$formula <- res1$call[[1]] <- NULL
res1$info$focal_name <- res1$info$reference_name <- NULL
res1$call$agg_fun <- NULL

res2$call$data <- res2$call$formula <- res2$call[[1]] <- NULL
res2$info$focal_name <- res2$info$reference_name <- NULL
res2$call$agg_fun <- NULL

res3$call$x <- res3$call$y <- res3$call[[1]] <- NULL
res3$info$focal_name <- res3$info$reference_name <- NULL

expect_identical(res1, res2)
expect_identical(res1, res3)

# Asymptotic test
res1 <- srt(data = data, formula = placebo ~ drug, distribution = "asymptotic")
res2 <- srt(
  data = data_tall,
  formula = time ~ treatment | subject,
  distribution = "asymptotic"
)
res3 <- srt2(x = data$placebo, y = data$drug, distribution = "asymptotic")

res1$call$data <- res1$call$formula <- res1$call[[1]] <- NULL
res1$info$focal_name <- res1$info$reference_name <- NULL
res1$call$agg_fun <- NULL

res2$call$data <- res2$call$formula <- res2$call[[1]] <- NULL
res2$info$focal_name <- res2$info$reference_name <- NULL
res2$call$agg_fun <- NULL

res3$call$x <- res3$call$y <- res3$call[[1]] <- NULL
res3$info$focal_name <- res3$info$reference_name <- NULL

expect_identical(res1, res2)
expect_identical(res1, res3)

# Exact test with CI
res1 <- srt(data = data, formula = placebo ~ drug, conf_level = 0.95)
res2 <- srt(
  data = data_tall,
  formula = time ~ treatment | subject,
  conf_level = 0.95
)
res3 <- srt2(x = data$placebo, y = data$drug, conf_level = 0.95)

res1$call$data <- res1$call$formula <- res1$call[[1]] <- NULL
res1$info$focal_name <- res1$info$reference_name <- NULL
res1$call$agg_fun <- NULL

res2$call$data <- res2$call$formula <- res2$call[[1]] <- NULL
res2$info$focal_name <- res2$info$reference_name <- NULL
res2$call$agg_fun <- NULL

res3$call$x <- res3$call$y <- res3$call[[1]] <- NULL
res3$info$focal_name <- res3$info$reference_name <- NULL

expect_identical(res1, res2)
expect_identical(res1, res3)

# Asymptotic test with CI
res1 <- srt(
  data = data,
  formula = placebo ~ drug,
  distribution = "asymptotic",
  conf_level = 0.95
)
res2 <- srt(
  data = data_tall,
  formula = time ~ treatment | subject,
  distribution = "asymptotic",
  conf_level = 0.95
)
res3 <- srt2(
  x = data$placebo,
  y = data$drug,
  distribution = "asymptotic",
  conf_level = 0.95
)

res1$call$data <- res1$call$formula <- res1$call[[1]] <- NULL
res1$info$focal_name <- res1$info$reference_name <- NULL
res1$call$agg_fun <- NULL

res2$call$data <- res2$call$formula <- res2$call[[1]] <- NULL
res2$info$focal_name <- res2$info$reference_name <- NULL
res2$call$agg_fun <- NULL

res3$call$x <- res3$call$y <- res3$call[[1]] <- NULL
res3$info$focal_name <- res3$info$reference_name <- NULL

expect_identical(res1, res2)
expect_identical(res1, res3)

#-------------------------------------------------------------------------------
# Bootstrap CI is successful
#-------------------------------------------------------------------------------
res <- srt(
  data = data,
  formula = placebo ~ drug,
  conf_level = 0.95,
  conf_method = "percentile"
)
expect_true(res$lower < 0 && res$lower > -5)
expect_true(res$upper > 0 && res$upper < 4)
expect_true(res$info$conf_method == "Percentile bootstrap")

#-------------------------------------------------------------------------------
# Argument checks return errors
#-------------------------------------------------------------------------------
expect_error(
  srt(data = 1, formula = placebo ~ drug),
  pattern = "Argument 'data' must be an object of class 'data.frame'.",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = data),
  pattern = "Argument 'formula' must be an object of class 'formula'.",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = placebo ~ drug, conf_level = "0.95"),
  pattern = "Argument 'conf_level' must be a number between 0 and 1.",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = placebo ~ drug, conf_level = 1.1),
  pattern = "Argument 'conf_level' must be a number between 0 and 1.",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = placebo ~ drug, alternative = "two"),
  pattern = "Argument 'alternative' must be one of 'two.sided', 'greater', or 'less'.",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = placebo ~ drug, mu = "0"),
  pattern = "Argument 'mu' must be a scalar vector of class 'numeric'.",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = placebo ~ drug, distribution = "Exact"),
  pattern = "Argument 'distribution' must be one of 'auto', 'exact', or 'asymptotic'.",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = placebo ~ drug, correct = 1),
  pattern = "Argument 'correct' must be a scalar vector of class 'logical'.",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = placebo ~ drug, zero_method = "Wilcoxon"),
  pattern = "Argument 'zero_method' must be one of 'wilcoxon' or 'pratt'.",
  fixed = TRUE
)
expect_error(
  srt(
    data = data_tall,
    formula = time ~ treatment | subject,
    agg_fun = "invalid"
  ),
  pattern = "Argument 'agg_fun' must be one of 'error', 'first', 'last', 'sum', 'mean', 'median', 'min', or 'max'.",
  fixed = TRUE
)
expect_error(
  srt(data = data_tall, formula = time ~ treatment | subject, digits_rank = 0),
  pattern = "Argument 'digits_rank' must be a scalar integer > 0 (may be Inf).",
  fixed = TRUE
)
expect_error(
  srt(data = data_tall, formula = time ~ treatment | subject, tol_root = 0),
  pattern = "Argument 'tol_root' must be a scalar numeric > 0.",
  fixed = TRUE
)

#-------------------------------------------------------------------------------
# Malformed formula returns errors
#-------------------------------------------------------------------------------
expect_error(
  srt(data = data, formula = placebo ~ drug2),
  pattern = "The formula included terms not found in the data: 'drug2'",
  fixed = TRUE
)
expect_error(
  srt(data = data, formula = placebo ~ drug + subject),
  pattern = "The formula must not have multiple 'group' components. For example, it must be of form 'y ~ x | z', 'y ~ x', or '~x'.",
  fixed = TRUE
)
expect_error(
  srt(data = data_tall, formula = time ~ treatment + xtra | subject),
  pattern = "The formula must not have multiple 'group' components. For example, it must be of form 'y ~ x | z', 'y ~ x', or '~x'.",
  fixed = TRUE
)

expect_warning(
  srt(data = data_tall, formula = time ~ treatment | xtra),
  pattern = "'xtra' was converted to a factor in formula 'time ~ treatment | xtra'.",
  fixed = TRUE
)
data_tall$xtra <- factor(data_tall$xtra)
expect_warning(
  srt(data = data_tall, formula = time ~ treatment | xtra),
  pattern = "The number of valid observations (non-missing and finite) must be at least 3.",
  fixed = TRUE
)
expect_error(
  srt(data = data_tall, formula = time ~ treatment | treatment),
  pattern = "Duplicate observations found for (treatment, treatment) cells",
  fixed = TRUE
)

data_tall$trt3 <- factor(rep(1:3, length.out = nrow(data_tall)))
expect_error(
  srt(data = data_tall, formula = time ~ trt3 | subject),
  pattern = "'trt3' must be a two-level factor in formula 'time ~ trt3 | subject'.",
  fixed = TRUE
)

data_tall$trt1 <- factor(rep(1, length.out = nrow(data_tall)))
expect_error(
  srt(data = data_tall, formula = time ~ trt1 | subject),
  pattern = "'trt1' must be a two-level factor in formula 'time ~ trt1 | subject'.",
  fixed = TRUE
)

#-------------------------------------------------------------------------------
# p-value matches Kornbrot 1990
#-------------------------------------------------------------------------------
# tall format expected result
expect_equivalent(
  srt(
    data = data_tall,
    formula = time ~ treatment | subject,
    alternative = "greater",
    distribution = "asymptotic",
    zero_method = "wilcoxon",
    correct = FALSE
  )$p_value,
  0.5693479,
  tolerance = 1e-6,
)

# wide format expected result
expect_equivalent(
  srt(
    data = data,
    formula = placebo ~ drug,
    alternative = "greater",
    distribution = "asymptotic",
    zero_method = "wilcoxon",
    correct = FALSE
  )$p_value,
  0.5693479,
  tolerance = 1e-6
)

# transformed wide format
expect_equivalent(
  srt(
    data = data,
    formula = I(60 / placebo) ~ I(60 / drug),
    alternative = "less",
    distribution = "asymptotic",
    zero_method = "wilcoxon",
    correct = FALSE
  )$p_value,
  0.04343004,
  tolerance = 1e-6
)

# transformed tall format
expect_equivalent(
  srt(
    data = data_tall,
    formula = I(60 / time) ~ treatment | subject,
    alternative = "less",
    distribution = "asymptotic",
    zero_method = "wilcoxon",
    correct = FALSE
  )$p_value,
  0.04343004,
  tolerance = 1e-6
)

#-------------------------------------------------------------------------------
# Same results after reordering
#-------------------------------------------------------------------------------
# Blocking handles random rows (same results after reordering)
data_tall2 <- data_tall[sample(seq_len(nrow(data_tall)), nrow(data_tall)), ]
res1 <- srt(
  data = data_tall,
  formula = time ~ treatment | subject,
  distribution = "asymptotic",
  zero_method = "wilcoxon",
  correct = FALSE
)
res2 <- srt(
  data = data_tall2,
  formula = time ~ treatment | subject,
  distribution = "asymptotic",
  zero_method = "wilcoxon",
  correct = FALSE
)
res1$call$data <- res2$call$data <- NULL
expect_identical(res1, res2)

# Same results after reordering
res1 <- srt(data = data, placebo ~ drug)
res2 <- srt(data = data[sample(nrow(data)), ], formula = placebo ~ drug)
res1$call$data <- NULL
res2$call$data <- NULL
expect_identical(res1, res2)

# Same results after reordering
tmp <- data
res1 <- srt2(tmp$placebo, tmp$drug)
tmp <- tmp[sample(nrow(tmp)), ]
res2 <- srt2(tmp$placebo, tmp$drug)
res1$call$data <- NULL
res2$call$data <- NULL
expect_identical(res1, res2)

#-------------------------------------------------------------------------------
# Return summaries are correct
# 1. Formula
#     - y ~ x
#     - y ~ x | block
# 2. Alternative
#     - two.sided
#     - greater
#     - less
# 3. Distribution
#     - exact
#     - asymptotic
#     - approximate
# 4. zero_method
#     - Wilcoxon
#     - Pratt
#-------------------------------------------------------------------------------
# Correct formula
res <- srt(
  data = data,
  formula = placebo ~ drug
)
res2 <- srt(
  data = data_tall,
  formula = time ~ treatment | subject
)
expect_equal(res$call$formula, placebo ~ drug)
expect_equal(res2$call$formula, time ~ treatment | subject)

# correct alternative
res <- srt(
  data = data,
  formula = placebo ~ drug,
  alternative = "two.sided"
)
res2 <- srt(
  data = data,
  formula = placebo ~ drug,
  alternative = "greater"
)
res3 <- srt(
  data = data,
  formula = placebo ~ drug,
  alternative = "less"
)
expect_equal(res$call$alternative, "two.sided")
expect_equal(res2$call$alternative, "greater")
expect_equal(res3$call$alternative, "less")

# correct method
res <- srt(
  data = data,
  formula = placebo ~ drug,
  distribution = "exact"
)
res2 <- srt(
  data = data,
  formula = placebo ~ drug,
  distribution = "asymptotic"
)
expect_equal(res$info$p_value_method, "Wilcoxon")
expect_equal(res2$info$p_value_method, "Asymptotic")

# correct zero-difference method
res <- srt(
  data = data,
  formula = placebo ~ drug,
  zero_method = "wilcoxon"
)
res2 <- srt(
  data = data,
  formula = placebo ~ drug,
  zero_method = "pratt"
)
expect_equal(
  res$call$zero_method,
  "wilcoxon"
)
expect_equal(
  res2$call$zero_method,
  "pratt"
)

#-------------------------------------------------------------------------------
# Confidence interval as expected
#-------------------------------------------------------------------------------
res <- srt(data = data, formula = placebo ~ drug)

# Confidence interval not requested
expect_null(res$lower)
expect_null(res$upper)
expect_null(res$info$conf_method)
expect_null(res$info$conf_level_achieved)
expect_equal(res$call$conf_level, 0)

# Confidence interval requested
res_ci <- srt(data = data, formula = placebo ~ drug, conf_level = 0.95)
expect_true(is.numeric(res_ci$lower) && !is.na(res_ci$lower))
expect_true(is.numeric(res_ci$upper) && !is.na(res_ci$upper))
expect_equal(
  res_ci$info$conf_method,
  "Inversion of exact Wilcoxon signed-rank test"
)
expect_equal(res_ci$info$conf_level_achieved, 0.951, tolerance = 0.01)
expect_equal(res_ci$call$conf_level, 0.95)

#-------------------------------------------------------------------------------
# Tall data, formula interface, and aggregation handles missing values
#-------------------------------------------------------------------------------
# orphaned group
data_tall2 <- data_tall
data_tall2 <- data_tall2[-1, ]
tmp <- srt(
  data = data_tall2,
  formula = time ~ treatment | subject
)
expect_equal(tmp$info$n_sample, 13)
expect_equal(tmp$info$n_analytic, 12)

# missing outcome
data_tall2 <- data_tall
data_tall2[1, "time"] <- NA
tmp <- srt(
  data = data_tall2,
  formula = time ~ treatment | subject
)
expect_equal(tmp$info$n_sample, 13)
expect_equal(tmp$info$n_analytic, 12)

# Missing group
data_tall2 <- data_tall
data_tall2[1, "treatment"] <- NA
tmp <- srt(
  data = data_tall2,
  formula = time ~ treatment | subject
)
expect_equal(tmp$info$n_sample, 13)
expect_equal(tmp$info$n_analytic, 12)

# Missing block
data_tall2 <- data_tall
data_tall2[1, "subject"] <- NA
tmp <- srt(
  data = data_tall2,
  formula = time ~ treatment | subject
)
expect_equal(tmp$info$n_sample, 13)
expect_equal(tmp$info$n_analytic, 12)

# Aggregation
data_tall2 <- data_tall
suppressWarnings(
  data_tall2[27, ] <- list(
    subject = 13,
    treatment = "drug",
    time = NA,
    xtra = 27
  )
)
suppressWarnings(
  data_tall2[28, ] <- list(
    subject = 12,
    treatment = "drug",
    time = NA,
    xtra = 28
  )
)
tmp <- srt(
  data = data_tall2,
  formula = time ~ treatment | subject,
  agg_fun = "first"
)
expect_equal(tmp$info$n_sample, 13)
expect_equal(tmp$info$n_analytic, 13)

# Error for unhandled duplicate group/block
dup <- data_tall
suppressWarnings(
  dup[27, ] <- list(subject = 13, treatment = "drug", time = NA, xtra = 27)
)
expect_error(
  srt(
    data = dup,
    formula = time ~ treatment | subject
  ),
  "Duplicate observations found for (subject, treatment)",
  fixed = TRUE
)

#-------------------------------------------------------------------------------
# mu shift as expected
# pseudomedian and CI bounds may not be the same after mu-shift...
#-------------------------------------------------------------------------------
set.seed(20251117)
data_wide <- data.frame(a = rnorm(10), b = rnorm(10))

res1 <- srt(data = data_wide, formula = a ~ b, mu = 0)
res2 <- srt(data = data_wide, formula = a ~ b, mu = 2)
expect_equal(res2$call$mu, 2)
expect_equal(res1$pseudomedian, res2$pseudomedian)

res1 <- srt(data = data_wide, formula = a ~ b, mu = 0, conf_level = 0.95)
res2 <- srt(data = data_wide, formula = a ~ b, mu = 2, conf_level = 0.95)
expect_equal(res2$call$mu, 2)
expect_equal(res1$pseudomedian, res2$pseudomedian)
expect_equal(res1$lower, res2$lower)
expect_equal(res1$upper, res2$upper)

res1 <- srt(
  data = data_wide,
  formula = a ~ b,
  mu = 0,
  distribution = "asymptotic"
)
res2 <- srt(
  data = data_wide,
  formula = a ~ b,
  mu = 2,
  distribution = "asymptotic"
)
expect_equal(res2$call$mu, 2)
expect_equal(res1$pseudomedian, res2$pseudomedian)

res1 <- srt(
  data = data_wide,
  formula = a ~ b,
  mu = 0,
  conf_level = 0.95,
  distribution = "asymptotic"
)
res2 <- srt(
  data = data_wide,
  formula = a ~ b,
  mu = 2,
  conf_level = 0.95,
  distribution = "asymptotic"
)
expect_equal(res2$call$mu, 2)
expect_equal(res1$pseudomedian, res2$pseudomedian)
expect_equal(res1$lower, res2$lower)
expect_equal(res1$upper, res2$upper)

#-------------------------------------------------------------------------------
# Invalid data
#-------------------------------------------------------------------------------
# All missing
data_na <- data.frame(a = NA_real_, b = NA_real_)
expect_warning(
  res <- srt(data = data_na, formula = a ~ b),
  "The number of valid observations (non-missing and finite) must be at least 3.",
  fixed = TRUE
)
expect_equal(res$pseudomedian, NA_real_)
expect_equal(res$p_value, NA_real_)
expect_equal(res$statistic, NA_real_, check.attributes = FALSE)

# All length 0
data_empty <- data.frame(a = numeric(0), b = numeric(0))
expect_warning(
  res <- srt(data = data_empty, formula = a ~ b),
  "The number of valid observations (non-missing and finite) must be at least 3.",
  fixed = TRUE
)
expect_equal(res$pseudomedian, NA_real_)
expect_equal(res$p_value, NA_real_)
expect_equal(res$statistic, NA_real_, check.attributes = FALSE)

## No zeros, all tied to same value (4)
x <- 5:8
y <- 1:4
expect_warning(
  res1 <- srt2(x, y, distribution = "exact", conf_level = 0.95),
  "Paired differences are constant/zero variance.",
  fixed = TRUE
)
expect_warning(
  res2 <- srt2(x, y, distribution = "asymptotic", conf_level = 0.95),
  "Paired differences are constant/zero variance.",
  fixed = TRUE
)

expect_equal(res1$pseudomedian, 4)
expect_equal(res1$p_value, NA_real_)
expect_equal(res1$statistic, NA_real_, check.attributes = FALSE)
expect_equal(res1$lower, NA_real_)
expect_equal(res1$upper, NA_real_)

expect_equal(res2$pseudomedian, 4)
expect_equal(res2$p_value, NA_real_)
expect_equal(res2$statistic, NA_real_, check.attributes = FALSE)
expect_equal(res2$lower, NA_real_)
expect_equal(res2$upper, NA_real_)

## All zeros
x <- c(1, 1, 1, 1)
y <- c(1, 1, 1, 1)
expect_warning(
  res1 <- srt2(x, y, distribution = "exact", conf_level = 0.95),
  "No paired differences left after removing zeros.",
  fixed = TRUE
)
expect_warning(
  res2 <- srt2(x, y, distribution = "asymptotic", conf_level = 0.95),
  "No paired differences left after removing zeros.",
  fixed = TRUE
)

expect_equal(res1$pseudomedian, NA_real_)
expect_equal(res1$p_value, NA_real_)
expect_equal(res1$statistic, NA_real_, check.attributes = FALSE)
expect_equal(res1$lower, NA_real_)
expect_equal(res1$upper, NA_real_)

expect_equal(res2$pseudomedian, NA_real_)
expect_equal(res2$p_value, NA_real_)
expect_equal(res2$statistic, NA_real_, check.attributes = FALSE)
expect_equal(res2$lower, NA_real_)
expect_equal(res2$upper, NA_real_)
