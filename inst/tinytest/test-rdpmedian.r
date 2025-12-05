library(rankdifferencetest)
library(tinytest)

set.seed(123)
data_wide <- data.frame(a = rnorm(10), b = rnorm(10))
data_tall <- data.frame(
  y = rnorm(20),
  x = factor(rep(c("grp1", "grp2"), each = 10)),
  z = factor(rep(1:10, times = 2))
)

#-------------------------------------------------------------------------------
# Successful run
#-------------------------------------------------------------------------------
lst_names <- c(
  "pseudomedian",
  "lower",
  "upper",
  "method",
  "info",
  "call"
)
info_names <- c(
  "n_sample",
  "n_analytic",
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
  "agg_fun"
)
call2_names <- c(
  "",
  "x",
  "y",
  "conf_level",
  "conf_method",
  "n_resamples"
)

# wide format
expect_silent(res <- rdpmedian(data = data_wide, formula = a ~ b))
expect_true(inherits(res, "rdpmedian"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call_names)

expect_silent(res <- rdpmedian2(x = data_wide$a, y = data_wide$b))
expect_true(inherits(res, "rdpmedian"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call2_names)

# tall format
expect_silent(
  res <- rdpmedian(data = data_tall, formula = y ~ x | z)
)
expect_true(inherits(res, "rdpmedian"))
expect_true(inherits(res, "list"))
expect_identical(names(res), lst_names)
expect_identical(names(res$info), info_names)
expect_identical(names(res$call), call_names)

# Check pseudomedian
expect_true(is.numeric(res$pseudomedian))

#-------------------------------------------------------------------------------
# Handles NA values
#-------------------------------------------------------------------------------
data_na <- data_wide
data_na$a[1] <- NA
expect_silent(rdpmedian(data = data_na, formula = a ~ b))

#-------------------------------------------------------------------------------
# rdpmedian is same as rdpmedian2
#-------------------------------------------------------------------------------
res1 <- rdpmedian(data = data_wide, formula = a ~ b)
res2 <- rdpmedian2(x = data_wide$a, y = data_wide$b)
res1$call$data <- res1$call$formula <- res1$call[[1]] <- NULL
res1$info$focal_name <- res1$info$reference_name <- NULL
res1$call$agg_fun <- NULL
res1$lower <- res1$upper <- NULL
res2$call$x <- res2$call$y <- res2$call[[1]] <- NULL
res2$info$focal_name <- res2$info$reference_name <- NULL
res2$lower <- res2$upper <- NULL
expect_equal(res1, res2, check.attributes = FALSE)

#-------------------------------------------------------------------------------
# Confidence interval
#-------------------------------------------------------------------------------
res_ci <- rdpmedian(data = data_wide, formula = a ~ b, conf_level = 0.95)
expect_true(is.numeric(res_ci$lower) && !is.na(res_ci$lower))
expect_true(is.numeric(res_ci$upper) && !is.na(res_ci$upper))

# Alternative conf_method options
res_cm <- rdpmedian(
  data = data_wide,
  formula = a ~ b,
  conf_method = "percentile"
)
expect_equal("percentile", res_cm$call$conf_method)
res_cm <- rdpmedian(data = data_wide, formula = a ~ b, conf_method = "bca")
expect_equal("bca", res_cm$call$conf_method)

# invalid conf_level
expect_error(rdpmedian(data = data_wide, formula = a ~ b, conf_level = 2))

# invalid conf_method
expect_error(rdpmedian(
  data = data_wide,
  formula = a ~ b,
  conf_method = "invalid"
))

# invalid n_resamples
expect_error(rdpmedian(data = data_wide, formula = a ~ b, n_resamples = 9))

#-------------------------------------------------------------------------------
# Same results after reordering
#-------------------------------------------------------------------------------
res1 <- rdpmedian(data = data_wide, formula = a ~ b)
res2 <- rdpmedian(data = data_wide[sample(nrow(data_wide)), ], formula = a ~ b)
res1$lower <- res1$upper <- NULL
res2$lower <- res2$upper <- NULL
res1$call$data <- NULL
res2$call$data <- NULL
expect_identical(res1, res2)
