library(tinytest)

set.seed(123)
data <- data.frame(y = rnorm(6), x = rnorm(6))
datablock <- data.frame(
  y = rnorm(6),
  x = gl(2, 3, labels = c("pre", "post")),
  xx = factor(c(rep("pre", 2), rep("post", 2), rep("post2", 2))),
  b = rep(1:3, times = 2)
)

expect_error(
  rankdifferencetest:::rdt_data(x = data, formula = ~x),
  pattern = "The formula is missing the left hand side."
)

expect_warning(
  rankdifferencetest:::rdt_data(x = datablock, formula = y ~ x | b),
  pattern = "'b' was converted to a factor in formula 'y ~ x | b'."
)

expect_error(
  rankdifferencetest:::rdt_data(x = datablock, formula = y ~ xx | x),
  pattern = "'xx' must be a two-level factor in formula 'y ~ xx | x'."
)

expect_error(
  rankdifferencetest:::rdt_data(x = datablock, formula = y ~ x | block),
  pattern = "The formula included terms not found in the data: 'block'"
)

expect_error(
  rankdifferencetest:::rdt_data(x = datablock, formula = y ~ x),
  pattern = "'x' must be numeric in formula 'y ~ x'."
)
