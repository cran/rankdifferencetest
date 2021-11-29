library(testit)

# Data
data <- kornbrot_table1
data_long <- reshape(
  data = kornbrot_table1,
  direction = "long",
  varying = c("placebo", "drug"),
  v.names = c("time"),
  idvar = "subject",
  times = c("placebo", "drug"),
  timevar = "treatment",
  new.row.names = seq_len(prod(length(c("placebo", "drug")), nrow(kornbrot_table1)))
)
data_long$subject <- factor(data_long$subject)
data_long$treatment <- factor(data_long$treatment, levels = c("placebo", "drug"))

#-------------------------------------------------------------------------------
# Argument checks return errors
#-------------------------------------------------------------------------------
assert(
  "rdt() errors if 'data' is not data.frame.",
  has_error(rdt(data = 1, formula = placebo ~ drug))
)
assert(
  "rdt() errors if 'formula' is not a formula.",
  has_error(rdt(data = data, formula = data))
)
assert(
  "rdt() errors if 'zero.method' is incorrect.",
  has_error(rdt(data = data, formula = placebo ~ drug, zero.method = "wilcoxon"))
)
assert(
  "rdt() errors if 'distribution' is incorrect.",
  has_error(rdt(data = data, formula = placebo ~ drug, distribution = "exactt"))
)
assert(
  "rdt() errors if 'alternative' is incorrect.",
  has_error(rdt(data = data, formula = placebo ~ drug, alternative = "lesser"))
)
assert(
  "rdt() errors if variable names not found in 'data'.",
  has_error(rdt(data = data, formula = placebo ~ drug2))
)
assert(
  "rdt() errors if y~x formula is misformed.",
  has_error(rdt(data = data, formula = placebo ~ drug + subject))
)
assert(
  "rdt() errors if y~x|block formula has extra variable.", {

  tmp <- data_long
  tmp$xtra <- seq_len(nrow(data_long))

  has_error(rdt(data = tmp, formula = time ~ treatment + xtra | subject))
})
assert(
  "rdt() errors if y~x|block formula uses wrong block.", {

  tmp <- data_long
  tmp$xtra <- seq_len(nrow(data_long))

  has_error(rdt(data = tmp, formula = time ~ treatment | xtra))
})
assert(
  "rdt() errors if y~x|block formula has more than 2 blocking factors.", {

  tmp <- data_long
  tmp$treatment <- factor(c("a", rep("b", 12), rep("c", 13)))

  has_error(rdt(data = tmp, formula = time ~ treatment | treatment))
})

#-------------------------------------------------------------------------------
# p-value matches Kornbrot 1990
#-------------------------------------------------------------------------------
assert(
  "rdt() returns expected p-value for wide format data.", {

  res <- rdt(
    data = data,
    formula = placebo ~ drug,
    alternative = "greater",
    distribution = "asymptotic"
  )
  isTRUE(all.equal(res$p.value, 0.1238529, tolerance = 1e-6))
})

assert(
  "rdt() returns expected p-value for transformed wide format data.", {

  res <- rdt(
    data = data,
    formula = I(60/placebo) ~ I(60/drug),
    alternative = "less",
    distribution = "asymptotic"
  )
  isTRUE(all.equal(res$p.value, 0.1238529, tolerance = 1e-6))
})

assert(
  "rdt() returns expected p-value for long format data.", {

  res <- rdt(
    data = data_long,
    formula = time ~ treatment | subject,
    alternative = "greater",
    distribution = "asymptotic"
  )
  isTRUE(all.equal(res$p.value, 0.1238529, tolerance = 1e-6))
})

assert(
  "rdt() returns expected p-value for transformed long format data.", {

  res <- rdt(
    data = data_long,
    formula = I(60/time) ~ treatment | subject,
    alternative = "less",
    distribution = "asymptotic"
  )
  isTRUE(all.equal(res$p.value, 0.1238529, tolerance = 1e-6))
})

#-------------------------------------------------------------------------------
# Blocking handles random rows
#-------------------------------------------------------------------------------
assert(
  "rdt() blocking correctly handles random rows.", {

  data_long2 <- data_long[sample(seq_len(nrow(data_long)), nrow(data_long)), ]

  res1 <- rdt(
    data = data_long,
    formula = time ~ treatment | subject
  )
  res2 <- rdt(
    data = data_long2,
    formula = time ~ treatment | subject
  )

  isTRUE(all.equal(res1$p.value, res2$p.value))
})

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
# 4. zero.method
#     - Wilcoxon
#     - Pratt
#-------------------------------------------------------------------------------
assert(
  "rdt() returns correct formula.", {

  # Formula
  res <- rdt(
    data = data,
    formula = placebo ~ drug
  )
  res2 <- rdt(
    data = data_long,
    formula = time ~ treatment | subject
  )

  res$formula == "placebo ~ drug" &
    res2$formula == "time ~ treatment | subject"
})

assert(
  "rdt() returns correct alternative.", {

  # Formula
  res <- rdt(
    data = data,
    formula = placebo ~ drug,
    alternative = "two.sided"
  )
  res2 <- rdt(
    data = data,
    formula = placebo ~ drug,
    alternative = "greater"
  )
  res3 <- rdt(
    data = data,
    formula = placebo ~ drug,
    alternative = "less"
  )

  res$alternative == "True location shift of ranks (placebo - drug) is not equal to 0" &
    res2$alternative == "True location shift of ranks (placebo - drug) is greater than 0" &
    res3$alternative == "True location shift of ranks (placebo - drug) is less than 0"
})

assert(
  "rdt() returns correct distribution.", {

  # Formula
  res <- rdt(
    data = data,
    formula = placebo ~ drug,
    distribution = "exact"
  )
  res2 <- rdt(
    data = data,
    formula = placebo ~ drug,
    distribution = "asymptotic"
  )
  res3 <- rdt(
    data = data,
    formula = placebo ~ drug,
    distribution = "approximate"
  )

  res$method == "Kornbrot's Rank Difference Test using the Exact Wilcoxon-Pratt Signed-Rank Test" &
    res2$method == "Kornbrot's Rank Difference Test using the Asymptotic Wilcoxon-Pratt Signed-Rank Test" &
    res3$method == "Kornbrot's Rank Difference Test using the Approximate Wilcoxon-Pratt Signed-Rank Test"
})

assert(
  "rdt() returns correct zero-difference method.", {

  # Formula
  res <- rdt(
    data = data,
    formula = placebo ~ drug,
    zero.method = "Wilcoxon"
  )
  res2 <- rdt(
    data = data,
    formula = placebo ~ drug,
    zero.method = "Pratt"
  )

  res$method == "Kornbrot's Rank Difference Test using the Asymptotic Wilcoxon Signed-Rank Test" &
    res2$method == "Kornbrot's Rank Difference Test using the Asymptotic Wilcoxon-Pratt Signed-Rank Test"
})

