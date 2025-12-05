#-------------------------------------------------------------------------------
# Below we verify the calculation of the test statistic.
# In conclusion, calculate test statistic components using the sum of the
# ranks with averages for ties. It produces the same result and is simpler.
#-------------------------------------------------------------------------------
f1 <- function(x, y, zero_method, correct) {
  call <- list(
    mu = 0,
    zero_method = zero_method,
    digits_rank = Inf,
    correct = correct,
    alternative = "two.sided"
  )
  data <- rankdifferencetest:::srt_data.numeric(x, y, "x", "y")
  data <- rankdifferencetest:::srt_ranks(data, call)
  data <- rankdifferencetest:::srt_method(
    data,
    test_class = "asymptotic"
  )
  data <- rankdifferencetest:::srt_statistic(data)
  data$statistic
}

f2 <- function(x, y, zero_method, correct) {
  call <- list(
    mu = 0,
    zero_method = zero_method,
    digits_rank = Inf,
    correct = correct,
    alternative = "two.sided"
  )
  data <- rankdifferencetest:::srt_data.numeric(x, y, "x", "y")
  data <- rankdifferencetest:::srt_ranks(data, call)
  data <- rankdifferencetest:::srt_method(
    data,
    test_class = "asymptotic"
  )

  ties <- rankdifferencetest:::ftabulate(data$ranks) - 1L

  if (zero_method == "wilcoxon") {
    wplus <- data$wplus
    n <- length(data$ranks)
    ewplus <- (n * (n + 1)) / 4
    vwplus <- (n * (n + 1) * (2 * n + 1)) / 24
    t <- sum(ties^3 - ties) / 48
    statistic <- (wplus - ewplus) / sqrt(vwplus - t)
  }

  if (zero_method == "pratt") {
    wplus <- data$wplus
    n <- data$n_analytic
    nz <- data$n_zeros
    ewplus <- ((n * (n + 1)) - (nz * (nz + 1))) / 4
    vwplus <- (n * (n + 1) * (2 * n + 1) - nz * (nz + 1) * (2 * nz + 1)) / 24
    t <- sum(ties^3 - ties) / 48
    statistic <- (wplus - ewplus) / sqrt(vwplus - t)
  }

  correction <- if (correct) {
    0.5 * sign(wplus - ewplus)
  } else {
    0
  }

  # Test statistic
  z <- (wplus - ewplus - correction) / sqrt(vwplus)
  statistic <- setNames(z, "Z")
  statistic
}

# no zeros no ties
x <- rnorm(100)
y <- rnorm(100)

# 50% zeros no ties
z <- rnorm(50)
x <- c(rnorm(50), z)
y <- c(rnorm(50), z)

# no zeros 30% single-ties
x <- c(rnorm(70))
y <- c(rnorm(70))
z <- sample(1:70, 30)
x <- c(x, x[z])
y <- c(y, y[z])
# sum(duplicated(rank(abs(x - y))))

# no zeros 50% 5-ties
x <- c(rnorm(60))
y <- c(rnorm(60))
z <- sample(1:60, 10)
x <- c(x, rep(x[z], 5))
y <- c(y, rep(y[z], 5))
# sum(duplicated(rank(abs(x - y))))

# 20% zeros and 20% 2-ties
z <- rnorm(20)
x <- c(rnorm(60), z)
y <- c(rnorm(60), z)

z <- sample(1:60, 10)
x <- c(x, rep(x[z], 2))
y <- c(y, rep(y[z], 2))
# diffs <- x - y
# sum(diffs == 0)
# sum(duplicated(rank(abs(diffs))[diffs != 0]))

df <- data.frame(
  y = c(x, y),
  group = factor(c(rep("x", length(x)), rep("y", length(y)))),
  block = factor(c(seq_along(x), seq_along(y)))
)

# wilcoxon and correct
f1(x, y, "wilcoxon", TRUE)
f2(x, y, "wilcoxon", TRUE)

# wilcoxon and not correct
f1(x, y, "wilcoxon", FALSE)
f2(x, y, "wilcoxon", FALSE)

coin::wilcoxsign_test(
  formula = y ~ group | block,
  data = df,
  paired = TRUE,
  zero.method = "Wilcoxon"
)

# pratt and correct
f1(x, y, "pratt", TRUE)
f2(x, y, "pratt", TRUE)

# pratt and not correct
f1(x, y, "pratt", FALSE)
f2(x, y, "pratt", FALSE)

coin::wilcoxsign_test(
  formula = y ~ group | block,
  data = df,
  paired = TRUE,
  zero.method = "Pratt"
)

# # benchmark
# bench::mark(
#   f1(x, y, "wilcoxon", TRUE),
#   f2(x, y, "wilcoxon", TRUE),
#   check = FALSE
# )
#
# bench::mark(
#   f1(x, y, "pratt", TRUE),
#   f2(x, y, "pratt", TRUE),
#   check = FALSE
# )
