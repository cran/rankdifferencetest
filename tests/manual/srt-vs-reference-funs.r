library(rankdifferencetest)
library(tinytest)

#-------------------------------------------------------------------------------
# Helpers
#-------------------------------------------------------------------------------
.srt_vs_wilcox.test <- function(s, w) {
  all(
    c(
      #as.numeric(s$statistic) == as.numeric(w$statistic),
      isTRUE(all.equal(
        s$p_value,
        w$p.value,
        tolerance = 1e-3,
        check.attributes = FALSE
      )),
      s$call$mu == w$null.value,
      s$alternative == w$alternative,
      isTRUE(all.equal(
        s$lower,
        w$conf.int[1],
        tolerance = 1e-3,
        check.attributes = FALSE
      )),
      isTRUE(all.equal(
        s$upper,
        w$conf.int[2],
        tolerance = 1e-3,
        check.attributes = FALSE
      )),
      isTRUE(all.equal(
        s$pseudomedian,
        as.numeric(w$estimate),
        tolerance = 1e-3
      ))
    )
  )
}

.srt_vs_wilcox.exact <- function(s, w) {
  all(
    c(
      #as.numeric(s$statistic) == as.numeric(w$statistic),
      isTRUE(all.equal(
        s$p_value,
        w$p.value,
        tolerance = 1e-3,
        check.attributes = FALSE
      )),
      s$call$mu == w$null.value,
      s$alternative == w$alternative,
      isTRUE(all.equal(
        s$lower,
        w$conf.int[1],
        tolerance = 1e-3,
        check.attributes = FALSE
      )),
      isTRUE(all.equal(
        s$upper,
        w$conf.int[2],
        tolerance = 1e-3,
        check.attributes = FALSE
      )),
      isTRUE(all.equal(
        s$pseudomedian,
        as.numeric(w$estimate),
        tolerance = 1e-3
      ))
    )
  )
}

# wilcoxsign_test differences are reversed from srt.
# wilcoxsign_test always uses correct=FALSE.
.srt_vs_wilcoxsign_test <- function(s, w) {
  all(
    c(
      if (s$call$distribution == "asymptotic") {
        isTRUE(all.equal(
          s$statistic,
          as.numeric(coin::statistic(w)),
          tolerance = 1e-3,
          check.attributes = FALSE
        ))
      } else {
        NULL
      },
      isTRUE(all.equal(
        s$p_value,
        as.numeric(coin::pvalue(w)),
        tolerance = 1e-3,
        check.attributes = FALSE
      )),
      isTRUE(all.equal(
        s$call$mu,
        w@nullvalue,
        tolerance = 1e-3,
        check.attributes = FALSE
      ))
    )
  )
}

.srt_vs_srt2 <- function(s, s2) {
  s$call <- NULL
  s$info$reference_name <- NULL
  s$info$focal_name <- NULL
  s2$call <- NULL
  s2$info$reference_name <- NULL
  s2$info$focal_name <- NULL
  isTRUE(all.equal(s, s2))
}

.vec2wide <- function(x, y) {
  x_name <- deparse(substitute(x))
  y_name <- deparse(substitute(y))
  d <- data.frame(x, y)
  names(d) <- c(x_name, y_name)
  d
}

.vec2tall <- function(x, y) {
  x_name <- deparse(substitute(x))
  y_name <- deparse(substitute(y))
  data.frame(
    value = c(y, x),
    group = factor(
      c(rep(y_name, length(y)), rep(x_name, length(y))),
      levels = c(y_name, x_name)
    ),
    block = factor(c(seq_along(y), seq_along(x)))
  )
}

.rows_to_list <- function(x) {
  res <- lapply(seq_len(nrow(x)), function(i) {
    vals <- as.list(x[i, , drop = FALSE])
    names(vals) <- names(x)
    out <- as.pairlist(vals)
    attributes(out) <- NULL
    names(out) <- names(x)
    out
  })

  if (!is.null(res[[1]]$data)) {
    res <- lapply(res, function(x) {
      x$data <- eval(x$data[[1]])
      x
    })
  } else {
    res <- lapply(res, function(x) {
      x$x <- eval(x$x[[1]])
      x$y <- eval(x$y[[1]])
      x
    })
  }

  if (!is.null(res[[1]]$formula)) {
    res <- lapply(res, function(x) {
      x$formula <- x$formula[[1]]
      x
    })
  }

  lapply(res, as.list)
}

#-------------------------------------------------------------------------------
# Simulate data
#-------------------------------------------------------------------------------
# Data
# no zeros no ties
n <- 40
set.seed(20251128)
x <- rnorm(n)
y <- rnorm(n)
d_tall_nznt <- .vec2tall(x, y)
d_wide_nznt <- .vec2wide(x, y)

# 50% zeros no ties
zprop <- 0.5
z <- rnorm(n * zprop)
x <- c(rnorm(n * (1 - zprop)), z)
y <- c(rnorm(n * (1 - zprop)), z)
d_tall_5znt <- .vec2tall(x, y)
d_wide_5znt <- .vec2wide(x, y)

# no zeros 30% single-ties
tprop <- 0.3
x <- c(rnorm(n * (1 - tprop)))
y <- c(rnorm(n * (1 - tprop)))
z <- sample(1:(n * (1 - tprop)), n * tprop)
x <- c(x, x[z])
y <- c(y, y[z])
# sum(duplicated(rank(abs(x - y))))
d_tall_nz3t <- .vec2tall(x, y)
d_wide_nz3t <- .vec2wide(x, y)

# no zeros 50% 5-ties
tprop <- 0.5
nties <- 5
x <- c(rnorm(n - (n * tprop) / nties))
y <- c(rnorm(n - (n * tprop) / nties))
z <- sample(1:(n - (n * tprop) / nties), (n * tprop) / nties)
x <- c(x, rep(x[z], nties))
y <- c(y, rep(y[z], nties))
# sum(duplicated(rank(abs(x - y))))
d_tall_nz5t <- .vec2tall(x, y)
d_wide_nz5t <- .vec2wide(x, y)

# 20% zeros and 40% 2-ties
zprop <- 0.2
tprop <- 0.4
nties <- 2

z <- rnorm(n * zprop)
x <- rnorm(n * (1 - (zprop + tprop)))
y <- rnorm(n * (1 - (zprop + tprop)))

t <- sample(seq_along(x), (n * tprop) / nties)
x <- c(x, z, rep(x[t], nties))
y <- c(y, z, rep(y[t], nties))
# diffs <- x - y
# sum(diffs == 0)
# sum(duplicated(rank(abs(diffs))[diffs != 0]))
d_tall_2z4t <- .vec2tall(x, y)
d_wide_2z4t <- .vec2wide(x, y)

# Reverse factor levels for wilcoxsign_test
d_tall_nzntr <- d_tall_nznt
d_tall_nzntr$group <- factor(
  d_tall_nzntr$group,
  levels = rev(levels(d_tall_nznt$group))
)
d_tall_5zntr <- d_tall_5znt
d_tall_5zntr$group <- factor(
  d_tall_5zntr$group,
  levels = rev(levels(d_tall_5znt$group))
)
d_tall_nz3tr <- d_tall_nz3t
d_tall_nz3tr$group <- factor(
  d_tall_nz3tr$group,
  levels = rev(levels(d_tall_nz3t$group))
)
d_tall_nz5tr <- d_tall_nz5t
d_tall_nz5tr$group <- factor(
  d_tall_nz5tr$group,
  levels = rev(levels(d_tall_nz5t$group))
)
d_tall_2z4tr <- d_tall_2z4t
d_tall_2z4tr$group <- factor(
  d_tall_2z4tr$group,
  levels = rev(levels(d_tall_2z4t$group))
)

#-------------------------------------------------------------------------------
# Argument grids
#-------------------------------------------------------------------------------
## srt
srt_args <- expand.grid(
  data = alist(
    d_wide_nznt,
    d_tall_nznt,
    d_wide_5znt,
    d_tall_5znt,
    d_wide_nz3t,
    d_tall_nz3t,
    d_wide_nz5t,
    d_tall_nz5t,
    d_wide_2z4t,
    d_tall_2z4t
  ),
  formula = c(~x, x ~ y, value ~ group | block),
  conf_level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  distribution = c("exact", "asymptotic"),
  correct = c(TRUE, FALSE),
  zero_method = c("wilcoxon", "pratt"),
  stringsAsFactors = FALSE
) |>
  as.data.frame()

ok <- !((startsWith(as.character(srt_args$data), "d_tall") &
  srt_args$formula == ~x) |
  (startsWith(as.character(srt_args$data), "d_wide") &
    srt_args$formula == "value ~ group | block") |
  (startsWith(as.character(srt_args$data), "d_tall") &
    srt_args$formula == "x ~ y") |
  srt_args$distribution == "exact" & srt_args$correct)
srt_args <- srt_args[ok, ]
rownames(srt_args) <- NULL
srt_args <- srt_args[
  order(
    srt_args$zero_method,
    srt_args$correct,
    srt_args$distribution,
    srt_args$alternative,
    decreasing = TRUE
  ),
]

## srt2
srt2_args <- expand.grid(
  x = alist(
    d_wide_nznt[["x"]],
    d_wide_5znt[["x"]],
    d_wide_nz3t[["x"]],
    d_wide_nz5t[["x"]],
    d_wide_2z4t[["x"]]
  ),
  y = alist(
    NULL,
    d_wide_nznt[["y"]],
    d_wide_5znt[["y"]],
    d_wide_nz3t[["y"]],
    d_wide_nz5t[["y"]],
    d_wide_2z4t[["y"]]
  ),
  conf_level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  distribution = c("exact", "asymptotic"),
  correct = c(TRUE, FALSE),
  zero_method = c("wilcoxon", "pratt"),
  stringsAsFactors = FALSE
) |>
  as.data.frame()

ok <- !((startsWith(as.character(srt2_args$x), "d_wide_2z4t") &
  !(startsWith(as.character(srt2_args$y), "d_wide_2z4t") |
    startsWith(as.character(srt2_args$y), "NULL"))) |
  (startsWith(as.character(srt2_args$x), "d_wide_5znt") &
    !(startsWith(as.character(srt2_args$y), "d_wide_5znt") |
      startsWith(as.character(srt2_args$y), "NULL"))) |
  (startsWith(as.character(srt2_args$x), "d_wide_nz3t") &
    !(startsWith(as.character(srt2_args$y), "d_wide_nz3t") |
      startsWith(as.character(srt2_args$y), "NULL"))) |
  (startsWith(as.character(srt2_args$x), "d_wide_nz5t") &
    !(startsWith(as.character(srt2_args$y), "d_wide_nz5t") |
      startsWith(as.character(srt2_args$y), "NULL"))) |
  (startsWith(as.character(srt2_args$x), "d_wide_nznt") &
    !(startsWith(as.character(srt2_args$y), "d_wide_nznt") |
      startsWith(as.character(srt2_args$y), "NULL"))) |
  srt2_args$distribution == "exact" & srt2_args$correct)
srt2_args <- srt2_args[ok, ]
rownames(srt2_args) <- NULL

# just double to match srt_args
srt2_nulls <- as.character(srt2_args$y) == "NULL"
srt2_args2 <- do.call(
  "rbind",
  replicate(2L, srt2_args[!srt2_nulls, ], simplify = FALSE)
)
srt2_args <- rbind(srt2_args[srt2_nulls, ], srt2_args2)
srt2_args <- srt2_args[
  order(
    srt2_args$zero_method,
    srt2_args$correct,
    srt2_args$distribution,
    srt2_args$alternative,
    decreasing = TRUE
  ),
]
srt2_nulls <- as.character(srt2_args$y) == "NULL"

# wilcox.test
wilcox.test_args <- expand.grid(
  x = alist(
    d_wide_nznt[["x"]],
    d_wide_5znt[["x"]],
    d_wide_nz3t[["x"]],
    d_wide_nz5t[["x"]],
    d_wide_2z4t[["x"]]
  ),
  y = alist(
    d_wide_nznt[["y"]],
    d_wide_5znt[["y"]],
    d_wide_nz3t[["y"]],
    d_wide_nz5t[["y"]],
    d_wide_2z4t[["y"]]
  ),
  conf.int = TRUE,
  alternative = c("two.sided", "greater", "less"),
  exact = c(TRUE, FALSE),
  correct = c(TRUE, FALSE),
  paired = TRUE,
  stringsAsFactors = FALSE
) |>
  as.data.frame()

ok <- !((startsWith(as.character(wilcox.test_args$x), "d_wide_nznt") &
  !startsWith(as.character(wilcox.test_args$y), "d_wide_nznt")) |
  (startsWith(as.character(wilcox.test_args$x), "d_wide_5znt") &
    !startsWith(as.character(wilcox.test_args$y), "d_wide_5znt")) |
  (startsWith(as.character(wilcox.test_args$x), "d_wide_nz3t") &
    !startsWith(as.character(wilcox.test_args$y), "d_wide_nz3t")) |
  (startsWith(as.character(wilcox.test_args$x), "d_wide_nz5t") &
    !startsWith(as.character(wilcox.test_args$y), "d_wide_nz5t")) |
  (startsWith(as.character(wilcox.test_args$x), "d_wide_2z4t") &
    !startsWith(as.character(wilcox.test_args$y), "d_wide_2z4t")) |
  wilcox.test_args$exact & wilcox.test_args$correct)
wilcox.test_args <- wilcox.test_args[ok, ]
rownames(wilcox.test_args) <- NULL

wilcox.test_args <- wilcox.test_args[
  order(
    wilcox.test_args$correct,
    wilcox.test_args$exact,
    wilcox.test_args$alternative,
    decreasing = TRUE
  ),
]

# wilcox.exact
wilcox.exact_args <- wilcox.test_args
wilcox.exact_args <- wilcox.exact_args[!wilcox.exact_args$correct, ]
wilcox.exact_args$correct <- NULL

# wilcoxsign_test
# wilcoxsign_test differences are reversed from srt.
# wilcoxsign_test always uses correct=FALSE.
wilcoxsign_test_args <- expand.grid(
  formula = c(value ~ group | block),
  data = alist(
    d_tall_nzntr,
    d_tall_5zntr,
    d_tall_nz3tr,
    d_tall_nz5tr,
    d_tall_2z4tr
  ),
  alternative = c("less", "greater", "two.sided"),
  distribution = c("exact", "asymptotic"),
  zero.method = c("Wilcoxon", "Pratt"),
  stringsAsFactors = FALSE
) |>
  as.data.frame()

wilcoxsign_test_args <- wilcoxsign_test_args[
  order(
    wilcoxsign_test_args$zero.method,
    wilcoxsign_test_args$distribution,
    wilcoxsign_test_args$alternative,
    decreasing = TRUE
  ),
]

#-------------------------------------------------------------------------------
# Run and compare tests
#-------------------------------------------------------------------------------
# srt vs srt2
srt_args_lst <- .rows_to_list(srt_args)
srt2_args_lst <- .rows_to_list(srt2_args)

srt_res <- lapply(seq_along(srt_args_lst), function(x) {
  do.call(srt, args = srt_args_lst[[x]])
})

srt2_res <- lapply(seq_along(srt2_args_lst), function(x) {
  do.call(srt2, args = srt2_args_lst[[x]])
})

res <- mapply(.srt_vs_srt2, srt_res, srt2_res)

expect_true(all(res))

# wilcox.test vs srt2 (matching asymptotic tests)
ok <- !wilcox.test_args$exact
wilcox.test_lst <- .rows_to_list(wilcox.test_args[ok, ])
wilcox.test_res <- lapply(seq_along(wilcox.test_lst), function(x) {
  do.call(wilcox.test, args = wilcox.test_lst[[x]])
})

ok <- !duplicated(srt2_args) &
  !srt2_nulls &
  srt2_args$zero_method == "wilcoxon" &
  srt2_args$distribution == "asymptotic"
srt2_args_lst <- .rows_to_list(srt2_args[ok, ])

srt2_res <- lapply(seq_along(srt2_args_lst), function(x) {
  do.call(srt2, args = srt2_args_lst[[x]])
})

res <- mapply(.srt_vs_wilcox.test, srt2_res, wilcox.test_res)

# Every `FALSE` occurs because wilcox.test forces `correct=FALSE` in its
# pseudomedian estimate. `srt` allows `correct` to be set by user.
expect_true(sum(res) / length(res) == 0.7)

# wilcox.exact vs srt2 (matching exact tests)
ok <- wilcox.exact_args$exact
wilcox.exact_lst <- .rows_to_list(wilcox.exact_args[ok, ])
wilcox.exact_res <- lapply(seq_along(wilcox.exact_lst), function(x) {
  do.call(exactRankTests::wilcox.exact, args = wilcox.exact_lst[[x]])
})

ok <- !duplicated(srt2_args) &
  !srt2_nulls &
  srt2_args$zero_method == "wilcoxon" &
  !srt2_args$correct &
  srt2_args$distribution == "exact"
srt2_args_lst <- .rows_to_list(srt2_args[ok, ])

srt2_res <- lapply(seq_along(srt2_args_lst), function(x) {
  do.call(srt2, args = srt2_args_lst[[x]])
})

res <- mapply(.srt_vs_wilcox.exact, srt2_res, wilcox.exact_res)

expect_true(all(res))

# wilcoxsign_test vs srt2
wilcoxsign_test_lst <- .rows_to_list(wilcoxsign_test_args)
wilcoxsign_test_res <- lapply(seq_along(wilcoxsign_test_lst), function(x) {
  do.call(coin::wilcoxsign_test, args = wilcoxsign_test_lst[[x]])
})

ok <- !duplicated(srt2_args) & !srt2_nulls & !srt2_args$correct
srt2_args_lst <- .rows_to_list(srt2_args[ok, ])

srt2_res <- lapply(seq_along(srt2_args_lst), function(x) {
  do.call(srt2, args = srt2_args_lst[[x]])
})

res <- mapply(.srt_vs_wilcoxsign_test, srt2_res, wilcoxsign_test_res)

expect_true(all(res))
