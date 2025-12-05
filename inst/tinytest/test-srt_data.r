library(tinytest)
library(rankdifferencetest)

#-------------------------------------------------------------------------------
# agg_fun handles duplicate group/block
#-------------------------------------------------------------------------------
# All NA
dup <- reshape(
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
dup$treatment <- factor(dup$treatment)
dup$subject <- factor(dup$subject)
dup[1, "time"] <- NA
dup[27, ] <- list(1, "placebo", 1)
dup[28, ] <- list(1, "placebo", 2)

res <- rankdifferencetest::srt_data(
  dup,
  time ~ treatment | subject,
  agg_fun = "first"
)
expect_equal(res$diffs[1], 1 - 2.9)
res <- rankdifferencetest::srt_data(
  dup,
  time ~ treatment | subject,
  agg_fun = "last"
)
expect_equal(res$diffs[1], 2 - 2.9)
res <- rankdifferencetest::srt_data(
  dup,
  time ~ treatment | subject,
  agg_fun = "sum"
)
expect_equal(res$diffs[1], 3 - 2.9)
res <- rankdifferencetest::srt_data(
  dup,
  time ~ treatment | subject,
  agg_fun = "mean"
)
expect_equal(res$diffs[1], 1.5 - 2.9)
res <- rankdifferencetest::srt_data(
  dup,
  time ~ treatment | subject,
  agg_fun = "median"
)
expect_equal(res$diffs[1], 1.5 - 2.9)
res <- rankdifferencetest::srt_data(
  dup,
  time ~ treatment | subject,
  agg_fun = "min"
)
expect_equal(res$diffs[1], 1 - 2.9)
res <- rankdifferencetest::srt_data(
  dup,
  time ~ treatment | subject,
  agg_fun = "max"
)
expect_equal(res$diffs[1], 2 - 2.9)
