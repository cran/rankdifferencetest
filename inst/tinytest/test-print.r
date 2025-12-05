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
# pmedian
#-------------------------------------------------------------------------------
print_pmedian1 <- capture.output(print(pmedian(
  data = data_wide,
  formula = ~b
)))
print_pmedian <- capture.output(print(pmedian(
  data = data_wide,
  formula = a ~ b
)))
print_pmedian2 <- capture.output(print(pmedian2(data_wide$a, data_wide$b)))
print_pmedian2_noci <- capture.output(print(pmedian2(
  data_wide$a,
  data_wide$b,
  conf_level = 0
)))

expect_equal(
  print_pmedian1,
  c(
    "",
    "Hodges-Lehmann estimator and percentile bootstrap confidence",
    "interval",
    "",
    "Observations: b",
    "",
    "Pseudomedian: 0.334",
    "95% CI: -0.473, 0.812"
  )
)
expect_equal(
  print_pmedian,
  c(
    "",
    "Hodges-Lehmann estimator and percentile bootstrap confidence",
    "interval",
    "",
    "Paired differences: a - b",
    "",
    "Pseudomedian: -0.0402",
    "95% CI: -0.714, 0.364"
  )
)
expect_equal(
  print_pmedian2,
  c(
    "",
    "Hodges-Lehmann estimator and percentile bootstrap confidence",
    "interval",
    "",
    "Paired differences: data_wide$a - data_wide$b",
    "",
    "Pseudomedian: -0.0402",
    "95% CI: -0.713, 0.364"
  )
)
expect_equal(
  print_pmedian2_noci,
  c(
    "",
    "Hodges-Lehmann estimator",
    "",
    "Paired differences: data_wide$a - data_wide$b",
    "",
    "Pseudomedian: -0.0402"
  )
)

#-------------------------------------------------------------------------------
# rdpmedian
#-------------------------------------------------------------------------------
print_rdpmedian <- capture.output(print(rdpmedian(
  data = data_wide,
  formula = a ~ b
)))
print_rdpmedian2 <- capture.output(print(rdpmedian2(data_wide$a, data_wide$b)))

expect_equal(
  print_rdpmedian,
  c(
    "",
    "Hodges-Lehmann estimator and percentile bootstrap confidence",
    "interval",
    "",
    "Paired differences of ranks: a - b",
    "",
    "Pseudomedian: -1",
    "95% CI: -7, 2"
  )
)
expect_equal(
  print_rdpmedian2,
  c(
    "",
    "Hodges-Lehmann estimator and percentile bootstrap confidence",
    "interval",
    "",
    "Paired differences of ranks: data_wide$a - data_wide$b",
    "",
    "Pseudomedian: -1",
    "95% CI: -7, 2.5"
  )
)

#-------------------------------------------------------------------------------
# rdt
#-------------------------------------------------------------------------------
print_rdt <- capture.output(print(rdt(
  data = data_wide,
  formula = a ~ b,
  conf_level = 0.95
)))
print_rdt2 <- capture.output(print(rdt2(
  data_wide$a,
  data_wide$b,
  conf_level = 0.95
)))

expect_equal(
  print_rdt,
  c(
    "",
    "exact Kornbrot-Wilcoxon rank difference test",
    "",
    "p = 0.516",
    "W+ = 21",
    "Paired differences of ranks: a - b",
    "Alternative hypothesis: True location shift of paired",
    "differences of ranks is not equal to 0",
    "",
    "Pseudomedian: -1.25",
    "94.8% CI: -7, 2.5"
  )
)
expect_equal(
  print_rdt2,
  c(
    "",
    "exact Kornbrot-Wilcoxon rank difference test",
    "",
    "p = 0.516",
    "W+ = 21",
    "Paired differences of ranks: data_wide$a - data_wide$b",
    "Alternative hypothesis: True location shift of paired",
    "differences of ranks is not equal to 0",
    "",
    "Pseudomedian: -1.25",
    "94.8% CI: -7, 2.5"
  )
)

#-------------------------------------------------------------------------------
# srt
#-------------------------------------------------------------------------------
print_srt1 <- capture.output(print(srt(
  data = data_wide,
  formula = ~b,
  conf_level = 0.95
)))
print_srt <- capture.output(print(srt(
  data = data_wide,
  formula = a ~ b,
  conf_level = 0.95
)))
print_srt2 <- capture.output(print(srt2(
  data_wide$a,
  data_wide$b,
  conf_level = 0.95
)))

expect_equal(
  print_srt1,
  c(
    "",
    "exact Wilcoxon signed-rank test",
    "",
    "p = 0.492",
    "W+ = 35",
    "Observations: b",
    "Alternative hypothesis: True location shift of observations",
    "is not equal to 0",
    "",
    "Pseudomedian: 0.334",
    "95% CI: -0.556, 0.949"
  )
)
expect_equal(
  print_srt,
  c(
    "",
    "exact Wilcoxon signed-rank test",
    "",
    "p = 0.625",
    "W+ = 22",
    "Paired differences: a - b",
    "Alternative hypothesis: True location shift of paired",
    "differences is not equal to 0",
    "",
    "Pseudomedian: -0.0402",
    "95% CI: -0.879, 0.559"
  )
)
expect_equal(
  print_srt2,
  c(
    "",
    "exact Wilcoxon signed-rank test",
    "",
    "p = 0.625",
    "W+ = 22",
    "Paired differences: data_wide$a - data_wide$b",
    "Alternative hypothesis: True location shift of paired",
    "differences is not equal to 0",
    "",
    "Pseudomedian: -0.0402",
    "95% CI: -0.879, 0.559"
  )
)
