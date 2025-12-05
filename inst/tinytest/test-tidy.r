library(rankdifferencetest)
library(tinytest)

set.seed(123)
data_wide <- data.frame(a = rnorm(10), b = rnorm(10))
data_tall <- data.frame(
  y = rnorm(20),
  x = factor(rep(c("grp1", "grp2"), each = 10)),
  z = factor(rep(1:10, times = 2))
)

pmedian_names <- c("estimate", "conf.low", "conf.high", "conf.level", "method")
pmedian_noci_names <- c("estimate", "method")
srt_names <- c(
  "estimate",
  "statistic",
  "p.value",
  "conf.low",
  "conf.high",
  "conf.level",
  "method",
  "alternative"
)
srt_noci_names <- c(
  "estimate",
  "statistic",
  "p.value",
  "method",
  "alternative"
)

#-------------------------------------------------------------------------------
# pmedian
#-------------------------------------------------------------------------------
df_pmedian1 <- tidy(pmedian(data = data_wide, formula = ~b))
df_pmedian <- tidy(pmedian(data = data_wide, formula = a ~ b))
df_pmedian2 <- tidy(pmedian2(data_wide$a, data_wide$b))

expect_equal(names(df_pmedian1), pmedian_names)
expect_equal(names(df_pmedian), pmedian_names)
expect_equal(names(df_pmedian2), pmedian_names)

df_pmedian <- tidy(pmedian(data = data_wide, formula = a ~ b, 0))
expect_equal(names(df_pmedian), pmedian_noci_names)

#-------------------------------------------------------------------------------
# rdpmedian
#-------------------------------------------------------------------------------
df_rdpmedian <- tidy(rdpmedian(data = data_wide, formula = a ~ b))
df_rdpmedian2 <- tidy(rdpmedian2(data_wide$a, data_wide$b))

expect_equal(names(df_rdpmedian), pmedian_names)
expect_equal(names(df_rdpmedian2), pmedian_names)

df_rdpmedian <- tidy(rdpmedian(data = data_wide, formula = a ~ b, 0))
expect_equal(names(df_rdpmedian), pmedian_noci_names)

#-------------------------------------------------------------------------------
# rdt
#-------------------------------------------------------------------------------
df_rdt <- tidy(rdt(
  data = data_wide,
  formula = a ~ b,
  conf_level = 0.95
))
df_rdt2 <- tidy(rdt2(data_wide$a, data_wide$b, conf_level = 0.95))

expect_equal(names(df_rdt), srt_names)
expect_equal(names(df_rdt2), srt_names)

df_rdt <- tidy(rdt(data = data_wide, formula = a ~ b))
expect_equal(names(df_rdt), srt_noci_names)

#-------------------------------------------------------------------------------
# srt
#-------------------------------------------------------------------------------
df_srt1 <- tidy(srt(
  data = data_wide,
  formula = ~b,
  conf_level = 0.95
))
df_srt <- tidy(srt(
  data = data_wide,
  formula = a ~ b,
  conf_level = 0.95
))
df_srt2 <- tidy(srt2(data_wide$a, data_wide$b, conf_level = 0.95))

expect_equal(names(df_srt1), srt_names)
expect_equal(names(df_srt), srt_names)
expect_equal(names(df_srt2), srt_names)

df_srt <- tidy(srt(data = data_wide, formula = a ~ b))
expect_equal(names(df_srt), srt_noci_names)
