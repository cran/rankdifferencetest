library(rankdifferencetest)

# srt() is faster than all reference implementations.
x <- 1:10
y <- 10:1

bench::mark(
  srt2_exact = srt2(x = x, y = y, distribution = "exact"),
  srt2_asymptotic = srt2(x = x, y = y, distribution = "asymptotic"),
  wilcox.test_asymptotic = wilcox.test(x = x, y = y),
  wilcox.exact_exact = exactRankTests::wilcox.exact(x = x, y = y),
  check = F
)

## A tibble: 4 × 13
#  expression                  min    median  `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
#  1 srt2_exact                169µs  179µs     5276.   25.97KB    10.5   2501     5      474ms
#  2 srt2_asymptotic           130µs  137µs     6804.   16.89KB    10.5   3243     5      477ms
#  3 wilcox.test_asymptotic    192µs  205µs     4818.    2.48KB    10.5   2300     5      477ms
#  4 wilcox.exact_exact        592µs  617µs     1620.   34.65KB     4.49   721     2      445ms

bench::mark(
  srt2_exact = srt2(x = x, y = y, distribution = "exact", conf_level = 0.95),
  srt2_asymptotic = srt2(
    x = x,
    y = y,
    distribution = "asymptotic",
    conf_level = 0.95
  ),
  wilcox.test_asymptotic = wilcox.test(x = x, y = y, conf.int = TRUE),
  wilcox.exact_exact = exactRankTests::wilcox.exact(
    x = x,
    y = y,
    conf.int = TRUE
  ),
  check = F
)

## A tibble: 4 × 13
#expression                   min    median  `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
#  1 srt2_exact               1.25ms 1.31ms      760.    44.9KB    15.0    354     7      466ms
#  2 srt2_asymptotic          1.92ms 1.99ms      471.     1.2KB    15.3    216     7      458ms
#  3 wilcox.test_asymptotic   4.96ms 5.09ms      196.   189.5KB    11.0     89     5      454ms
#  4 wilcox.exact_exact       3.17ms 3.22ms      309.   131.9KB     8.51   145     4      470ms

data <- kornbrot_table1

bench::mark(
  srt2_exact = srt(
    data = data,
    formula = placebo ~ drug,
    distribution = "exact",
    zero_method = "pratt",
    correct = FALSE
  ),
  srt2_asymptotic = srt(
    data = data,
    formula = placebo ~ drug,
    distribution = "asymptotic",
    zero_method = "pratt",
    correct = FALSE
  ),
  wilcoxsign_test_exact = coin::wilcoxsign_test(
    formula = placebo ~ drug,
    data = data,
    distribution = "exact"
  ),
  wilcoxsign_test_asymptotic = coin::wilcoxsign_test(
    formula = placebo ~ drug,
    data = data,
    distribution = "asymptotic"
  ),
  check = F
)

## A tibble: 4 × 13
#expression                min     median    `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
#  1 srt2_exact           159.83µs 173.57µs     5674.     458KB     4.62  2455     2      433ms
#  2 srt2_asymptotic      155.29µs 165.76µs     5858.    16.9KB     6.25  2810     3      480ms
#  3 wilcoxsign_test_exa…   4.18ms   4.32ms      227.    80.8MB     6.44   106     3      466ms
#  4 wilcoxsign_test_asy…   3.42ms    3.5ms      281.   228.7KB     6.43   131     3      467ms
