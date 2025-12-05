
# rankdifferencetest

## Overview

rankdifferencetest provides a modified Wilcoxon signed-rank test which
produces consistent and meaningful results for ordinal or monotonically
transformed data. The procedure was described by Diana Kornbrot in
<https://doi.org/10.1111/j.2044-8317.1990.tb00939.x>

## Installation

``` r
# Install from CRAN
install.packages("rankdifferencetest")

# Or the development version from bitbucket
remotes::install_bitbucket("bklamer/rankdifferencetest")
```

## Examples

``` r
# See ?rdt and ?kornbrot_table1
library(rankdifferencetest)

res <- rdt(
  data = kornbrot_table1,
  formula = placebo ~ drug,
  conf_level = 0.95,
  alternative = "two.sided",
  distribution = "asymptotic",
  zero_method = "wilcoxon",
  correct = TRUE
)
res
#> 
#> asymptotic Kornbrot-Wilcoxon rank difference test
#> 
#> p = 0.262
#> Z = 1.12
#> Paired differences of ranks: placebo - drug
#> Alternative hypothesis: True location shift of paired
#> differences of ranks is not equal to 0
#> 
#> Pseudomedian: 1.5
#> 95% CI: -2.5, 6
```

## License

Copyright: Brett Klamer - 2021 - GPL-3
(<https://opensource.org/license/gpl-3-0>)
