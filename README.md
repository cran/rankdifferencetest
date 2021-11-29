# rankdifferencetest

`rankdifferencetest` provides a modified Wilcoxon signed-rank test which produces 
consistent and meaningful results for ordinal or monotonically transformed data.
The procedure was described by Diana Kornbrot in <https://doi.org/10.1111/j.2044-8317.1990.tb00939.x>

## Installation

To install from CRAN

```{r}
install.packages("rankdifferencetest")
```

Or, check the code on bitbucket at <https://bitbucket.org/bklamer/rankdifferencetest>.

## Examples

```{r}
# See ?rdt
library(rankdifferencetest)

res <- rdt(
  data = kornbrot_table1,
  formula = placebo ~ drug,
  alternative = "two.sided",
  distribution = "asymptotic"
)
res
```

## License

Copyright: Brett Klamer - 2021 - MIT (http://opensource.org/licenses/MIT)
