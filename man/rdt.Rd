% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdt.r
\name{rdt}
\alias{rdt}
\title{Rank difference test}
\usage{
rdt(
  data,
  formula,
  zero.method = c("Pratt", "Wilcoxon"),
  distribution = c("asymptotic", "approximate", "exact"),
  alternative = c("two.sided", "greater", "less"),
  return = c("data.frame", "coin"),
  ...
)
}
\arguments{
\item{data}{A data.frame}

\item{formula}{A formula of either form:
\describe{
\item{y ~ x | block}{
For use when \code{data} is in long format. \code{y} is the numeric outcome,
\code{x} is the binary grouping factor, and \code{block} is the
observation-level grouping factor. \code{x} and \code{block} must be factors.
}
\item{y ~ x}{
For use when \code{data} is in wide format. \code{y} is the first measurement
and \code{x} is the second measurement on the same observation.
}
}
The differences in ranks are calculated as \code{y - x} or, when \code{x} is a
factor, the first factor level minus the second.}

\item{zero.method}{A string for the method used to handle differences equal to
zero: \code{"Pratt"} (default) or \code{"Wilcoxon"}.}

\item{distribution}{A string for the method used to calculate the conditional
null distribution of the test statistic: Asymptotic distribution
\code{"asymptotic"} (default), Monte Carlo resampling \code{"approximate"},
or the exact distribution \code{"exact"}.}

\item{alternative}{A string for the alternative hypothesis: \code{"two.sided"}
(default), \code{"greater"}, or \code{"less"}.}

\item{return}{A string for the return object: \code{"data.frame"} (default) or
\code{"coin"}.}

\item{...}{Further arguments passed to \code{coin::\link[coin]{wilcoxsign_test}}.}
}
\value{
If \code{return = "data.frame"}, a data.frame with columns \code{p.value}, \code{z.statistic},
\code{formula}, \code{alternative}, and \code{method}.

If \code{return = "coin"}, an object inheriting from class \link[coin]{IndependenceTest}.
}
\description{
Performs Kornbrot's rank difference test. It's a modified Wilcoxon signed-rank
test which produces consistent and meaningful results for ordinal or monotonically
transformed data.
}
\details{
For ordinal scale data, the Wilcoxon signed-rank test results in subtraction
of those ordinal scale values. However, this subtraction is not meaningful on
the ordinal scale. In addition, any monotone transformation of the data will
result in different signed ranks, thus different p-values. However, ranking
the original data allows for meaningful addition and subtraction of ranks and
preserves ranks over monotonic transformation. Kornbrot described the rank
difference test for this reason.

Kornbrot recommends that the rank difference procedure be used in preference
to the Wilcoxon signed-rank test in all paired comparison designs where the data
are not both of interval scale type and of known distribution type. The rank
difference test preserves good power compared to Wilcoxon's signed-rank test,
is more powerful than the sign test, and has the benefit of being a true
distribution-free test.

The procedure for Wilcoxon's signed-rank test is as follows:
\enumerate{
\item Calculate differences for each paired observation.
\item Remove differences equal to zero.
\item Order the absolute differences from smallest to largest.
\item Assign ranks \eqn{1, \dots, n} with average rank for ties.
\item Calculate W+ = sum of the ranks for positive differences. The sum of W+ and
W- is \eqn{n(n+1)/2}, so either can be calculated from the other.
\item Choose the smaller of W+ and W- as the test statistic W.
\item Since the test statistic is the smaller of W+ and W-, the critical region
is the left tail of the distribution. W is distributed approximately normal
with mean \eqn{mu = (n(n+1))/4} and variance
\eqn{sigma^2 = (Tn(n+1)(2n+1))/24}, where T is a correction for ties and
\eqn{T = 1-(sum(t^3-t)/(N^3-N))}, summed over all ties, where t is the
length of a tie. The continuity corrected mean \eqn{mu = ((n(n+1))/4)+0.5}.
}

The procedure for Kornbrot's rank difference test is as follows:
\enumerate{
\item Combine all 2n observations.
\item Assign ranks \eqn{1, \dots, 2n} with average rank for ties.
\item Perform the Wilcoxon signed-rank test using the paired ranks.
}

The test statistic for the rank difference test (D) is not exactly equal to the
test statistic (W) of the naive rank-transformed Wilcoxon signed-rank test
(the latter being implemented in \code{rdt()}). Using W should result in a
conservative estimate for D, and they approach in distribution as the sample
size increases. \insertCite{kornbrot1990;textual}{rankdifferencetest}
discusses methods for calculating D when n<7 and 8<n<=20.

\code{zero.method = "Pratt"} uses the method by Pratt (1959), which first
rank-transforms the absolute differences (including zeros) and then removes
the ranks corresponding to zero-differences. \code{zero.method = "Wilcoxon"} uses
the method by Wilcoxon (1950), which first removes the zero-differences and
then rank-transforms the remaining absolute differences.
}
\examples{
#----------------------------------------------------------------------------
# rdt() example
#----------------------------------------------------------------------------
library(rankdifferencetest)

# Use example data from Kornbrot (1990)
data <- kornbrot_table1

# Create long-format data for demonstration purposes
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
# Subject and treatment must be factors. The ordering of the treatment factor
# will determine the difference (placebo - drug).
data_long$subject <- factor(data_long$subject)
data_long$treatment <- factor(data_long$treatment, levels = c("placebo", "drug"))

# Recreate analysis and results from section 7.1 in Kornbrot (1990)
## The p-value shown in Kornbrot (1990) was continuity corrected. rdt() does
## not apply a continuity correction, so the p-value here will be slightly
## lower. It does match the uncorrected p-value shown in footnote on page 246.
rdt(
  data = data,
  formula = placebo ~ drug,
  alternative = "two.sided",
  distribution = "asymptotic"
)$p.value/2
rdt(
  data = data_long,
  formula = time ~ treatment | subject,
  alternative = "two.sided",
  distribution = "asymptotic"
)$p.value/2

# The same outcome is seen after transforming time to rate.
## The rate transformation inverts the rank ordering.
data$placebo_rate <- 60 / data$placebo
data$drug_rate <- 60 / data$drug
data_long$rate <- 60 / data_long$time

rdt(
  data = data,
  formula = placebo_rate ~ drug_rate,
  alternative = "two.sided",
  distribution = "asymptotic"
)$p.value/2
rdt(
  data = data_long,
  formula = rate ~ treatment | subject,
  alternative = "two.sided",
  distribution = "asymptotic"
)$p.value/2

# In contrast to the rank difference test, the Wilcoxon signed-rank test
# produces differing results. See table 1 and table 2 in Kornbrot (1990).
wilcox.test(
  x = data$placebo,
  y = data$drug,
  paired = TRUE,
  exact = TRUE,
  alternative = "two.sided"
)$p.value/2
wilcox.test(
  x = data$placebo_rate,
  y = data$drug_rate,
  paired = TRUE,
  exact = TRUE,
  alternative = "two.sided"
)$p.value/2

}
\references{
\insertRef{kornbrot1990}{rankdifferencetest}

\insertRef{pratt1959}{rankdifferencetest}

\insertRef{wilcoxon1950}{rankdifferencetest}
}
