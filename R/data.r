#' Alertness example data
#'
#' An example dataset as seen in table 1 from
#' \insertCite{kornbrot1990;textual}{rankdifferencetest}. The time/problem
#' was recorded for each subject under placebo and drug conditions for the
#' purpose of measuring 'alertness'.
#'
#' The table 1 values appear to be rounded, thus results do not match exactly
#' with further calculations in
#' \insertCite{kornbrot1990;textual}{rankdifferencetest}.
#'
#' @format A data frame with 13 rows and 3 variables:
#' \describe{
#'   \item{subject}{Subject identifier}
#'   \item{placebo}{The time required to complete a task under the placebo condition}
#'   \item{drug}{The time required to complete a task under the drug condition}
#' }
#'
#' @source
#' \insertRef{kornbrot1990}{rankdifferencetest}
"kornbrot_table1"
