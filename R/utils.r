#' Walsh averages
#'
#' Calculates the set of Walsh (all pairwise) averages.
#' Missing values are silently excluded.
#'
#' The set of Walsh (all pairwise) averages is defined as
#'
#' \deqn{
#' w_{ij} = \frac{x_i + x_j}{2}, \quad \text{for } 1 \leq i \leq j \leq n.
#' }
#'
#' with length \eqn{\frac{n(n+1)}{2}}.
#'
#' The Hodges-Lehman pseudomedian is computed as the median of the Walsh averages.
#'
#' @param x
#' (numeric)\cr
#' The numeric vector used to calculate the Walsh averages.
#' Missing values are silently excluded.
#'
#' @param sort
#' (Scalar logical: `FALSE`)\cr
#' Whether or not the vector of walsh averages should be sorted.
#'
#' @returns
#' numeric
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # walsh() example
#' #----------------------------------------------------------------------------
#'
#' walsh(1:3)
#' walsh(c(1:3, NA))
#'
#' fmedian(walsh(1:3))
#'
#' @noRd
walsh <- function(x, sort = FALSE) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.numeric(x)) {
    stop("Argument 'x' must be a numeric vector.")
  }
  if (!is.logical(sort) || length(sort) != 1L) {
    stop("Argument 'sort' must be a scalar logical.")
  }

  #-----------------------------------------------------------------------------
  # Calculate walsh averages
  #-----------------------------------------------------------------------------
  if (length(x) == 0L) {
    return(NA_real_)
  }
  if (anyNA(x)) {
    x <- x[!is.na(x)]
  }

  walsh <- fouter(x, x, `+`)
  walsh <- walsh[lower.tri(walsh, diag = TRUE)] / 2
  if (sort) {
    walsh <- sort(walsh)
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  walsh
}
