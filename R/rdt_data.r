#' @title
#' Get data for the Kornbrot rank difference test.
#'
#' @description
#' Gets data needed to compute Kornbrot rank difference results.
#' The returned list is designed to be reused by higher-level functions.
#'
#' @inheritParams rdt
#'
#' @param x
#' (numeric)\cr
#' Numeric vector or data.frame of data.
#' If numeric, differences of ranks correspond with `x - y`.
#' Pairs with missing values are silently dropped.
#'
#' @param x_name
#' (Scalar character)\cr
#' Name of `x` variable.
#'
#' @param y_name
#' (Scalar character)\cr
#' Name of `y` variable.
#'
#' @param ...
#' Unused additional arguments.
#'
#' @returns
#' list
#'
#' @importFrom stats complete.cases
#' @rdname rdt_data
#' @keywords internal
#' @export
rdt_data <- function(x, ...) {
  UseMethod("rdt_data")
}

#' @rdname rdt_data
#' @keywords internal
#' @export
rdt_data.numeric <- function(x, y, x_name, y_name, ...) {
  if (!is.numeric(x)) {
    stop("Argument 'x' must be an object of class 'numeric'.")
  }
  if (!is.numeric(y)) {
    stop("Argument 'y' must be an object of class 'numeric'.")
  }
  if (length(x) != length(y)) {
    stop("Arguments 'x' and 'y' must be the same length.")
  }

  n_sample <- length(x)

  # ranks allow us to keep infinites
  if (anyNA(x) || anyNA(y)) {
    ok <- complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]
  }

  n_analytic <- length(y)

  # Note: Kornbrot orders opposite direction: rank(-xtfrm(data)).
  ranks <- rank(c(x, y), na.last = NA)
  diffs <- ranks[seq_len(n_analytic)] - ranks[-seq_len(n_analytic)]

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    method = "Kornbrot-Wilcoxon",
    diffs_original = diffs,
    diffs = diffs,
    focal_name = x_name,
    reference_name = y_name,
    n_sample = n_sample,
    n_analytic = n_analytic
  )
}

#' @rdname rdt_data
#' @keywords internal
#' @export
rdt_data.data.frame <- function(x, formula, agg_fun = "error", ...) {
  #-----------------------------------------------------------------------------
  # formula to data
  #-----------------------------------------------------------------------------
  lst <- get_dep_data(
    data = x,
    formula = formula,
    single_term = FALSE,
    agg_fun = agg_fun
  )
  data <- lst$data

  #-----------------------------------------------------------------------------
  # prepare return
  #-----------------------------------------------------------------------------
  if (length(data) != 2L) {
    msg <- paste0(
      "Check argument 'data' and/or 'formula'.\n",
      "The data preparation step created a ",
      ncol(data),
      " column data frame.\n",
      "However, it must have exactly 2 columns."
    )
    stop(msg)
  }

  # y is focal; x is reference
  y <- data[[1L]]
  x <- data[[2L]]

  # Argument check
  if (!(is.numeric(y) && is.numeric(x))) {
    stop(
      "Data must be numeric."
    )
  }

  n_sample <- length(y)

  # ranks allow us to keep infinite
  if (anyNA(y) || anyNA(x)) {
    ok <- complete.cases(y, x)
    y <- y[ok]
    x <- x[ok]
  }

  n_analytic <- length(y)

  # Note: Kornbrot orders opposite direction: rank(-xtfrm(data)).
  ranks <- rank(c(y, x), na.last = NA)
  diffs <- ranks[seq_len(n_analytic)] - ranks[-seq_len(n_analytic)]

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    method = "Kornbrot-Wilcoxon",
    diffs_original = diffs,
    diffs = diffs,
    focal_name = lst$focal_group,
    reference_name = lst$ref_group,
    n_sample = n_sample,
    n_analytic = n_analytic
  )
}
