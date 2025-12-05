#' @title
#' Get data for the signed-rank test.
#'
#' @description
#' Gets data needed to compute signed-rank results.
#' The returned list is designed to be reused by higher-level functions.
#'
#' @inheritParams srt
#'
#' @param x
#' (numeric)\cr
#' Numeric vector or data.frame of data.
#' Values with non-finite values (infinite or missing) are silently dropped.
#'
#' @param x_name
#' (Scalar character)\cr
#' Name of `x` variable.
#'
#' @param y_name
#' (Scalar character or `NULL`)\cr
#' Name of `y` variable.
#' If `y = NULL` then `y_name = NULL`.
#'
#' @param ...
#' Unused additional arguments.
#'
#' @returns
#' list
#'
#' @keywords internal
#' @export
srt_data <- function(x, ...) {
  UseMethod("srt_data")
}

#' @rdname srt_data
#' @keywords internal
#' @export
srt_data.numeric <- function(x, y = NULL, x_name, y_name = NULL, ...) {
  if (!is.null(y)) {
    if (!is.numeric(y)) {
      stop("Argument 'y' must be an object of class 'numeric'.")
    }
    if (length(x) != length(y)) {
      stop("Arguments 'x' and 'y' must be the same length.")
    }
    n_sample <- length(x)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok] - y[ok]
    y <- NULL
  } else {
    n_sample <- length(x)
    x <- x[is.finite(x)]
  }

  n_analytic <- length(x)

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    method = "Wilcoxon",
    diffs_original = x,
    diffs = x,
    focal_name = x_name,
    reference_name = y_name,
    n_sample = n_sample,
    n_analytic = n_analytic
  )
}

#' @rdname srt_data
#' @keywords internal
#' @export
srt_data.data.frame <- function(x, formula, agg_fun = "error", ...) {
  #-----------------------------------------------------------------------------
  # formula to data
  #-----------------------------------------------------------------------------
  lst <- get_dep_data(
    data = x,
    formula = formula,
    single_term = TRUE,
    agg_fun = agg_fun
  )
  data <- lst$data

  #-----------------------------------------------------------------------------
  # prepare data
  #-----------------------------------------------------------------------------
  diffs <- if (length(data) == 2L) {
    data[[1L]] - data[[2L]]
  } else if (length(data) == 1L) {
    data[[1L]]
  } else {
    stop(
      paste0(
        "Check argument 'data' and/or 'formula'.\n",
        "The data preparation step created a ",
        ncol(data),
        " column data frame.\n",
        "However, it must only have 1 or 2 columns."
      )
    )
  }

  n_sample <- length(diffs)
  diffs <- diffs[is.finite(diffs)]
  n_analytic <- length(diffs)

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    method = "Wilcoxon",
    diffs_original = diffs,
    diffs = diffs,
    focal_name = lst$focal_group,
    reference_name = lst$ref_group,
    n_sample = n_sample,
    n_analytic = n_analytic
  )
}
