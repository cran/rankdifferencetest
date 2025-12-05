# Standalone bkbase functions:
# rm_na
# is_formula
# fdeparse
# fmean
# fmedian
# ftabulate
# parse_formula
# tall2wide
# ftall2wide
# get_dep_data
# fouter

#' Remove missing values
#'
#' Remove missing values from an atomic vector using a pipe-compatible function.
#'
#' @param x (vector)\cr
#' The vector to remove missing values from.
#'
#' @returns
#' vector
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # rm_na() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' rm_na(1:3)
#' rm_na(c(1:3, NA))
#'
#' @noRd
rm_na <- function(x) {
  if (anyNA(x)) {
    x[!is.na(x)]
  } else {
    x
  }
}

#' Is object a formula
#'
#' Returns `TRUE` if object is a formula.
#'
#' @param x
#' The object to be tested.
#'
#' @returns
#' Scalar logical
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # is_formula() examples
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' is_formula(1:4)
#' is_formula("y ~ x")
#' is_formula(y ~ x)
#'
#' @noRd
is_formula <- function(x) {
  inherits(x, "formula")
}

#' Fast deparse
#'
#' A convenience function between [base::deparse()] and [base::deparse1()].
#'
#' @param expr
#' (expression)\cr
#' Any R expression.
#'
#' @param width.cutoff
#' (integer: 500L)\cr
#' The cutoff at which line-breaking occurs. See [base::deparse1()] if you need more.
#'
#' @param ...
#' Further arguments passed to [base::deparse()]
#'
#' @returns
#' character vector
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # fdeparse() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' fdeparse(y ~ x)
#'
#' @noRd
fdeparse <- function(expr, width.cutoff = 500L, ...) {
  deparse(expr = expr, width.cutoff = width.cutoff, ...)
}

#' Fast mean
#'
#' Calculates the arithmetic mean. A fast, unsafe alternative to [base::mean()].
#' Missing values are automatically excluded before the mean is calculated.
#'
#' @param x (numeric)\cr
#' The numeric vector used to calculate the mean.
#'
#' @returns
#' Scalar numeric
#'
#' @seealso
#' [base::mean()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # fmean() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' fmean(1:3)
#' fmean(c(1:3, NA))
#' fmean(NA)
#' fmean(numeric(0))
#'
#' @noRd
fmean <- function(x) {
  if (anyNA(x)) {
    x <- x[!is.na(x)]
  }
  sum(x) / length(x)
}

#' Fast median
#'
#' Calculates the median. A fast, unsafe alternative to [stats::median()].
#' Missing values are automatically excluded before the median is calculated.
#'
#' @param x (numeric)\cr
#' The numeric vector used to calculate the median.
#'
#' @param is_sorted (Scalar logical: `FALSE`)\cr
#' Whether or not `x` is already sorted.
#' Only set `is_sorted = TRUE` when `x` is guaranteed to already be sorted.
#'
#' @returns
#' Scalar numeric
#'
#' @seealso
#' [stats::median()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # fmedian() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' fmedian(1:3)
#' fmedian(c(1:3, NA))
#' fmedian(NA)
#' fmedian(numeric(0))
#'
#' @noRd
fmedian <- function(x, is_sorted = FALSE) {
  if (anyNA(x)) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (n == 0L) {
    return(NA_real_)
  }

  # (n + 1L) %/% 2L is slower, but returns integer
  half <- as.integer(ceiling(n / 2))

  if (is_sorted) {
    if (n %% 2L == 1L) {
      x[half]
    } else {
      fmean(x[half + 0L:1L])
    }
  } else {
    if (n %% 2L == 1L) {
      sort.int(x, partial = half)[half]
    } else {
      fmean(sort.int(x, partial = half + 0L:1L)[half + 0L:1L])
    }
  }
}

#' Fast tabulate
#'
#' Calculates the number of times each unique element occurs.
#' A fast, generalized alternative to [base::tabulate()].
#'
#' Missing values are automatically removed.
#'
#' [base::tabulate()] only counts positive integers. `ftabulate()` will count any nonmissing element type.
#'
#' @param x (Atomic vector)\cr
#' The atomic vector for tabulation.
#'
#' @param is_sorted (Scalar logical: `FALSE`)\cr
#' Whether or not `x` is already sorted. Only set `is_sorted = TRUE` when `x` is guaranteed to already be sorted.
#'
#' @returns
#' Integer vector
#'
#' @seealso
#' [base::tabulate()], [base::table()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # ftabulate() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' x <- sample(1:10, size = 100, replace = TRUE)
#' ftabulate(x)
#' tabulate(x)
#' table(x)
#'
#' @noRd
ftabulate <- function(x, is_sorted = FALSE) {
  if (anyNA(x)) {
    x <- x[!is.na(x)]
  }
  if (is_sorted) {
    # Modified from base::rle() and base::diff() v4.4.1 by R Core Team, GPL-2
    n <- length(x)
    y <- x[-1L] != x[-n]
    i <- c(0L, which(y), n)
    i[-1L] - i[-length(i)]
  } else {
    tabulate(match(x, unique(x)))
  }
}

#' Parse formula
#'
#' Parse a formula into its terms `lhs ~ rhs` where `rhs` can have two components separated by `|`.
#' For example, `y ~ x | z` has `rhs = ~ x | z` where `x` is typically described as the grouping term and `z` is the blocking term.
#' A fast, unsafe alternative to [modeltools::ParseFormula()].
#'
#' @param formula
#' (formula)\cr
#' The formula that should be parsed.
#'
#' @param specials
#' (character: `NULL`)\cr
#' A character vector of functions in `formula` marked as special in the `terms` object.
#' @param data
#' (data frame)\cr
#' Used to infer the meaning of the special symbol `.` in `formula`.
#' Can be a `data.frame` or `list`.
#'
#' @param terms
#' (scalar logical: `FALSE`)\cr
#' Whether or not to parse formula using [stats::terms.formula()].
#' Defaults to `FALSE`.
#'
#' @returns
#' list
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # parse_formula() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' parse_formula(formula("y ~ x | z"))
#' parse_formula(formula("y ~ x"))
#' parse_formula(formula("~x"))
#'
#' @importFrom stats terms.formula
#' @noRd
parse_formula <- function(
  formula,
  specials = NULL,
  data = NULL,
  terms = FALSE
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!inherits(x = formula, what = "formula")) {
    stop("Argument 'formula' must be an object of class 'formula'.")
  }
  if (!is.null(specials)) {
    if (!is.character(specials)) {
      stop("Argument 'specials' must be a character vector.")
    }
  }
  if (!is.null(data)) {
    if (!is.list(data)) {
      stop("Argument 'data' must be an object of class 'list' or 'data.frame'.")
    }
  }
  if (!is.logical(terms)) {
    stop("Argument 'terms' must be a scalar logical.")
  }

  #-----------------------------------------------------------------------------
  # Terms
  #-----------------------------------------------------------------------------
  # terms.formula expands the dot operator (y ~ .) for the data.frame context.
  # It also deduplicates formula parameters and normalizes interaction terms.
  if (!is.null(specials) || !is.null(data) || terms) {
    formula <- terms.formula(x = formula, specials = specials, data = data)
    attributes(formula) <- NULL
  }

  # length==3 if `lhs ~ rhs`
  # length==2 if `~rhs`
  if (length(formula) == 3L) {
    lhs <- formula[c(1L, 2L)]
    rhs <- formula[c(1L, 3L)]
  } else if (length(formula) == 2L) {
    lhs <- NULL
    rhs <- formula
  }

  rhs_group <- rhs
  rhs_block <- rhs

  if (!is.null(rhs) && length(rhs[[2L]]) > 1L) {
    if (fdeparse(rhs[[2]][[1]]) == "|") {
      rhs_group[[2]] <- rhs[[2]][[2]]
      rhs_block[[2]] <- rhs[[2]][[3]]
    } else {
      rhs_block <- NULL
    }
  } else {
    rhs_block <- NULL
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    formula = formula,
    lhs = lhs,
    rhs = rhs,
    rhs_group = rhs_group,
    rhs_block = rhs_block
  )
}

#' Reshape data from tall to wide
#'
#' Reshapes data from tall to wide.
#' Designed for the case when you apply the formula `y ~ x | z` on data with structure:
#' - `y` is a numeric vector for the outcome.
#' - `x` is a factor representing the grouping variable.
#' - `z` is a factor representing the blocking variable (subject/item index).
#'
#' - Rows with missing values (`NA`s) in `x` or `z` will be dropped.
#' - Rows with missing values (`NA`s) in `y` are kept.
#' - Built in aggregating functions will return `NA` if all `y` values in duplicated conditions are `NA`.
#' - Built in aggregating function all apply `na.rm = TRUE`.
#'
#' ```{r, eval=FALSE}
#' # Analogous to
#' reshape(
#'   data = data,
#'   direction = "wide",
#'   idvar = "z",
#'   timevar = "x",
#'   v.names = "y"
#' )
#' ```
#'
#' @param data
#' (data.frame)\cr
#' The tall dataframe to convert to wide format.
#'
#' @param formula
#' (formula)\cr
#' The formula to be parsed and evaluated on the data. Must be in form
#' `y ~ x | z`.
#'
#' @param y
#' (Scalar character)\cr
#' The name of the numeric outcome variable.
#'
#' @param x
#' (Scalar character)\cr
#' The name of the factor grouping variable.
#'
#' @param z
#' (Scalar character)\cr
#' The name of the factor blocking variable.
#'
#' @param agg_fun
#' (Scalar character or function: `"error"`)\cr
#' Used for aggregating duplicate cases of grouping/blocking combinations when formula of form `y ~ x | z` is used.
#' `"error"` (default) will return an error if duplicate grouping/blocking combinations are encountered.
#' Select one of `"first"`, `"last"`, `"sum"`, `"mean"`, `"median"`, `"min"`, or `"max"` for built in aggregation handling (each applies `na.rm = TRUE`).
#' Or define your own function.
#' For example, `myfun <- function(x) {as.numeric(quantile(x, 0.75, na.rm = TRUE))}`.
#'
#' @returns
#' data.frame where column order matches `levels(x)` and row order matches `levels(z)`.
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # tall2wide() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' n <- 10
#'
#' data <- data.frame(
#'   y = rnorm(n*2),
#'   x = factor(rep(letters[1:2], each = n)),
#'   z = factor(rep(seq_len(n), times = 2))
#' )
#'
#' data
#'
#' ftall2wide(data = data, y = "y", x = "x", z = "z")
#'
#' tall2wide(data = data, formula = y ~ x | z)
#'
#' reshape(
#'   data = data,
#'   direction = "wide",
#'   idvar = "z",
#'   timevar = "x",
#'   v.names = "y"
#' )
#'
#' @importFrom stats median
#' @importFrom stats reshape
#'
#' @noRd
tall2wide <- function(data, formula, agg_fun = "error") {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be an object of class 'data.frame'.")
  }
  if (!is_formula(formula)) {
    stop("Argument 'formula' must be an object of class 'formula'.")
  }

  pf <- parse_formula(formula = formula)

  lhs <- pf$lhs
  rhs_group <- pf$rhs_group
  rhs_block <- pf$rhs_block

  if (is.null(lhs) || is.null(rhs_group) || is.null(rhs_block)) {
    stop("Argument 'formula' must have form as 'y ~ x | z'.")
  }

  lhs <- fdeparse(lhs[[2]])
  rhs_group <- fdeparse(rhs_group[[2]])
  rhs_block <- fdeparse(rhs_block[[2]])

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  ftall2wide(
    data = data,
    y = lhs,
    x = rhs_group,
    z = rhs_block,
    agg_fun = agg_fun
  )
}

#' @noRd
ftall2wide <- function(data, y, x, z, agg_fun = "error") {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be an object of class 'data.frame'.")
  }
  if (!is.character(y) || length(y) != 1L) {
    stop("Argument 'y' must be a scalar character.")
  }
  if (!is.character(x) || length(x) != 1L) {
    stop("Argument 'x' must be a scalar character.")
  }
  if (!is.character(z) || length(z) != 1L) {
    stop("Argument 'z' must be a scalar character.")
  }

  if (is.character(agg_fun)) {
    if (
      !(length(agg_fun) == 1L &&
        agg_fun %in%
          c("error", "first", "last", "sum", "mean", "median", "min", "max"))
    ) {
      msg <- "Argument 'agg_fun' must be one of 'error' 'first', 'last', 'sum', 'mean', 'median', 'min', or 'max'."
      stop(msg)
    }
  } else {
    if (!is.function(agg_fun)) {
      stop("Argument 'agg_fun' must be a function or scalar character.")
    }
  }

  if (!all(c(y, x, z) %in% names(data))) {
    missing <- setdiff(c(y, x, z), names(data))
    stop(paste0(
      "Column(s) ",
      paste0(sQuote(missing, "'"), collapse = ", "),
      " not found in 'data'."
    ))
  }

  yv <- data[[y]]
  xv <- data[[x]]
  zv <- data[[z]]

  if (!is.numeric(yv)) {
    stop("Column '", y, "' must be numeric.")
  }
  if (!is.factor(xv)) {
    xv <- as.factor(xv)
  }
  if (!is.factor(zv)) {
    zv <- as.factor(zv)
  }

  # Drop rows with NA in x or z, keep y=NA as-is
  keep <- !(is.na(xv) | is.na(zv))
  if (!all(keep)) {
    yv <- yv[keep]
    xv <- droplevels(xv[keep])
    zv <- droplevels(zv[keep])
  }

  nx <- nlevels(xv)
  nz <- nlevels(zv)

  # If nothing left, return an empty wide frame with the intended structure
  if (length(yv) == 0L || nx == 0L || nz == 0L) {
    out <- as.data.frame(matrix(numeric(0), nrow = 0L, ncol = nx))
    colnames(out) <- levels(xv)
    out[[z]] <- factor(character(0), levels = levels(zv))
    return(out[c(z, levels(xv))])
  }

  # Integer indices
  zi <- as.integer(zv)
  xi <- as.integer(xv)

  # Single-key encoding
  # Key maps each (z, x) pair to 1..(nz*nx). Use row-major or column-major
  # consistently.
  if (as.double(nx) * as.double(nz) > .Machine$integer.max) {
    msg <- "`nx * nz` exceeds 32-bit integer range; consider chunking or alternative solution."
    stop(msg)
  }
  key <- xi + nx * (zi - 1L)

  # Detect duplicates
  has_dups <- anyDuplicated(key) > 0L

  # Pre-allocate output matrix (no dimnames yet)
  mat <- matrix(NA_real_, nrow = nz, ncol = nx)

  if (!has_dups) {
    mat[cbind(zi, xi)] <- yv
  }

  if (has_dups && is.character(agg_fun)) {
    # Helper to fill a set of (unique) key positions with values
    .fill_by_keys <- function(keys_u, vals) {
      if (length(keys_u) == 0L) {
        return(invisible(NULL))
      }
      rz <- ((keys_u - 1L) %/% nx) + 1L
      rx <- ((keys_u - 1L) %% nx) + 1L
      mat[cbind(rz, rx)] <<- vals
    }

    switch(
      agg_fun,
      error = {
        # Generate duplicate summary
        cnt <- tabulate(key, nbins = nz * nx)
        dup_pos <- which(cnt > 1L)
        # Map back to z/x indices
        dup_zi <- ((dup_pos - 1L) %/% nx) + 1L
        dup_xi <- ((dup_pos - 1L) %% nx) + 1L
        # Build up to 5 examples
        n_show <- min(5L, length(dup_pos))
        examples <- paste0(
          "(",
          levels(zv)[dup_zi[seq_len(n_show)]],
          ", ",
          levels(xv)[dup_xi[seq_len(n_show)]],
          ")"
        )
        stop(
          "Duplicate observations found for (",
          z,
          ", ",
          x,
          ") cells: ",
          paste(examples, collapse = ", "),
          if (length(dup_pos) > 5L) " ..." else "",
          ". Provide 'agg_fun' to combine duplicates (e.g., mean/sum/first/last/...)."
        )
      },
      first = {
        # first with na.rm=TRUE (if all NA, then return NA)
        sel <- !is.na(yv)
        if (any(sel)) {
          k2 <- key[sel]
          y2 <- yv[sel]
          idx <- !duplicated(k2)
          .fill_by_keys(k2[idx], y2[idx])
        }
      },
      last = {
        # last with na.rm=TRUE (if all NA, then return NA)
        sel <- !is.na(yv)
        if (any(sel)) {
          k2 <- key[sel]
          y2 <- yv[sel]
          idx <- !duplicated(k2, fromLast = TRUE)
          .fill_by_keys(k2[idx], y2[idx])
        }
      },
      sum = {
        # rowsum with na.rm=TRUE
        s <- rowsum(yv, key, reorder = FALSE, na.rm = TRUE)
        keys_u <- as.integer(rownames(s))
        .fill_by_keys(keys_u, as.numeric(s))
      },
      mean = {
        # mean with na.rm=TRUE (if all NA, then return NA)
        s <- rowsum(yv, key, reorder = FALSE, na.rm = TRUE)
        nn <- rowsum(
          as.integer(!is.na(yv)),
          key,
          reorder = FALSE,
          na.rm = FALSE
        )
        m <- as.numeric(s) / as.numeric(nn)
        m[!is.finite(m)] <- NA_real_ # handles nn==0 (NaN/Inf)
        keys_u <- as.integer(rownames(s))
        .fill_by_keys(keys_u, m)
      },
      median = {
        # median with na.rm=TRUE (if all NA, then return NA)
        sel <- !is.na(yv)
        if (any(sel)) {
          k2 <- key[sel]
          y2 <- yv[sel]

          # Sort by (key, value) so each group's values are contiguous and ordered
          o <- order(k2, y2)
          k2o <- k2[o]
          y2o <- y2[o]

          # Run-length encode to get group sizes and starts
          rl <- rle(k2o)$lengths
          cs <- cumsum(rl)
          starts <- cs - rl + 1L

          # Unique keys in the same order as rle (keys are sorted)
          keys_u <- k2o[starts]
          # keys_u <- unique(k2o)

          # Compute per-group medians
          odd <- (rl %% 2L) == 1L
          med <- numeric(length(rl))

          if (any(odd)) {
            # Middle element for odd-sized groups
            idx <- starts[odd] + (rl[odd] - 1L) %/% 2L
            med[odd] <- y2o[idx]
          }
          if (any(!odd)) {
            # Average of the two middle elements for even-sized groups
            idx1 <- starts[!odd] + (rl[!odd] %/% 2L - 1L)
            idx2 <- idx1 + 1L
            med[!odd] <- (y2o[idx1] + y2o[idx2]) / 2
          }

          # Fill the matrix at the (z, x) positions for these keys
          .fill_by_keys(keys_u, med)
        }
      },
      min = {
        # min with na.rm=TRUE (if all NA, then return NA)
        sel <- !is.na(yv)
        if (any(sel)) {
          k2 <- key[sel]
          y2 <- yv[sel]

          # Sort by (key, value) so the first per key is the minimum
          o <- order(k2, y2)
          k2o <- k2[o]
          y2o <- y2[o]

          # First occurrence of each key after sorting
          idx <- !duplicated(k2o)
          keys_u <- k2o[idx]
          minvals <- y2o[idx]

          # Fill the matrix at the (z, x) positions for these keys
          .fill_by_keys(keys_u, minvals)
        }
      },
      max = {
        # max with na.rm=TRUE (if all NA, then return NA)
        sel <- !is.na(yv)
        if (any(sel)) {
          k2 <- key[sel]
          y2 <- yv[sel]

          # Sort by (key, value) so the first per key is the minimum
          o <- order(k2, y2)
          k2o <- k2[o]
          y2o <- y2[o]

          # Last occurrence of each key after sorting
          idx <- !duplicated(k2o, fromLast = TRUE)
          keys_u <- k2o[idx]
          maxvals <- y2o[idx]

          # Fill the matrix at the (z, x) positions for these keys
          .fill_by_keys(keys_u, maxvals)
        }
      }
    )
  }

  # Aggregator via tapply
  if (has_dups && is.function(agg_fun)) {
    ta <- tapply(yv, list(zv, xv), FUN = agg_fun)
    ta[is.nan(ta)] <- NA_real_ # Normalize NaN to NA
    # Coerce to numeric matrix (tapply may return an array)
    mat <- suppressWarnings(matrix(as.numeric(ta), nrow = nz, ncol = nx))
  }

  # Prepare data.frame
  colnames(mat) <- levels(xv)
  out <- as.data.frame.matrix(
    mat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out[[z]] <- factor(levels(zv), levels = levels(zv))

  # Return
  out[c(z, levels(xv))]
}

#' Get dependent test data
#'
#' Get data for one-sample or dependent two-sample tests when defined by a formula of form:
#' - `y ~ x | z`
#' - `y ~ x`
#' - `~x`
#'
#' ## Dependent two-sample tests
#'
#' For dependent two-sample tests with data in tall format, the formula with form `y ~ x | z` is defined by
#'
#' - `y` is a numeric vector for the outcome.
#' - `x` is a factor with two levels representing the binary group variable.
#' - `z` is a factor with `n` levels representing the blocking variable (subject/item index).
#'
#' For factor `x`, the first level will be the reference group.
#' For example, `levels(data$x) <- c("pre", "post")` will result in the difference `post - pre`.
#'
#' For dependent two-sample tests with data in wide format, the formula with form `y ~ x` is defined by
#'
#' - `y` is a numeric vector for the group of interest.
#' - `x` is a numeric vector for the reference group.
#'
#' Differences are calculated as `data$y - data$x`.
#'
#' ## One-sample (or paired) tests
#'
#' For one-sample (or paired) tests, the formula has form `~x` where
#'
#' - `x` is a numeric vector for the outcome (or differences).
#'
#' @param data
#' (data.frame)\cr
#' The dataframe.
#'
#' @param formula
#' (formula)\cr
#' The formula to be parsed and evaluated on the data.
#'
#' @param single_term
#' (Scalar logical: `TRUE`)\cr
#' Whether or not a formula of form `~x` is allowed.
#'
#' @param agg_fun
#' (Scalar character or function: `"error"`)\cr
#' Used for aggregating duplicate cases of grouping/blocking combinations when formula of form `y ~ x | z` is used.
#' `"error"` (default) will return an error if duplicate grouping/blocking combinations are encountered.
#' Select one of `"first"`, `"last"`, `"sum"`, `"mean"`, `"median"`, `"min"`, or `"max"` for built in aggregation handling (each applies `na.rm = TRUE`).
#' Or define your own function.
#' For example, `myfun <- function(x) {as.numeric(quantile(x, 0.75, na.rm = TRUE))}`.
#'
#' @returns
#' list
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # get_dep_data() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' # Dependent two-sample tests
#' ## Long format data with 2 levels
#' data <- data.frame(
#'   y = rnorm(6),
#'   x = factor(rep(letters[1:2], each = 3)),
#'   z = factor(rep(1:3, times = 2))
#' )
#' formula <- y ~ x | z
#' get_dep_data(data = data, formula = formula)
#'
#' ## Wide format data
#' data <- data.frame(
#'   y = rnorm(3),
#'   x = rnorm(3)
#' )
#' formula <- y ~ x
#' get_dep_data(data = data, formula = formula)
#'
#' # Formula for single term
#' data <- data.frame(
#'   x = rnorm(3)
#' )
#' formula <- ~x
#' get_dep_data(data = data, formula = formula)
#'
#' @noRd
get_dep_data <- function(
  data,
  formula,
  single_term = TRUE,
  agg_fun = "error"
) {
  #-----------------------------------------------------------------------------
  # Check arguments
  #-----------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be an object of class 'data.frame'.")
  }
  if (!is_formula(formula)) {
    stop("Argument 'formula' must be an object of class 'formula'.")
  }
  if (!(is.logical(single_term) && length(single_term) == 1L)) {
    stop("Argument 'single_term' must be a scalar logical.")
  }

  if (is.character(agg_fun)) {
    if (
      !(length(agg_fun) == 1L &&
        agg_fun %in%
          c("error", "first", "last", "sum", "mean", "median", "min", "max"))
    ) {
      msg <- "Argument 'agg_fun' must be one of 'error', 'first', 'last', 'sum', 'mean', 'median', 'min', or 'max'."
      stop(msg)
    }
  } else {
    if (!is.function(agg_fun)) {
      stop("Argument 'agg_fun' must be a function or scalar character.")
    }
  }

  data <- as.list.data.frame(data)

  #-----------------------------------------------------------------------------
  # Parse formula
  #-----------------------------------------------------------------------------
  # For example:
  # rhs_group is a formula
  # rhs_group_char is a character vector with possible multiple elements for I()
  # deparse(rhs_group) is a string for the concatenated character vector
  pf <- parse_formula(formula = formula)

  lhs <- pf$lhs
  rhs_group <- pf$rhs_group
  rhs_block <- pf$rhs_block

  lhs_char <- if (is.null(lhs)) NULL else as.character(lhs[[2L]])
  rhs_group_char <- if (is.null(rhs_group)) {
    NULL
  } else {
    as.character(rhs_group[[2L]])
  }
  rhs_block_char <- if (is.null(rhs_block)) {
    NULL
  } else {
    as.character(rhs_block[[2L]])
  }

  #-----------------------------------------------------------------------------
  # Validate formula
  #-----------------------------------------------------------------------------
  forms <- if (single_term) {
    "'y ~ x | z', 'y ~ x', or '~x'"
  } else {
    "'y ~ x | z' or 'y ~ x'"
  }

  if (is.null(rhs_group)) {
    stop(
      paste0(
        "The formula is missing the group variable on the right hand side. For example, it must be of form ",
        forms,
        "."
      )
    )
  }
  if (is.null(lhs) && !single_term) {
    stop(
      "The formula is missing the left hand side. For example, it must be of form ",
      forms,
      "."
    )
  }

  # Check for malformed formula and insert new data created by `I()`
  # To avoid issues applying get_dep_data() twice, rename to `.I()` and ignore if found.
  if (length(lhs_char) > 1L) {
    if (lhs_char[1] == "I") {
      lhs_char <- paste0(".", fdeparse(lhs[[2]]))
      data[[lhs_char]] <- as.numeric(eval(lhs[[2]], data))
    } else if (lhs_char[1] == ".I") {
      # Nothing
    } else {
      stop(
        "The formula must not have multiple left hand side components. For example, it must be of form ",
        forms,
        "."
      )
    }
  }

  if (length(rhs_group_char) > 1L) {
    if (rhs_group_char[1] == "I") {
      rhs_group_char <- paste0(".", fdeparse(rhs_group[[2]]))
      data[[rhs_group_char]] <- as.numeric(eval(rhs_group[[2]], data))
    } else if (rhs_group_char[1] == ".I") {
      # Nothing
    } else {
      stop(
        "The formula must not have multiple 'group' components. For example, it must be of form ",
        forms,
        "."
      )
    }
  }

  if (length(rhs_block_char) > 1L) {
    if (rhs_block_char[1] == "I") {
      rhs_block_char <- paste0(".", fdeparse(rhs_block[[2]]))
      data[[rhs_block_char]] <- as.numeric(eval(rhs_block[[2]], data))
    } else if (rhs_block_char[1] == ".I") {
      # Nothing
    } else {
      stop(
        "The formula must not have multiple 'block' components. For example, it must be of form ",
        forms,
        "."
      )
    }
  }

  # Terms must be in data
  vars <- c(lhs_char, rhs_group_char, rhs_block_char)
  if (!all(vars %in% names(data))) {
    missing_vars <- vars[!vars %in% names(data)]
    stop(paste0(
      "The formula included terms not found in the data: ",
      paste0(sQuote(missing_vars, "'"), collapse = ", ")
    ))
  }

  # lhs must be numeric
  if (!is.null(lhs)) {
    if (!is.numeric(data[[lhs_char]])) {
      stop(paste0(
        sQuote(lhs_char, "'"),
        " must be numeric in formula '",
        fdeparse(pf$formula),
        "'."
      ))
    }
  }

  # rhs should be 2-level factor or numeric
  if (!is.null(rhs_block)) {
    # Group
    if (!is.factor(data[[rhs_group_char]])) {
      data[[rhs_group_char]] <- as.factor(data[[rhs_group_char]])
      warning(paste0(
        sQuote(rhs_group_char, "'"),
        " was converted to a factor in formula '",
        fdeparse(pf$formula),
        "'."
      ))
    }
    if (nlevels(data[[rhs_group_char]]) != 2L) {
      stop(paste0(
        sQuote(rhs_group_char, "'"),
        " must be a two-level factor in formula '",
        fdeparse(pf$formula),
        "'."
      ))
    }

    # Block
    if (!is.factor(data[[rhs_block_char]])) {
      data[[rhs_block_char]] <- as.factor(data[[rhs_block_char]])
      warning(paste0(
        sQuote(rhs_block_char, "'"),
        " was converted to a factor in formula '",
        fdeparse(pf$formula),
        "'."
      ))
    }
  } else {
    if (!is.numeric(data[[rhs_group_char]])) {
      stop(paste0(
        sQuote(rhs_group_char, "'"),
        " must be numeric in formula '",
        fdeparse(pf$formula),
        "'."
      ))
    }
  }

  #-----------------------------------------------------------------------------
  # Prepare data
  #-----------------------------------------------------------------------------
  if (is.null(rhs_block)) {
    if (is.null(lhs)) {
      ret <- data[rhs_group_char]
      focal_group <- rhs_group_char
      ref_group <- NULL
    } else {
      ret <- data[c(lhs_char, rhs_group_char)]
      focal_group <- lhs_char
      ref_group <- rhs_group_char
    }
  } else {
    ret <- data[c(lhs_char, rhs_group_char, rhs_block_char)]
    attr(ret, "names") <- names(ret)
    attr(ret, "row.names") <- .set_row_names(length(ret[[1L]]))
    attr(ret, "class") <- "data.frame"

    ret <- ftall2wide(
      data = ret,
      y = lhs_char,
      x = rhs_group_char,
      z = rhs_block_char,
      agg_fun = agg_fun
    )

    lvls <- levels(data[[rhs_group_char]])
    ref_group <- lvls[1L]
    focal_group <- lvls[2L]
    ret <- ret[c(focal_group, ref_group)]
  }

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------
  list(
    data = ret,
    focal_group = focal_group,
    ref_group = ref_group,
    parse_formula = pf
  )
}

#' Fast outer
#'
#' A fast, unsafe alternative to [base::outer()].
#'
#' @param x
#' (atomic vector)\cr
#' The first argument for function `fun`.
#'
#' @param y
#' (atomic vector)\cr
#' The second argument for function `fun`.
#'
#' @param fun
#' (function)\cr
#' The function to apply to `x` and `y`.
#'
#' @returns
#' matrix
#'
#' @seealso
#' [base::outer()]
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # fouter() example
#' #----------------------------------------------------------------------------
#' library(bkbase)
#'
#' x <- c(1, 2)
#' y <- c(3, 4)
#'
#' fouter(x, y, `+`)
#'
#' @noRd
fouter <- function(x, y, fun) {
  nx <- length(x)
  ny <- length(y)
  y <- rep(y, rep.int(nx, ny))
  x <- rep(x, times = ceiling(length(y) / nx))
  res <- fun(x, y)
  dim(res) <- c(nx, ny)
  res
}
