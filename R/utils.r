# Extract data based on the model formula.
#' @importFrom modeltools ModelEnvFormula has
#' @importFrom stats na.omit
model2data <- function(formula, data) {
  # ModelEnvFormula returns the possibly transformed data and includes the
  # transformation in the variable name.
  md <- ModelEnvFormula(
    formula = formula,
    data = data,
    na.action = na.omit,
    designMatrix = FALSE,
    responseMatrix = FALSE
  )

  if(has(md, "input")) {
    x <- md@get("input")
  } else {
    stop("Check the formula used in rdt(). It's missing the right hand side!")
  }

  if(has(md, "response")) {
    y <- md@get("response")
  } else {
    stop("Check the formula used in rdt(). It's missing the left hand side!")
  }

  if(has(md, "blocks")) {
    block <- md@get("blocks")
  } else {
    block <- NULL
  }

  list(
    y = y,
    x = x,
    block = block
  )
}
