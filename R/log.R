#' Logarithm Function
#'
#' \code{log} returns the logarithm of the input
#'
#' This version of the logarithm function allows both numeric and complex inputs, including negative numeric values.  If the output of the
#' logarithm has no complex part then the output is given as a numeric value.
#'
#' @param x An input value (numeric/complex scalar or vector)
#' @param base The base for the logarithm (a positive scalar value)
#' @return The logarithm of the input
#'
#' @examples
#' log(1)
#' log(-1)
#'
log <- function(x, base = exp(1)) {
  LOG <- base::log(as.complex(x), base = base)
  if (all(Im(LOG) == 0)) { LOG <- Re(LOG) }
  LOG }
