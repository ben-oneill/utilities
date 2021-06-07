#' Logarithm Function
#'
#' \code{log} returns the logarithm of the input
#'
#' This version of the logarithm function allows both numeric and complex inputs, including negative numeric values.  If the output of the
#' logarithm has no complex part then the output is given as a numeric value.  It also allows the user to generate the gradient and Hessian.
#'
#' @param x An input value (numeric/complex scalar or vector)
#' @param base The base for the logarithm (a positive scalar value)
#' @param gradient Logical; if \code{TRUE} the output will include a \code{'gradient'} attribute
#' @param hessian Logical; if \code{TRUE} the output will include a \code{'hessian'} attribute
#' @return The logarithm of the input
#'
#' @examples
#' log(1)
#' log(-1)
#' log10(-10, TRUE, TRUE)

log <- function(x, base = exp(1), gradient = FALSE, hessian = FALSE) {
  LOG <- base::log(as.complex(x), base = base)
  ADJ <- 1/base::log(base)
  if (all(Im(ADJ) == 0)) { ADJ <- Re(ADJ) }
  if (gradient) { attr(LOG, 'gradient') <- ADJ*(1/x) }
  if (hessian)  { attr(LOG, 'hessian')  <- -ADJ/(x^2) }
  LOG }

#' @rdname log
log2 <- function(x, gradient = FALSE, hessian = FALSE) {
  log(x, 2, gradient, hessian) }

#' @rdname log
log10 <- function(x, gradient = FALSE, hessian = FALSE) {
  log(x, 10, gradient, hessian) }
