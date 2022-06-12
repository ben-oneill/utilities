#' Logarithm Function
#'
#' \code{log} returns the logarithm of the input
#'
#' The logarithm function in base \code{R} accomodates complex numbers but it does not accomodate negative values (which is strange).  The
#' pressent version of the logarithm function allows both numeric and complex inputs, including negative numeric values.  For negative inputs
#' this function gives the princpal complex logarithm of the input value.  If the output of the logarithm has no complex part then the output
#' is given as a numeric value.  This function also allows the user to generate the gradient and Hessian.
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

  if(length(x) == 0) return(c())

  #Compute natural logarithm and adjustment (allowing for negative/complex inputs)
  LOG <- base::log(as.complex(x))
  ADJ <- 1/base::log(as.complex(base))
  LOG <- ADJ*LOG

  #Deal with special cases
  LOG[(x == 0)]    <- as.complex(-Inf)
  LOG[(x == Inf)]  <- as.complex(Inf)
  LOG[(x == -Inf)] <- complex(real = Inf, imaginary = pi)

  #Convert back to numeric if there are no complex values
  if (all(Im(LOG) == 0, na.rm=TRUE)) { LOG <- Re(LOG) }
  if (all(Im(ADJ) == 0, na.rm=TRUE)) { ADJ <- Re(ADJ) }

  #Add gradient and Hessian
  if (gradient) { attr(LOG, 'gradient') <- ADJ*(1/x) }
  if (hessian)  { attr(LOG, 'hessian')  <- -ADJ/(x^2) }

  #Return value
  LOG }

#' @rdname log
log2 <- function(x, gradient = FALSE, hessian = FALSE) {
  log(x, 2, gradient, hessian) }

#' @rdname log
log10 <- function(x, gradient = FALSE, hessian = FALSE) {
  log(x, 10, gradient, hessian) }
