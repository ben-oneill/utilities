#' Softmax and inverse-softmax functions
#'
#' \code{softmax} returns the value of the softmax function.
#' \code{softmaxinv} returns the value of the inverse-softmax function.
#'
#' The softmax function is a bijective function that maps a real vector with length \code{m-1} to a probability vector
#' with length \code{m} with all non-zero probabilities.  The softmax function is useful in a wide range of probability
#' and statistical applications.  The present functions define the softmax function and its inverse, both with a tuning
#' parameter.
#'
#' @param eta A numeric vector input
#' @param lambda Tuning parameter (a single positive value)
#' @param gradient Logical; if \code{TRUE} the output will include a \code{'gradient'} attribute
#' @param hessian Logical; if \code{TRUE} the output will include a \code{'hessian'} attribute
#' @return Value of the softmax function
#'
#' @examples
#' softmax(5:7)
#' softmaxinv(softmax(5:7))

softmax <- function(eta, lambda = 1, gradient = FALSE, hessian = FALSE) {

  stopifnot(requireNamespace("matrixStats", quietly = TRUE))

  #Compute the softmax function
  m     <- length(eta)+1
  DEN   <- matrixStats::logSumExp(c(lambda*eta, 0))
  LSOFT <- c(lambda*eta, 0) - DEN
  SOFT  <- exp(LSOFT)

  #Add the gradient (if required)
  if (gradient) {
    D1 <- array(0, dim = c(m, m-1))
    for (k in 1:m)     {
      for (i in 1:(m-1)) {
        D1[k,i] <- ((k == i) - SOFT[k])*SOFT[i] } }
    attr(SOFT, 'gradient') <- lambda*D1 }

  #Add the Hessian (if required)
  if (hessian) {
    D2 <- array(0, dim = c(m, m-1, m-1))
    for (k in 1:m)     {
      for (i in 1:(m-1)) {
        I1 <- (k == i)
        for (j in 1:(m-1)) {
          I2 <- (k == j)
          I3 <- (i == j)
          D2[k,i,j] <- SOFT[k]*(I1*I2 + 2*SOFT[i]*SOFT[j] - I1*SOFT[j] - I2*SOFT[i] - I3*SOFT[i]) } } }
    attr(SOFT, 'hessian') <- (lambda^2)*D2 }

  #Return the output
  SOFT }

#' @param p A probability vector (i.e., numeric vector of non-negative values that sum to one)
#' @rdname softmax
softmaxinv <- function(p, lambda = 1, gradient = FALSE, hessian = FALSE) {

  #Compute the inverse-softmax function
  m <- length(p)
  if (m > 1) {
    SOFTINV <- (base::log(p) - base::log(p[m]))[1:(m-1)]/lambda } else {
    SOFTINV <- numeric(0) }

  #Add the gradient (if required)
  if (gradient) {
    D1 <- array(0, dim = c(m-1, m))
    if (m > 1) {
    for (k in 1:(m-1))   {
      for (i in 1:(m-1)) {
        D1[k,i] <-  1/p[i] }
        D1[k,m] <- -1/p[m] } }
    attr(SOFTINV, 'gradient') <- D1/lambda }

  #Add the Hessian (if required)
  if (hessian) {
    D2 <- array(0, dim = c(m-1, m, m))
    if (m > 1) {
    for (k in 1:(m-1))   {
      for (i in 1:(m-1)) {
        D2[k, i, i] <- -1/p[i]^2 }
        D2[k, m, m] <-  1/p[m]^2 } }
    attr(SOFTINV, 'hessian') <- D2/(lambda^2) }

  #Return the output
  SOFTINV }

