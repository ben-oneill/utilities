#' Softmax function
#'
#' \code{softmax} returns the value of the softmax function
#'
#' The softmax function is a bijective function that maps a real vector with length \code{m-1} to a probability vector
#' with length \code{m} with all non-zero probabilities.  The softmax function is useful in a wide range of probability
#' and statistical applications.
#'
#' @usage \code{softmax}
#' @param eta A numeric vector input
#' @param gradient Logical; if \code{TRUE} the output will include a \code{'gradient'} attribute
#' @param hessian Logical; if \code{TRUE} the output will include a \code{'hessian'} attribute
#' @return Value of the softmax function

  softmax <- function(eta, gradient = FALSE, hessian = FALSE) {

    #Compute the softmax function
    m     <- length(eta)+1
    DEN   <- matrixStats::logSumExp(c(eta, 0))
    LSOFT <- c(eta, 0) - DEN
    SOFT  <- exp(LSOFT)

    #Add the gradient (if required)
    if (gradient) {
      D1 <- array(0, dim = c(m, m-1))
      for (k in 1:m)     {
        for (i in 1:(m-1)) {
          D1[k,i] <- ((k == i) - SOFT[k])*SOFT[i] } }
      attr(SOFT, 'gradient') <- D1 }

    #Add the Hessian (if required)
    if (hessian) {
      D2 <- array(0, dim = c(m, m-1, m-1))
      for (k in 1:m)     {
        for (i in 1:(m-1)) {
          I1 <- (k == i)
          for (j in 1:(m-1)) {
            I2 <- (k == j)
            I3 <- (i == j)
            D2[k,i,j] <- I1*I2 - I3*SOFT[k] - (I1+I2)*SOFT[i]*SOFT[j] +
              2*SOFT[i]*SOFT[j]*SOFT[k] } } }
      attr(SOFT, 'hessian') <- D2 }

    #Return the output
    SOFT }

