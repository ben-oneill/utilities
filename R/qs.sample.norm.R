#' Quantile-stratified sampling from a normal distribution
#'
#' \code{qs.sample.norm} returns a quantile-stratified sample from a normal distribution.
#'
#' This function generates a quantile-stratified sample of the stipulated size from the normal distribution (univariate
#' or multivariate).   Quantile-stratified sampling is a method of sampling from a known probability distribution
#' that is analogous to sampling-without-replacement; it induces negative correlation between the sample values and the
#' empirical quantiles adhere more closely to the true quantile function than for an IID sample.  Further details on
#' this sampling method and its properties can be found in the following paper:
#' 
#' O'Neill, B. (2025) Quantile-stratified sampling for multivariate normal distributions and other multivariate
#' distributions.
#' 
#' @usage \code{qs.sample.norm(n, mean, var, layers = NULL)}
#' @param n The number of sample values to be generated (a non-negative integer)
#' @param mean The mean vector for the normal distribution
#' @param var The variance matrix for the normal distribution
#' @param layers Optional vector giving the number of sample values in each layer of the sample
#' @return A quantile-stratified sample of size n generated from the normal distribution with the specified mean 
#' and variance

qs.sample.norm <- function(n, mean, var, layers = NULL) {
  
  #Check input n
  if (!is.vector(n))                         stop('Error: Input n should be a vector')
  if (!is.numeric(n))                        stop('Error: Input n should be a numeric value')
  if (length(n) != 1)                        stop('Error: Input n should be a single numeric value')
  nn <- as.integer(n)
  if (nn != n)                               stop('Error: Input n should be an integer')
  if (nn < 0)                                stop('Error: Input n should be a non-negative integer')
  
  #Check input mean
  if (is.matrix(mean)) {
    if (ncol(mean) != 1)                     stop('Error: Input mean should be a vector or a matrix with a single column') }
  if ((!is.vector(mean))&(!is.matrix(mean))) stop('Error: Input mean should be a vector or a matrix with a single column')
  if (!is.numeric(mean))                     stop('Error: Input mean should be numeric')
  MEAN <- matrix(mean, ncol = 1)
  DIM  <- length(mean)
  
  #Check input var
  if (is.matrix(var)) {
    if (nrow(var) != DIM)                    stop('Error: Input var should have the same number of rows as length of mean')
    if (ncol(var) != DIM)                    stop('Error: Input var should have the same number of columns as length of mean') 
    if (!is.numeric(mean))                   stop('Error: Input var should be a numeric matrix')
    if (any(var != t(var)))                  stop('Error: Input var should be a symmetric matrix')
    VAR <- var }
  if ((!is.matrix(var))&(is.vector(var))) {
    if (DIM != 1)                            stop('Error: Input var should be a matrix with rows and columns equal to the length of mean')
    if (length(var) != 1)                    stop('Error: Input var should be a matrix')
    if (!is.numeric(mean))                   stop('Error: Input var should be a numeric vector or matrix')
    VAR <- as.matrix(var) }
  if ((!is.matrix(var))&(!is.vector(var)))   stop('Error: Input var must be a vector or matrix')
  VAR.EIGEN <- eigen(VAR)
  if (min(VAR.EIGEN$values < 0))             stop('Error: Input var must be non-negative definite')
  
  #Check input layers
  if (is.null(layers)) {
    LAYERS <- nn 
  } else {
    if (!is.vector(layers))                  stop('Error: Input layers should be a vector')
    if (!is.numeric(layers))                 stop('Error: Input layers should be a numeric vector')
    LAYERS <- as.integer(layers)
    if (any(LAYERS != layers))                   stop('Error: Input layers should contain integers')
    if (min(LAYERS) < 0)                         stop('Error: Input layers should contain positive integers')
    if (sum(LAYERS) != nn)                       stop('Error: Input layers should contain values that add up to input n')
    LAYERS <- LAYERS[LAYERS > 0] }
  
  #Deal with special case where n = 0
  if (nn == 0) return(numeric(0))
  
  #Generate standard deviation matrix
  A   <- diag(VAR.EIGEN$values)
  V   <- VAR.EIGEN$vectors
  STD <- V %*% sqrt(A) %*% t(V)
  
  #Generate QS sample
  OUT <- matrix(0, nrow = DIM, ncol = n)
  rownames(OUT) <- sprintf('x[%s]', 1:DIM)
  colnames(OUT) <- sprintf('SIM[%s]', 1:n)
  RR <- sqrt(qs.sample(n, Q = qchisq, layers = LAYERS, df = DIM))
  ZZ <- matrix(rnorm(n*DIM), nrow = DIM, ncol = n)
  for (i in 1:n) {
    OUT[,i] <- MEAN + RR[i]*(STD %*% ZZ[,i])/sqrt(sum(ZZ[,i]^2)) }
  
  #Give output
  OUT }
