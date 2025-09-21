#' Quasi-quantile-stratified sampling (from a known distribution)
#'
#' \code{qs.sample.quasi} returns a quasi-quantile-stratified sample from a specified quantile function.
#'
#' This function generates a quasi-quantile-stratified sample of the stipulated size using the input random generator and 
#' density function for the sampling distribution.  Quantile-stratified sampling is a method of sampling from a known 
#' probability distribution that is analogous to sampling-without-replacement; it induces negative correlation between 
#' the sample values and the empirical quantiles adhere more closely to the true quantile function than for an IID sample.  
#' Quasi-quantile-stratified sampling approximates this method by generating independent samples and using their order with
#' respect to their density to estimate the quantiles to perform quantile-stratified sampling.  Further details on this 
#' sampling method and its properties can be found in the following papers:
#' 
#' O'Neill, B. (2025) One-dimensional quantile-stratified sampling and its applications in statistical simulations.
#'  
#' The user must specify the random generation function \code{r} and density function \code{d} for the sampling distribution and 
#' and the \code{multiplier} used for generating the quasi-quantile-stratified sample.  The user can also simulate using layered 
#' quasi-quantile-stratified sampling by giving a vector of sizes for each of the layers (which must add up to the total sample size).
#' 
#' WARNING: Quasi-quantile-stratified sampling assumes a continuous quantile function and the sampling method is generally 
#' inappropriate/incorrect if the quantile function is non-continuous.  The present function will not be able to detect 
#' whether or not the input quantile function is for a continuous distribution or not so the user should be aware of this.  
#' 
#' @usage \code{qs.sample.quasi(n, multiplier, r, d, size.arg = 'n', input.arg = 'x', layers = NULL, ...)}
#' @param n The number of sample values to be generated (a non-negative integer)
#' @param multiplier The multiplier used for quasi-quantile-stratified sampling (must be a positive integer)
#' @param r The random generation function for the sampling distribution (must be a function)
#' @param d The density function for the sampling distribution (must be a function)
#' @param size.arg The name of the input sample size argument in the random generation function r
#' @param input.arg The name of the input argument in the density function d
#' @param layers Optional vector giving the number of sample values in each layer of the sample
#' @param ... Distribution parameters to be passed through to the quantile function Q
#' @return A quantile-stratified sample of size n generated from the quantile function Q

qs.sample.quasi <- function(n, multiplier, r, d, size.arg = 'n', input.arg = 'x', layers = NULL, ...) {
  
  #Check input n
  if (!is.vector(n))                         stop('Error: Input n should be a vector')
  if (!is.numeric(n))                        stop('Error: Input n should be a numeric value')
  if (length(n) != 1)                        stop('Error: Input n should be a single numeric value')
  nn <- as.integer(n)
  if (nn != n)                               stop('Error: Input n should be an integer')
  if (nn < 0)                                stop('Error: Input n should be a non-negative integer')
  
  #Check input multiplier
  if (!is.vector(n))                         stop('Error: Input multiplier should be a vector')
  if (!is.numeric(n))                        stop('Error: Input multiplier should be a numeric value')
  if (length(n) != 1)                        stop('Error: Input multiplier should be a single numeric value')
  nn <- as.integer(n)
  if (nn != n)                               stop('Error: Input multiplier should be an integer')
  if (nn <= 0)                               stop('Error: Input multiplier should be a positive integer')
  
  #Check input input.arg and size.arg
  if (!is.vector(input.arg))                 stop('Error: Input input.arg should be a vector')
  if (!is.character(input.arg))              stop('Error: Input input.arg should be a character string')
  if (length(input.arg) != 1)                stop('Error: Input input.arg should be a single character string')
  if (!is.vector(size.arg))                  stop('Error: Input size.arg should be a vector')
  if (!is.character(size.arg))               stop('Error: Input size.arg should be a character string')
  if (length(size.arg) != 1)                 stop('Error: Input size.arg should be a single character string')
  
  #Check inputs r and d
  if (!is.function(r))                       stop('Error: Input r should be a function')
  if (!(size.arg %in% names(formals(r))))    stop('Error: Input r should have an argument given by size.arg')
  if (!is.function(d))                       stop('Error: Input d should be a function')
  if (!(input.arg %in% names(formals(d))))   stop('Error: Input d should have an argument given by input.arg')
  
  #Check input layers
  if (is.null(layers)) {
    LAYERS <- nn 
  } else {
    if (!is.vector(layers))                  stop('Error: Input layers should be a vector')
    if (!is.numeric(layers))                 stop('Error: Input layers should be a numeric vector')
    LAYERS <- as.integer(layers)
    if (any(LAYERS != layers))               stop('Error: Input layers should contain integers')
    if (min(LAYERS) < 0)                     stop('Error: Input layers should contain positive integers')
    if (sum(LAYERS) != nn)                   stop('Error: Input layers should contain values that add up to input n')
    LAYERS <- LAYERS[LAYERS > 0] }
  
  #Deal with special case where n = 0
  if (nn == 0) return(numeric(0))
  
  #Create quasi-QS samples for each layer
  M <- length(LAYERS)
  SAMPLES <- vector(mode = 'list', length = M)
  for (m in 1:M) {
    
    #Generate quasi-quantile-stratified sample from r using d for ordering (using distribution parameters in ... input)
    #Method generates (nn x multiplier) random values from r and orders them by the log-density d(...,log = TRUE).  
    #Method samples one value in each consecutive set of values with size = multiplier, yielding nn sampled values.
    #The vector RRR contains the randomly generated values from the sample distribution, ordered by increasing density
    #The vector assigned to SAMPLES[[m]] is the quasi-quantile sample for this layer
    nnn  <- LAYERS[m]
    ARGSr <- list(...)
    ARGSr[[size.arg]]  <- nn*multiplier
    ARGSr$log.p        <- NULL
    RR <- do.call(r, args = ARGSr)
    ARGSd <- list(...)
    ARGSd[[input.arg]] <- RR
    ARGSd$log          <- TRUE
    DD <- do.call(d, args = ARGSd)
    RRR <- RR[order(DD)]
    
    #Generate the quasi-quantile-stratified sample for this layer
    UNIF <- floor(multiplier*runif(nn))
    SAMPLES[[m]] <- RRR[(1:nn)*multiplier - UNIF] }
  
  #Combine and randomly permute sample layers
  PERM <- sample(1:nn, replace = FALSE)
  OUT <- unlist(SAMPLES)[PERM]
  
  #Give output
  OUT }
