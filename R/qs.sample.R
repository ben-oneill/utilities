#' Quantile-stratified sampling (from a known distribution)
#'
#' \code{qs.sample} returns a quantile-stratified sample from a specified quantile function.
#'
#' This function generates a quantile-stratified sample of the stipulated size using the input quantile function for
#' the sampling distribution.  Quantile-stratified sampling is a method of sampling from a known probability distribution
#' that is analogous to sampling-without-replacement; it induces negative correlation between the sample values and the
#' empirical quantiles adhere more closely to the true quantile function than for an IID sample.  Further details on
#' this sampling method and its properties can be found in the following paper:
#' 
#' O'Neill, B. (2025) Quantile-stratified sampling and quantile-stratified importance sampling.
#' 
#' The user must specify the quantile function Q for the sampling distribution and the argument in this quantile function
#' that represents the probability input.  For most quantile functions programmed in \code{R} the probability input has the 
#' argument name \code{'p'} which is used as the default.  (Other argument names can be specified if the quantile function uses
#' a different variable name.)  The user can also simulate using layered quantile-stratified sampling by giving a vector
#' of sizes for each of the layers (which must add up to the total sample size).
#' 
#' WARNING: Quantile-stratified sampling assumes a continuous quantile function  and the sampling method is generally 
#' inappropriate/incorrect if the quantile function is non-continuous.  The present function will not be able to detect 
#' whether or not the input quantile function is for a continuous distribution or not so the user should be aware of this.  
#' 
#' @usage \code{qs.sample(n, Q, prob.arg = 'p', layers = NULL, ...)}
#' @param n The number of sample values to be generated (a non-negative integer)
#' @param Q The quantile function of the sampling distribution (must be a function)
#' @param prob.arg The name of the probability argument in the quantile function Q
#' @param layers Optional vector giving the number of sample values in each layer of the sample
#' @param ... Distribution parameters to be passed through to the quantile function Q
#' @return A quantile-stratified sample of size n generated from the quantile function Q

qs.sample <- function(n, Q, prob.arg = 'p', layers = NULL, ...) {
  
  #Check input n
  if (!is.vector(n))                         stop('Error: Input n should be a vector')
  if (!is.numeric(n))                        stop('Error: Input n should be a numeric value')
  if (length(n) != 1)                        stop('Error: Input n should be a single numeric value')
  nn <- as.integer(n)
  if (nn != n)                               stop('Error: Input n should be an integer')
  if (nn < 0)                                stop('Error: Input n should be a non-negative integer')
  
  #Check input prob.arg
  if (!is.vector(prob.arg))                  stop('Error: Input prob.arg should be a vector')
  if (!is.character(prob.arg))               stop('Error: Input prob.arg should be a character string')
  if (length(prob.arg) != 1)                 stop('Error: Input prob.arg should be a single character string')
  
  #Check input Q
  if (!is.function(Q))                       stop('Error: Input Q should be a function')
  if (!(prob.arg %in% names(formals(Q))))    stop('Error: Input Q should have an argument given by prob.arg')
  
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
  
  #Create QS samples for each layer
  M <- length(LAYERS)
  SAMPLES <- vector(mode = 'list', length = M)
  for (m in 1:M) {
    
    #Generate stratified sample from quantile function Q (using distribution parameters in ... input)
    #Method uses inverse-transformation sampling with PB as a quantile-stratified sample from uniform distribution
    #Note: To avoid calculation error, we do not allow specification of lower.tail or log.p in ... inputs for Q
    nnn <- LAYERS[m]
    PP <- runif(nnn)
    PB <- (1:nnn-PP)/nnn
    ARGS <- list(...)
    ARGS[[prob.arg]] <- PB
    ARGS$lower.tail  <- NULL
    ARGS$log.p       <- NULL
    QQ <- do.call(Q, args = ARGS)
    
    #Check validity of blocked sample (checking quantile function Q)
    if(!is.vector(QQ))                        stop('Error: Input Q should give vector output')
    if(!is.numeric(QQ))                       stop('Error: Input Q should give numeric output')
    if (nnn > 1) {
      for(r in 2:nnn) {
        if (QQ[r] < QQ[r-1])                    stop('Error: Input Q should be non-decreasing') } }
    
    #Add sample to sample list
    SAMPLES[[m]] <- QQ }
  
  #Combine and randomly permute sample layers
  PERM <- sample(1:nn, replace = FALSE)
  OUT <- unlist(SAMPLES)[PERM]
  
  #Give output
  OUT }
