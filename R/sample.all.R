#' All Sampling Variations/Permutations
#'
#' \code{sample.all} returns a matrix of all sampling variations/permutations from a set of integers
#'
#' This function computes all sample vectors of size \code{size} composed of the elements \code{1, ..., n}, either with or without replacement
#' of elements.  If \code{size = n} and \code{replace = TRUE} then the list of all sample vectors corresponds to a list of all permutations of
#' the integers \code{1, ..., n}.
#'
#' @param n Number of integers to sample from
#' @param size Length of the sample vectors
#' @param replace Logical value; if \code{FALSE} the sampling is without replacement; if \code{TRUE} the sampling is with replacement
#' @param prob Probability vector giving the sampling probability for each element (must be a probability vector with length \code{n})
#' @return A matrix of all permutations of the elements \code{1, ..., n} (rows of the matrix give the permutations)
#'
#' @examples
#' sample.all(n = 4, replace = FALSE)
sample.all <- function(n, size = n, replace = FALSE, prob = NULL) {

  #Check input n
  if (!is.vector(n))                                   { stop('Error: Input n should be a positive integer') }
  if (!is.numeric(n))                                  { stop('Error: Input n should be a positive integer') }
  if (length(n) != 1)                                  { stop('Error: Input n should be a single positive integer') }
  if (as.integer(n) != n)                              { stop('Error: Input n should be a positive integer') }
  if (n < 1)                                           { stop('Error: Input n should be a positive integer') }

  #Check input size
  if (!is.vector(size))                                { stop('Error: Input size should be a positive integer') }
  if (!is.numeric(size))                               { stop('Error: Input size should be a positive integer') }
  if (length(size) != 1)                               { stop('Error: Input size should be a single positive integer') }
  if (as.integer(size) != size)                        { stop('Error: Input size should be a positive integer') }
  if (size < 1)                                        { stop('Error: Input size should be a positive integer') }

  #Check input replace
  if (!is.vector(replace))                             { stop('Error: Input replace should be a logical value') }
  if (!is.logical(replace))                            { stop('Error: Input replace should be a logical value') }
  if (length(replace) != 1)                            { stop('Error: Input replace should be a logical value') }

  #Check input prob
  PROB <- FALSE
  if (!is.null(prob)) {
    PROB <- TRUE
    if (!is.vector(prob))                              { stop('Error: Input prob should be a probability vector with length n') }
    if (length(prob) != n)                             { stop('Error: Input prob should be a probability vector with length n') }
    if (min(prob) < 0)                                 { stop('Error: Input prob should be a probability vector with length n') }
    SUM  <- sum(prob)
    if (SUM != 1) { prob <- prob/SUM } }

  #Check matrix size and give warnings
  if (replace)  { ELEMENTS <- (n^size)*(size + PROB)  }
  if (!replace) { ELEMENTS <- (factorial(n)/factorial(n-min(n, size)))*(min(n, size) + PROB) }
  if (ELEMENTS >= 2^30)                                { stop('Error: There are too many elements for the output matrix') }
  if (ELEMENTS >= 2^25) {
    PROCEED <- readline(paste0('This function will create a large matrix containing ', ELEMENTS,
                               ' elements --- Are you sure you want to proceed?  Y/N \n'))
    while (!(PROCEED %in% c('Y', 'y', 'N', 'n'))) {
      PROCEED <- readline('    Invalid response --- Are you sure you want to proceed?  Y/N \n') }
      if (PROCEED %in% c('N', 'n'))                    { stop('Function terminated at user request') } }

  #Show all sample vectors with replacement
  if (replace) {

    #Create output matrix
    N   <- n^size
    OUT <- matrix(0, nrow = N, ncol = size)
    rownames(OUT) <- sprintf('Sample[%s]', 1:N)
    colnames(OUT) <- sprintf('X[%s]', 1:size)
    if (PROB) { OUT <- cbind(OUT, Probability = 1) }

    #Fill in the output matrix
    for (i in 1:size) { OUT[, i] <- rep(1:n, each = n^(size-i)) }

    #Fill in the probabilities
    if (PROB) {
    for (i in 1:N) {
      OUT[i, size+1] <- prod(prob[OUT[i, 1:size]]) } } }

  #Show all sample vectors without replacement (permutations)
  if (!replace) {

    #Adjust size if it is too large
    size <- min(n, size)

    #Create output and index matrices
    N   <- factorial(n)/factorial(n-size)
    OUT <- matrix(0, nrow = N, ncol = size)
    rownames(OUT) <- sprintf('Sample[%s]', 1:N)
    colnames(OUT) <- sprintf('X[%s]', 1:size)
    if (PROB) { OUT <- cbind(OUT, Probability = 1) }
    for (i in 1:size)   { OUT[, i] <- rep(1:(n-i+1), each = factorial(n-i)/factorial(n-size)) }

    #Fill in the output matrix
    for (i in 1:N) {
      VALS <- 1:n
      for (j in 1:size) {
        index    <- OUT[i,j]
        OUT[i,j] <- VALS[index]
        if (PROB) { OUT[i, size+1] <- OUT[i, size+1]*prob[OUT[i, j]]/sum(prob[VALS]) }
        VALS     <- VALS[-index] } } }

  #Return output
  OUT }
