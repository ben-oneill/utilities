#' All Sampling Variations/Permutations
#'
#' \code{sample.all} returns a matrix of all sampling variations/permutations from a set of integers
#'
#' This function computes all sample vectors of size \code{size} composed of the elements \code{1, ..., n}, either with or without replacement
#' of elements.  If \code{size = n} and \code{replace = TRUE} then the list of all sample vectors corresponds to a list of all permutations of
#' the integers \code{1, ..., n}.
#'
#' @usage \code{sample.all}
#' @param n Number of integers to sample from
#' @param size Length of the sample vectors
#' @param replace Logical value; if \code{FALSE} the sampling is without replacement; if \code{TRUE} the sampling is with replacement
#' @return A matrix of all permutations of the elements \code{1, ..., n} (rows of the matrix give the permutations)

sample.all <- function(n, size = n, replace = FALSE) {

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

  #Show all sample vectors with replacement
  if (replace) {

    #Create output matrix
    N   <- n^size
    OUT <- matrix(0, nrow = N, ncol = size)

    #Fill in the output matrix
    for (i in 1:size) { OUT[, i] <- rep(1:n, each = n^(size-i)) } }

  #Show all sample vectors without replacement (permutations)
  if (!replace) {

    #Adjust size if it is too large
    size <- min(n, size)

    #Create output and index matrices
    N   <- factorial(n)/factorial(n-size)
    OUT <- matrix(0, nrow = N, ncol = size)
    IND <- matrix(0, nrow = N, ncol = size)
    for (i in 1:size) { IND[, i] <- rep(1:(n-i+1), each = factorial(n-i)/factorial(n-size)) }

    #Fill in the output matrix
    for (i in 1:N) {
      VALS <- 1:n
      for (j in 1:size) {
        index    <- IND[i,j]
        OUT[i,j] <- VALS[index]
        VALS     <- VALS[-index] } } }

  #Return output
  OUT }
