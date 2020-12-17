#' Sample Kurtosis
#'
#' \code{kurtosis} returns the sample kurtosis of a data vector/matrix
#'
#' This function computes the sample kurtosis for a data vector or matrix.  For a vector input the function returns a single value for the
#' sample kurtosis of the data.  For a matrix input the function treats each column as a data vector and returns a vector of values for the
#' sample kurtosis of each of these datasets.  The function can compute different types of kurtosis statistics using the \code{kurt.type} input.
#'
#' @usage \code{kurtosis}
#' @param x A data vector/matrix
#' @param kurt.type The type of kurtosis statistic used ('Moment', 'Fisher Pearson' or 'Adjusted Fisher Pearson')
#' @param kurt.excess Logical value; if \code{TRUE} the function gives the excess kurtosis (instead of raw kurtosis)
#' @param na.rm Logical value; if \code{TRUE} the function removes \code{NA} values
#' @return The sample kurtosis of the data vector/matrix

kurtosis <- function(x, kurt.type = NULL, kurt.excess = FALSE, na.rm = FALSE) {

  #Check inputs
  if ((!is.vector(x))&(!is.matrix(x)))              { stop('Error: Input x should be a data vector or matrix') }
  if (!is.numeric(x))                               { stop('Error: Input x should be numeric') }
  if (missing(kurt.type)) { kurt.type <- 'Moment' }
  TYPES <- c('Moment', 'Fisher Pearson', 'Adjusted Fisher Pearson', 'Minitab', 'Excel', 'SPSS', 'SAS', 'Stata')
  if (!(kurt.type %in% TYPES))                      { stop('Error: Input kurt.type not recognised') }
  if (!is.vector(kurt.excess))                      { stop('Error: Input kurt.excess should be a single logical value') }
  if (!is.logical(kurt.excess))                     { stop('Error: Input kurt.excess should be a single logical value') }
  if (length(kurt.excess) != 1)                     { stop('Error: Input kurt.excess should be a single logical value') }
  if (!is.vector(na.rm))                            { stop('Error: Input na.rm should be a single logical value') }
  if (!is.logical(na.rm))                           { stop('Error: Input na.rm should be a single logical value') }
  if (length(na.rm) != 1)                           { stop('Error: Input na.rm should be a single logical value') }

  #Set kurt adjustment function
  kurt.adj <- function(n) {
    B <- 1
    if (kurt.type %in% c('Adjusted Fisher Pearson', 'Minitab', 'Excel', 'SPSS', 'SAS')) {
      B <- (n+1)*n^2/((n-1)*(n-2)*(n-3)) }
    if (kurt.type %in% c('Fisher Pearson', 'Stata')) {
      B <- (n/(n-1))^2 }
    B }

  #Set excess adjustment function
  excess.adj <- function(n) {
    C <- -3*kurt.excess
    if (kurt.type %in% c('Adjusted Fisher Pearson', 'Minitab', 'Excel', 'SPSS', 'SAS')) {
      C <- -3*kurt.excess*(n-1)^2/((n-2)*(n-3)) }
    C }

  #Compute the sample kurtosis
  if (is.vector(x)) { x <- as.matrix(x, ncols = 1) }
  m   <- ncol(x)
  OUT <- numeric(m)
  for (i in 1:m) {
    if (na.rm) { xx  <- x[!is.na(x[,i]),i] } else { xx <- x[,i] }
    n   <- length(xx)
    SS  <- sum((xx - mean(xx))^2)
    SQ  <- sum((xx - mean(xx))^4)
    OUT[i] <- ifelse(n >= 2, kurt.adj(n)*n*SQ/SS^2 + excess.adj(n), NA) }

  #Give output
  OUT }





