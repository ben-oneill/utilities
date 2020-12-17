#' Sample Skewness
#'
#' \code{skewness} returns the sample skewness of a data vector/matrix
#'
#' This function computes the sample skewness for a data vector or matrix.  For a vector input the function returns a single value for the
#' sample skewness of the data.  For a matrix input the function treats each column as a data vector and returns a vector of values for the
#' sample skewness of each of these datasets.  The function can compute different types of skewness statistics using the \code{skew.type} input.
#'
#' @usage \code{skewness}
#' @param x A data vector/matrix
#' @param skew.type The type of skewness statistic used ('Moment', 'Fisher Pearson' or 'Adjusted Fisher Pearson')
#' @param na.rm Logical value; if \code{TRUE} the function removes \code{NA} values
#' @return The sample skewness of the data vector/matrix

skewness <- function(x, skew.type = NULL, na.rm = FALSE) {

  #Check inputs
  if ((!is.vector(x))&(!is.matrix(x)))              { stop('Error: Input x should be a data vector or matrix') }
  if (!is.numeric(x))                               { stop('Error: Input x should be numeric') }
  if (missing(skew.type)) { skew.type <- 'Moment' }
  TYPES <- c('Moment', 'Fisher Pearson', 'Adjusted Fisher Pearson', 'Minitab', 'Excel', 'SPSS', 'SAS', 'Stata')
  if (!(skew.type %in% TYPES))                      { stop('Error: Input skew.type not recognised') }
  if (!is.vector(na.rm))                            { stop('Error: Input na.rm should be a single logical value') }
  if (!is.logical(na.rm))                           { stop('Error: Input na.rm should be a single logical value') }
  if (length(na.rm) != 1)                           { stop('Error: Input na.rm should be a single logical value') }

  #Set skew adjustment
  skew.adj <- function(n) {
    A <- 1
    if (skew.type %in% c('Adjusted Fisher Pearson', 'Minitab', 'Excel', 'SPSS', 'SAS')) { A <- (n^2)/((n-1)*(n-2)) }
    if (skew.type %in% c('Fisher Pearson', 'Stata')) { A <- (n/(n-1))^(3/2) }
    A }

  #Compute the sample skewness
  if (is.vector(x)) { x <- as.matrix(x, ncols = 1) }
  m   <- ncol(x)
  OUT <- numeric(m)
  for (i in 1:m) {
    if (na.rm) { xx  <- x[!is.na(x[,i]),i] } else { xx <- x[,i] }
    n   <- length(xx)
    SS  <- sum((xx - mean(xx))^2)
    SC  <- sum((xx - mean(xx))^3)
    OUT[i] <- ifelse(n >= 2, skew.adj(n)*sqrt(n)*SC/SS^(3/2), NA) }

  #Give output
  OUT }
