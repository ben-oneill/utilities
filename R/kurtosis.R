#' Sample Kurtosis
#'
#' \code{kurtosis} returns the sample kurtosis of a data vector/matrix
#'
#' This function computes the sample kurtosis for a data vector or matrix.  For a vector input the function returns a single value for the
#' sample kurtosis of the data.  For a matrix input the function treats each column as a data vector and returns a vector of values for the
#' sample kurtosis of each of these datasets.  The function can compute different types of kurtosis statistics using the \code{kurt.type} input.
#'
#' @param x A data vector/matrix
#' @param kurt.type The type of kurtosis statistic used ('Moment', 'Fisher Pearson' or 'Adjusted Fisher Pearson')
#' @param kurt.excess Logical value; if \code{TRUE} the function gives the excess kurtosis (instead of raw kurtosis)
#' @param na.rm Logical value; if \code{TRUE} the function removes \code{NA} values
#' @return The sample kurtosis of the data vector/matrix
#'
#' @examples
#' kurtosis(rnorm(1000))
#' kurtosis(rexp(1000))
kurtosis <- function(x, kurt.type = NULL, kurt.excess = FALSE, na.rm = FALSE) {

  #Check inputs
  if ((!is.vector(x))&(!is.matrix(x)))              { stop('Error: Input x should be a data vector or matrix') }
  if (!is.numeric(x))                               { stop('Error: Input x should be numeric') }

  #Check input kurt.type
  #Types are as follows:
  #  b = Moment (Minitab)
  #  g = Fisher Pearson (moments package in R, STATA)
  #  G = Adjusted Fisher Pearson (Excel, SPSS, SAS)
  if (missing(kurt.type)) { kurt.type <- 'Fisher Pearson' }
  TYPES <- c('Moment', 'Fisher Pearson', 'Adjusted Fisher Pearson', 'b', 'g', 'G', 'Minitab', 'Excel', 'SPSS', 'SAS', 'Stata')
  if (!(kurt.type %in% TYPES))                      { stop('Error: Input kurt.type not recognised') }

  #Check input kurt.excess
  if (!is.vector(kurt.excess))                      { stop('Error: Input kurt.excess should be a single logical value') }
  if (!is.logical(kurt.excess))                     { stop('Error: Input kurt.excess should be a single logical value') }
  if (length(kurt.excess) != 1)                     { stop('Error: Input kurt.excess should be a single logical value') }

  #Check input na.rm
  if (!is.vector(na.rm))                            { stop('Error: Input na.rm should be a single logical value') }
  if (!is.logical(na.rm))                           { stop('Error: Input na.rm should be a single logical value') }
  if (length(na.rm) != 1)                           { stop('Error: Input na.rm should be a single logical value') }

  #Set kurt adjustment function
  #Default type with no adjustment is 'Fisher Pearson'
  kurt.adj <- function(n) {
    B <- 1
    if (kurt.type %in% c('Moment', 'b', 'Minitab')) {
      B <- (n/(n-1))^2 }
    if (kurt.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
      B <- (n+1)*n^2/((n-1)*(n-2)*(n-3)) }
    B }

  #Set excess adjustment function
  #Default type with no adjustment is 'Fisher Pearson'
  excess.adj <- function(n) {
    C <- -3*kurt.excess
    if (kurt.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
      C <- -3*kurt.excess*(n-1)^2/((n-2)*(n-3)) }
    C }

  #Compute the sample kurtosis
  if (is.vector(x)) { x <- as.matrix(x, ncols = 1) }
  m   <- ncol(x)
  OUT <- numeric(m)
  for (i in 1:m) {
    if (na.rm) { xx  <- x[!is.na(x[,i]),i] } else { xx <- x[,i] }
    n  <- length(xx)
    MM <- rep(0, n)
    SS <- rep(0, n)
    SC <- rep(0, n)
    SQ <- rep(0, n)
    MM[1] <- xx[1]
    if (n > 1) {
    for (k in 2:n) {
      MM[k] <- ((k-1)*MM[k-1] + xx[k])/k
      DD    <- MM[k-1] - xx[k]
      SS[k] <- SS[k-1] + ((k-1)/k)*DD^2
      SC[k] <- SC[k-1] + ((3*SS[k-1])/k)*DD - ((k-1)*(k-2)/k^2)*DD^3
      SQ[k] <- SQ[k-1] + ((4*SC[k-1])/k)*DD + ((6*SS[k-1])/k^2)*DD^2 + ((k-1)*(1+(k-1)^3)/k^4)*DD^4 } }
    SS <- SS[n]
    SC <- SC[n]
    SQ <- SQ[n]
    OUT[i] <- ifelse(n >= 2, kurt.adj(n)*n*SQ/SS^2 + excess.adj(n), NA) }

  #Give output
  OUT }
