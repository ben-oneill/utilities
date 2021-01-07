#' Sample Skewness
#'
#' \code{skewness} returns the sample skewness of a data vector/matrix
#'
#' This function computes the sample skewness for a data vector or matrix.  For a vector input the function returns a single value for the
#' sample skewness of the data.  For a matrix input the function treats each column as a data vector and returns a vector of values for the
#' sample skewness of each of these datasets.  The function can compute different types of skewness statistics using the \code{skew.type} input.
#'
#' @param x A data vector/matrix
#' @param skew.type The type of skewness statistic used ('Moment', 'Fisher Pearson' or 'Adjusted Fisher Pearson')
#' @param na.rm Logical value; if \code{TRUE} the function removes \code{NA} values
#' @return The sample skewness of the data vector/matrix
#'
#' @examples
#'
#' skewness(rnorm(1000))
#' skewness(rexp(1000))
#'
skewness <- function(x, skew.type = NULL, na.rm = FALSE) {

  #Check input x
  if ((!is.vector(x))&(!is.matrix(x)))              { stop('Error: Input x should be a data vector or matrix') }
  if (!is.numeric(x))                               { stop('Error: Input x should be numeric') }

  #Check input skew.type
  #Types are as follows:
  #  b = Moment (Minitab)
  #  g = Fisher Pearson (moments package in R, STATA)
  #  G = Adjusted Fisher Pearson (Excel, SPSS, SAS)
  if (missing(skew.type)) { skew.type <- 'Fisher Pearson' }
  TYPES <- c('Moment', 'Fisher Pearson', 'Adjusted Fisher Pearson', 'b', 'g', 'G', 'Minitab', 'Excel', 'SPSS', 'SAS', 'Stata')
  if (!(skew.type %in% TYPES))                      { stop('Error: Input skew.type not recognised') }

  #Check input na.rm
  if (!is.vector(na.rm))                            { stop('Error: Input na.rm should be a single logical value') }
  if (!is.logical(na.rm))                           { stop('Error: Input na.rm should be a single logical value') }
  if (length(na.rm) != 1)                           { stop('Error: Input na.rm should be a single logical value') }

  #Set skew adjustment
  #Default type with no adjustment is 'Fisher Pearson'
  skew.adj <- function(n) {
    A <- 1
    if (skew.type %in% c('Moment', 'b', 'Minitab')) {
      A <- (n/(n-1))^(3/2) }
    if (skew.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
      A <- (n^2)/((n-1)*(n-2)) }
    A }

  #Compute the sample skewness
  if (is.vector(x)) { x <- as.matrix(x, ncols = 1) }
  m   <- ncol(x)
  OUT <- numeric(m)
  for (i in 1:m) {
    if (na.rm) { xx  <- x[!is.na(x[,i]),i] } else { xx <- x[,i] }
    n  <- length(xx)
    MM <- rep(0, n)
    SS <- rep(0, n)
    SC <- rep(0, n)
    MM[1] <- xx[1]
    if (n > 1) {
    for (k in 2:n) {
      MM[k] <- ((k-1)*MM[k-1] + xx[k])/k
      DD    <- MM[k-1] - xx[k]
      SS[k] <- SS[k-1] + ((k-1)/k)*DD^2
      SC[k] <- SC[k-1] + ((3*SS[k-1])/k)*DD - ((k-1)*(k-2)/k^2)*DD^3 } }
    SS <- SS[n]
    SC <- SC[n]
    OUT[i] <- ifelse(n >= 2, skew.adj(n)*sqrt(n)*SC/SS^(3/2), NA) }

  #Give output
  OUT }
