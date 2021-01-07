#' Sample Moments
#'
#' \code{moments} returns the sample moments of a data vector/matrix
#'
#' This function computes the sample moments for a data vector, matrix or list (sample mean, sample variance, sample skewness and sample kurtosis).
#' For a vector input the function returns a single value for each sample moment of the data.  For a matrix or list input the function treats each
#' column/element as a data vector and returns a matrix of values for the sample moments of each of these datasets.  The function can compute
#' different types of skewness and kurtosis statistics using the \code{skew.type}, \code{kurt.type} and \code{kurt.excess} inputs.  (For details
#' on the different types of skewness and kurtosis statistics, see Joanes and Gill 1998.)
#'
#' @param x A data vector/matrix/list
#' @param skew.type The type of kurtosis statistic used ('Moment', 'Fisher Pearson' or 'Adjusted Fisher Pearson')
#' @param kurt.type The type of kurtosis statistic used ('Moment', 'Fisher Pearson' or 'Adjusted Fisher Pearson')
#' @param kurt.excess Logical value; if \code{TRUE} the function gives the excess kurtosis (instead of raw kurtosis)
#' @param na.rm Logical value; if \code{TRUE} the function removes \code{NA} values
#' @param include.sd Logical value; if \code{TRUE} the output includes a column for the sample standard deviation (if needed)
#' @return A data frame containing the sample moments of the data vector/matrix
#'
#' @examples
#' #Create some subgroups of mock data and a pooled dataset
#' set.seed(1)
#' N    <- c(28, 44, 51)
#' SUB1 <- rnorm(N[1])
#' SUB2 <- rnorm(N[2])
#' SUB3 <- rnorm(N[3])
#' DATA <- list(Subgroup1 = SUB1, Subgroup2 = SUB2, Subgroup3 = SUB3)
#' POOL <- c(SUB1, SUB2, SUB3)
#'
#' #Compute sample moments for subgroups and pooled data
#' MOMENTS <- moments(DATA)
#' POOLMOM <- moments(POOL)
#'
#' #Compute pooled moments via sample decomposition
#' sample.decomp(moments = MOMENTS)

moments <- function(x, skew.type = NULL, kurt.type = NULL, kurt.excess = FALSE, na.rm = TRUE, include.sd = FALSE) {

  #Extract the data name
  DATANAME <- deparse(substitute(x))

  #Check input x
  if (is.list(x))  {
    m <- length(x)
    for (i in 1:m) {
      if (!is.numeric(x[[i]]))                      { stop(paste0('Error: Element ', i, ' of input list x should be numeric')) } } }
  if (!is.list(x)) {
    if (!is.numeric(x))                             { stop('Error: Input x should be numeric') } }

  #Check inputs skew.type and kurt.type
  #Types are as follows:
  #  b = Moment (Minitab)
  #  g = Fisher Pearson (moments package in R, STATA)
  #  G = Adjusted Fisher Pearson (Excel, SPSS, SAS)
  TYPES <- c('Moment', 'Fisher Pearson', 'Adjusted Fisher Pearson', 'b', 'g', 'G', 'Minitab', 'Excel', 'SPSS', 'SAS', 'Stata')
  if (missing(skew.type)) { skew.type <- 'Fisher Pearson' }
  if (missing(kurt.type)) { kurt.type <- 'Fisher Pearson' }
  if (!(skew.type %in% TYPES))                      { stop('Error: Input skew.type not recognised') }
  if (!(kurt.type %in% TYPES))                      { stop('Error: Input kurt.type not recognised') }

  #Check input kurt.excess
  if (!is.vector(kurt.excess))                      { stop('Error: Input kurt.excess should be a single logical value') }
  if (!is.logical(kurt.excess))                     { stop('Error: Input kurt.excess should be a single logical value') }
  if (length(kurt.excess) != 1)                     { stop('Error: Input kurt.excess should be a single logical value') }

  #Check input na.rm
  if (!is.vector(na.rm))                            { stop('Error: Input na.rm should be a single logical value') }
  if (!is.logical(na.rm))                           { stop('Error: Input na.rm should be a single logical value') }
  if (length(na.rm) != 1)                           { stop('Error: Input na.rm should be a single logical value') }

  #Check input include.sd
  if (!is.vector(include.sd))                       { stop('Error: Input include.sd should be a single logical value') }
  if (!is.logical(include.sd))                      { stop('Error: Input include.sd should be a single logical value') }
  if (length(include.sd) != 1)                      { stop('Error: Input include.sd should be a single logical value') }

  #Set skew and kurt adjustments
  #Default type with no adjustment is 'Fisher Pearson'
  skew.adj <- function(n) {
    A <- 1
    if (skew.type %in% c('Moment', 'b', 'Minitab')) {
      A <- (n/(n-1))^(3/2) }
    if (skew.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
      A <- (n^2)/((n-1)*(n-2)) }
    A }
  kurt.adj <- function(n) {
    B <- 1
    if (kurt.type %in% c('Moment', 'b', 'Minitab')) {
      B <- (n/(n-1))^2 }
    if (kurt.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
      B <- (n+1)*n^2/((n-1)*(n-2)*(n-3)) }
    B }
  excess.adj <- function(n) {
    C <- -3*kurt.excess
    if (kurt.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
      C <- -3*kurt.excess*(n-1)^2/((n-2)*(n-3)) }
    C }

  #Convert data to list
  if (is.list(x)) {
    DATA  <- x
    NAMES <- names(x) }
  if ((!is.list(x))&(is.vector(x))) {
    DATA  <- lapply(1, function(i) x)
    NAMES <- DATANAME
    names(DATA) <- NAMES }
  if ((!is.list(x))&(is.matrix(x))) {
    DATA  <- lapply(seq_len(ncol(x)), function(i) x[,i])
    NAMES <- colnames(x)
    names(DATA) <- NAMES }

  #Create output data frame
  m <- length(DATA)
  if (include.sd) {
    OUT <- data.frame(n = rep(0, m), sample.mean = rep(0, m), sample.sd = rep(0, m), sample.var = rep(0, m),
                      sample.skew = rep(0, m), sample.kurt = rep(0, m), NAs = rep(0, m)) } else {
    OUT <- data.frame(n = rep(0, m), sample.mean = rep(0, m), sample.var = rep(0, m),
                      sample.skew = rep(0, m), sample.kurt = rep(0, m), NAs = rep(0, m)) }
  if (m == 1) { rownames(OUT) <- DATANAME } else {
    if (is.null(NAMES)) { rownames(OUT) <- paste0(DATANAME, sprintf('[%s]', 1:m)) } else { rownames(OUT) <- NAMES } }
  class(OUT) <- c('moments', 'data.frame')
  attr(OUT, 'skew.type')   <- skew.type
  attr(OUT, 'kurt.type')   <- kurt.type
  attr(OUT, 'kurt.excess') <- kurt.excess

  #Compute the sample moments
  for (i in 1:m) {
    NAS <- sum(is.na(DATA[[i]]))
    if (na.rm) { xx <- DATA[[i]][!is.na(DATA[[i]])] } else { xx <- DATA[[i]] }
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
    MM <- MM[n]
    SS <- SS[n]
    SC <- SC[n]
    SQ <- SQ[n]
    sample.mean <- MM
    sample.var  <- ifelse(n >= 2, SS/(n-1), NA)
    sample.skew <- ifelse(n >= 2, skew.adj(n)*sqrt(n)*SC/SS^(3/2), NA)
    sample.kurt <- ifelse(n >= 2, kurt.adj(n)*n*SQ/SS^2 + kurt.excess, NA)

    #Add to output
    OUT$n[i] <- n
    OUT$sample.mean[i] <- sample.mean
    if (include.sd) { OUT$sample.sd[i] <- sqrt(sample.var) }
    OUT$sample.var[i]  <- sample.var
    OUT$sample.skew[i] <- sample.skew
    OUT$sample.kurt[i] <- sample.kurt
    OUT$NAs[i]         <- NAS }

  #Give output
  OUT }
