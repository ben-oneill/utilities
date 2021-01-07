#' Sample decomposition
#'
#' \code{sample.decomp} returns the data-frame of sample statistics for sample groups and their pooled sample
#'
#' It is often useful to take a set of sample groups with known sample statistics and aggregate these into a single pooled sample and find the
#' sample statistics of the pooled sample.  Likewise, it is sometimes useful to take a set of sample groups and a pooled group with known sample
#' statistics and determine the statistics of the other group required to complete the pooled sample.  Both of these tasks can be accomplished
#' using decomposition formulae for the sample size, sample mean and sample variance (or sample standard deviation).  This function implements
#' either of these two decomposition methods to find the sample statistics of the pooled sample or the other group remaining to obtain the pooled
#' sample.  The user inputs vectors for the sample size, sample mean and sample variance (or sample standard deviation).  By default the groups
#' are taken to be separate groups and the function computes the sample statistics for the pooled sample  However, the user can input the number
#' pooled sample as the input \code{pooled}; in this case that group is treated as the pooled sample and the function computes the other sample
#' group required to obtain this pooled sample.  The function returns a data-frame showing the sample statistics for all the groups including
#' the pooled sample.
#'
#' @param moments A data-frame of moments (an object of class 'moments')
#' @param n A vector of sample sizes
#' @param sample.mean A vector of sample means
#' @param sample.sd A vector of sample standard deviations
#' @param sample.var A vector of sample variances
#' @param sample.skew A vector of sample skewness
#' @param sample.kurt A vector of sample kurotsis
#' @param names A vector of names for the sample groups
#' @param pooled The number of the pooled group (if the pooled group is already present)
#' @param skew.type The type of skewness statistic used ('Moment', 'Fisher Pearson' or 'Adjusted Fisher Pearson')
#' @param kurt.type The type of kurtosis statistic used ('Moment', 'Fisher Pearson' or 'Adjusted Fisher Pearson')
#' @param kurt.excess Logical value; if \code{TRUE} the sample kurtosis is the excess kurtosis (instead of the raw kurtosis)
#' @param include.sd Logical value; if \code{TRUE} the output includes a column for the sample standard deviation (if needed)
#' @return A data-frame of all groups showing their sample sizes and sample moments
#'
#' @seealso \code{\link{moments}}
sample.decomp <- function(moments = NULL, n = NULL,
                          sample.mean = NULL, sample.sd = NULL, sample.var = NULL, sample.skew = NULL, sample.kurt = NULL,
                          names = NULL, pooled = NULL, skew.type = NULL, kurt.type = NULL, kurt.excess = NULL, include.sd = FALSE) {

  #Check input moments
  if (!is.null(moments)) {
    if (!('moments' %in% class(moments)))           { stop('Error: Input moments must be a moments object') }
    if (is.null(n))           { n           <- moments$n           } else { stop('Error: Input moments or descriptive statistics but not both') }
    if (is.null(sample.mean)) { sample.mean <- moments$sample.mean } else { stop('Error: Input moments or descriptive statistics but not both') }
    if ('sample.sd' %in% colnames(moments)) {
    if (is.null(sample.sd))   { sample.sd   <- moments$sample.sd   } else { stop('Error: Input moments or descriptive statistics but not both') } }
    if (is.null(sample.var))  { sample.var  <- moments$sample.var  } else { stop('Error: Input moments or descriptive statistics but not both') }
    if (is.null(sample.skew)) { sample.skew <- moments$sample.skew } else { stop('Error: Input moments or descriptive statistics but not both') }
    if (is.null(sample.kurt)) { sample.kurt <- moments$sample.kurt } else { stop('Error: Input moments or descriptive statistics but not both') }
    if (is.null(skew.type))   { skew.type   <- attributes(moments)$skew.type   } else { stop('Error: Input moments or skew.type but not both') }
    if (is.null(kurt.type))   { kurt.type   <- attributes(moments)$kurt.type   } else { stop('Error: Input moments or kurt.type but not both') }
    if (is.null(kurt.excess)) { kurt.excess <- attributes(moments)$kurt.excess } else { stop('Error: Input moments or kurt.excess but not both') }
    if (is.null(names))       { names       <- rownames(moments) } }

  #Check input n
  if (is.null(n))                                   { stop('Error: You must input n') }
  if (!is.vector(n))                                { stop('Error: Input n should be a vector of positive integers') }
  if (!is.numeric(n))                               { stop('Error: Input n should be a vector of positive integers') }
  if (any(as.integer(n) != n))                      { stop('Error: Input n should be a vector of positive integers') }
  if (min(n) < 1)                                   { stop('Error: Input n should be a vector of positive integers') }
  N <- length(n)

  #Check which moments are specified
  MOMENTS <- rep(FALSE, 4)
  if (!is.null(sample.mean)) { MOMENTS[1] <- TRUE }
  if (!is.null(sample.var))  { MOMENTS[2] <- TRUE }
  if (!is.null(sample.sd))   { MOMENTS[2] <- TRUE }
  if (!is.null(sample.skew)) { MOMENTS[3] <- TRUE }
  if (!is.null(sample.kurt)) { MOMENTS[4] <- TRUE }

  #Set maximum specified moment
  MAXMOM <- 0
  if (MOMENTS[1])              { MAXMOM <- 1 }
  if (prod(MOMENTS[1:2] == 1)) { MAXMOM <- 2 }
  if (prod(MOMENTS[1:3] == 1)) { MAXMOM <- 3 }
  if (prod(MOMENTS[1:4] == 1)) { MAXMOM <- 4 }

  #Check input sample.mean
  if (MOMENTS[1]) {
    if (!is.vector(sample.mean))                    { stop('Error: Input sample.mean should be a numeric vector (if specified)') }
    if (!is.numeric(sample.mean))                   { stop('Error: Input sample.mean should be a numeric vector (if specified)') }
    if (length(sample.mean) != N)                   { stop('Error: Input sample.mean must have the same length as n (if specified)') } }

  #Check inputs sample.sd and sample.var
  if (!is.null(sample.var)) {
    if (!is.vector(sample.var))                     { stop('Error: Input sample.var should be a numeric vector (if specified)') }
    if (!is.numeric(sample.var))                    { stop('Error: Input sample.var should be a numeric vector (if specified)') }
    if (length(sample.var) != N)                    { stop('Error: Input sample.var must have the same length as n (if specified)') }
    if (min(sample.var) < 0)                        { stop('Error: Values in sample.var cannot be negative') } }
  if (!is.null(sample.sd)) {
    if (!is.vector(sample.sd))                      { stop('Error: Input sample.sd should be a numeric vector (if specified)') }
    if (!is.numeric(sample.sd))                     { stop('Error: Input sample.sd should be a numeric vector (if specified)') }
    if (length(sample.sd) != N)                     { stop('Error: Input sample.sd must have the same length as n (if specified)') }
    if (min(sample.sd) < 0)                         { stop('Error: Values in sample.sd cannot be negative') } }
  if ((!is.null(sample.sd))&(!is.null(sample.var))) {
    ERROR     <- sum((sample.sd^2 - sample.var)^2)
    THRESHOLD <- (1e-10)*min(sample.sd)
    if (ERROR > THRESHOLD) {
      stop('Error: You may specify sample.sd or sample.var but not both') } else {
      warning('You have specified sample.sd and sample.var --- only sample.var is used for calculations') } }
  if ((!is.null(sample.sd))&(is.null(sample.var))) { sample.var <- sample.sd^2 }

  #Check input sample.skew
  if (MOMENTS[3]) {
    if (!is.vector(sample.skew))                    { stop('Error: Input sample.skew should be a numeric vector (if specified)') }
    if (!is.numeric(sample.skew))                   { stop('Error: Input sample.skew should be a numeric vector (if specified)') }
    if (length(sample.skew) != N)                   { stop('Error: Input sample.skew must have the same length as n (if specified)') } }

  #Check input sample.kurt
  if (MOMENTS[4]) {
    if (!is.vector(sample.kurt))                    { stop('Error: Input sample.kurt should be a numeric vector (if specified)') }
    if (!is.numeric(sample.kurt))                   { stop('Error: Input sample.kurt should be a numeric vector (if specified)') }
    if (length(sample.kurt) != N)                   { stop('Error: Input sample.kurt must have the same length as n (if specified)') } }

  #Check input names
  if (!is.null(names)) {
    if (!is.vector(names))                          { stop('Error: Input names should be a character vector (if specified)') }
    if (!is.character(names))                       { stop('Error: Input names should be a character vector (if specified)') }
    if (length(names) != N)                         { stop('Error: Input names must have the same length as n (if specified)') }
    if (length(unique(names)) != N)                 { stop('Error: Input names should have unique names for the sample groups') } }

  #Check input pooled
  if (!is.null(pooled)) {
    if (!is.vector(pooled))                         { stop('Error: Input pooled should be a single integer (if specified)') }
    if (!is.numeric(pooled))                        { stop('Error: Input pooled should be a single integer (if specified)') }
    if (length(pooled) != 1)                        { stop('Error: Input pooled should be a single integer (if specified)') }
    if (!(pooled %in% 1:N))                         { stop('Error: Input pooled must specify a sample group by number (if specified)') }
    if ('--other--' %in% names)                     { warning('The name --other-- is included in your names vector --- this is confusing') }
    if (n[pooled] <= sum(n[-pooled]))               { stop('Error: The size of the pooled group should be larger than the total size of the other groups') } }

  #Check inputs skew.type, kurt.type and kurt.excess
  #Types are as follows:
  #  b = Moment (Minitab)
  #  g = Fisher Pearson (moments package in R, STATA)
  #  G = Adjusted Fisher Pearson (Excel, SPSS, SAS)
  if (MAXMOM >= 3) {
    TYPES <- c('Moment', 'Fisher Pearson', 'Adjusted Fisher Pearson', 'b', 'g', 'G', 'Minitab', 'Excel', 'SPSS', 'SAS', 'Stata')
    if (is.null(skew.type)) { skew.type <- 'Fisher Pearson' }
    if (!(skew.type %in% TYPES))                    { stop('Error: Input skew.type not recognised') } }
  if (MAXMOM >= 4) {
    if (is.null(kurt.type)) { kurt.type <- 'Fisher Pearson' }
    if (is.null(kurt.excess)) { kurt.excess <- FALSE }
    if (!(kurt.type %in% TYPES))                    { stop('Error: Input kurt.type not recognised') }
    if (!is.vector(kurt.excess))                    { stop('Error: Input kurt.excess should be a single logical value') }
    if (!is.logical(kurt.excess))                   { stop('Error: Input kurt.excess should be a single logical value') }
    if (length(kurt.excess) != 1)                   { stop('Error: Input kurt.excess should be a single logical value') } }

  #Check input include.sd
  if (!is.vector(include.sd))                       { stop('Error: Input include.sd should be a single logical value') }
  if (!is.logical(include.sd))                      { stop('Error: Input include.sd should be a single logical value') }
  if (length(include.sd) != 1)                      { stop('Error: Input include.sd should be a single logical value') }

  #Set skew and kurt adjustments
  #Default type with no adjustment is 'Fisher Pearson'
  if (MAXMOM >= 3) {
    skew.adj <- function(n) {
      A <- 1
      if (skew.type %in% c('Moment', 'b', 'Minitab')) {
        A <- ((n-1)/n)^(3/2) }
      if (skew.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
        A <- sqrt(n*(n-1))/(n-2) }
      A } }
  if (MAXMOM >= 4) {
    kurt.adj <- function(n) {
      B <- 1
      if (kurt.type %in% c('Moment', 'b', 'Minitab')) {
        B <- ((n-1)/n)^2 }
      if (kurt.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
        B <- (n+1)*(n-1)/((n-2)*(n-3)) }
      B }
    excess.adj <- function(n) {
      C <- -3*kurt.excess
      if (kurt.type %in% c('Adjusted Fisher Pearson', 'G', 'Excel', 'SPSS', 'SAS')) {
        C <- -3*kurt.excess*(n-1)^2/((n-2)*(n-3)) }
      C } }

  #Set output data frame when pooled is NULL
  #In this case all the groups are sample groups yet to be pooled
  if (is.null(pooled)) {

    #Set output data frame
    pool.n <- sum(n)
    OUT <- data.frame(n = c(n, pool.n))
    if (!is.null(names)) { rownames(OUT)[1:N] <- names }
    rownames(OUT)[N+1] <- '--pooled--'

    #Compute the pooled sample mean
    if (MAXMOM >= 1) {
      pool.mean <- sum(n*sample.mean)/pool.n
      deviation <- sample.mean - pool.mean
      OUT$sample.mean <- c(sample.mean, pool.mean) }

    #Compute the pooled sample variance
    if (MAXMOM >= 2) {
      SS        <- (n-1)*sample.var
      pool.SS   <- sum(SS) + sum(n*deviation^2)
      pool.var  <- pool.SS/(pool.n-1)
      if (include.sd) { OUT$sample.sd  <- c(sqrt(sample.var), sqrt(pool.var)) }
      OUT$sample.var <- c(sample.var, pool.var) }

    #Compute the pooled sample skewness
    if (MAXMOM >= 3) {
      SC        <- sample.skew*(SS^(3/2))/(skew.adj(n)*sqrt(n))
      pool.SC   <- sum(SC) + 3*sum(SS*deviation) + sum(n*deviation^3)
      pool.skew <- skew.adj(pool.n)*sqrt(pool.n)*pool.SC/pool.SS^(3/2)
      OUT$sample.skew <- c(sample.skew, pool.skew) }

    #Compute the pooled sample kurtosis
    if (MAXMOM >= 4) {
      SQ        <- (sample.kurt - excess.adj(n))*SS^2/(kurt.adj(n)*n)
      pool.SQ   <- sum(SQ) + 4*sum(SC*deviation) + 6*sum(SS*deviation^2) + sum(n*deviation^4)
      pool.kurt <- kurt.adj(pool.n)*pool.n*pool.SQ/pool.SS^2 + excess.adj(pool.n)
      OUT$sample.kurt <- c(sample.kurt, pool.kurt) } }

  #Set output data frame when pooled is specified
  #In this case the specified grouup is the pooled group and we are computing another group
  if (!is.null(pooled)) {

    #Set output data frame
    pool.n  <- n[pooled]
    part.n  <- sum(n[-pooled])
    other.n <- pool.n - part.n
    OUT <- data.frame(n = c(n, other.n))
    if (!is.null(names)) { rownames(OUT)[1:N] <- names }
    rownames(OUT)[pooled] <- '--pooled--'
    rownames(OUT)[N+1]    <- '--other--'

    #Compute the other sample mean
    if (MAXMOM >= 1) {
      pool.mean  <- sample.mean[pooled]
      part.mean  <- sum(n[-pooled]*sample.mean[-pooled])/part.n
      other.mean <- (pool.n*pool.mean - part.n*part.mean)/other.n
      deviation  <- sample.mean[-pooled] - part.mean
      OUT$sample.mean <- c(sample.mean, other.mean) }

    #Compute the other sample variance
    if (MAXMOM >= 2) {
      SS        <- (n[-pooled]-1)*sample.var[-pooled]
      part.SS   <- sum(SS) + sum(n[-pooled]*deviation^2)
      pool.SS   <- (pool.n-1)*sample.var[pooled]
      other.SS  <- pool.SS - part.SS - (part.n*pool.n/other.n)*(part.mean - pool.mean)^2
      other.var <- other.SS/(other.n-1)
      if (include.sd) { OUT$sample.sd  <- c(sqrt(sample.var), sqrt(other.var)) }
      OUT$sample.var <- c(sample.var, other.var) }

    #Compute the other sample skewness
    if (MAXMOM >= 3) {
      SC        <- sample.skew[-pooled]*SS^(3/2)/(skew.adj(n[-pooled])*sqrt(n[-pooled]))
      part.SC   <- sum(SC) + 3*sum(SS*deviation) + sum(n[-pooled]*deviation^3)
      pool.SC   <- sample.skew[pooled]*pool.SS^(3/2)/(skew.adj(pool.n)*sqrt(pool.n))
      other.SC  <- pool.SC - part.SC - (3*(pool.n*part.SS - part.n*pool.SS)/(other.n))*(part.mean - pool.mean) -
                   ((pool.n*part.n)*(pool.n + part.n)/(other.n)^2)*(part.mean - pool.mean)^3
      other.skew <- skew.adj(other.n)*sqrt(other.n)*other.SC/other.SS^(3/2)
      OUT$sample.skew <- c(sample.skew, other.skew) }

    #Compute the other sample kurtosis
    if (MAXMOM >= 4) {
      SQ        <- (sample.kurt[-pooled] - excess.adj(n[-pooled]))*SS^2/(kurt.adj(n[-pooled])*n[-pooled])
      part.SQ   <- sum(SQ) + 4*sum(SC*deviation) + 6*sum(SS*deviation^2) + sum(n[-pooled]*deviation^4)
      pool.SQ   <- (sample.kurt[pooled] - excess.adj(pool.n))*pool.SS^2/(kurt.adj(pool.n)*pool.n)
      other.SQ  <- pool.SQ - part.SQ - 4*((pool.n*other.SC - other.n*pool.SC)/(pool.n - other.n))*(other.mean - pool.mean) -
                   6*((pool.n^2*other.SS - other.n^2*pool.SS)/(pool.n - other.n)^2)*(other.mean - pool.mean)^2 -
                   ((other.n*pool.n^3 + other.n^2*pool.n^2 + other.n^3*pool.n)/(pool.n - other.n)^3)*(other.mean - pool.mean)^4
      other.kurt <- kurt.adj(other.n)*other.n*other.SQ/other.SS^2 + excess.adj(other.n)
      OUT$sample.kurt <- c(sample.kurt, other.kurt) }

    #Put pooled group at the end
    ORDER <- 1:(N+1)
    ORDER <- ORDER[-pooled]
    ORDER[N+1] <- pooled
    OUT <- OUT[ORDER, ] }

  #Add attributes
  attr(OUT, 'skew.type')   <- skew.type
  attr(OUT, 'kurt.type')   <- kurt.type
  attr(OUT, 'kurt.excess') <- kurt.excess

  #Return output
  OUT }
