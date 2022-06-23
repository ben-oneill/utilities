#' Kernel Density Estimator
#'
#' \code{KDE} returns the probability function for the kernel density estimator
#'
#' The kernel density estimator for a set of input data is obtained by taking a mixture distribution
#' consisting of a (possibly weighted) combination of kernels.  In this function we compute the KDE
#' using the kernel of the T-distribution; the function can also estimate a discretised version of the
#' KDE (taken over the integers) if required.  The degrees-of-freedom and the bandwidth for the KDE can
#' be specified in the inputs; if the bandwidth is not specified then it are estimated using the methods
#' set out in Sheather and Jones (1991) used in the \code{stats::density} function.  The output of the
#' function is a list of class \code{kde} that contains the probability functions for the KDE and
#' associated information.  The output object can be plotted to show the density function for the KDE.
#'
#' Note: The function has an option \code{to.environment} to allow the user to load the probability
#' functions to the global environment or another specified environment.  If this is set to \code{TRUE}
#' then the probability functions are loaded to the specified environment in addition to appearing as
#' elements of the output; there is a message informing the user if existing objects in the global
#' environment were overwritten.  If the functions are not loaded to the environment then the user
#' can use the function \code{KDE.load} to load them later from the produced object.
#'
#' @param data Input data for the kernel density estimator (a numeric vector)
#' @param weights Weights for the kernel density estimator (a numeric vector with the same length as the data)
#' @param bandwidth Bandwidth for the KDE; if \code{NULL} it is estimated
#' @param df Degrees-of-freedom for the T-distribution
#' @param density.name Name of the KDE distribution; used for naming of the probability functions (a character string)
#' @param value.name Name of the values in the data; used for naming the plot of the KDE
#' @param discrete Logical; if \code{TRUE} the function produces a discrete KDE over the integers
#' @param discrete.warn Logical; if \code{TRUE} the function gives a warning if non-discrete data is used to produce a discrete KDE
#' @param to.environment Logical; if \code{TRUE} the probability functions are attached to the global environment
#' @param envir The environment where the probability functions are loaded (if \code{to.environment} is \code{TRUE})
#' @return A \code{kde} object containing the probability functions for the kernel density estimator
#' @examples
#' k <- KDE(rnorm(500))
#' print(k)
#' plot(k)
#' KDE.load(k, environment()); ls()

KDE <- function (data, weights = NULL, bandwidth = NULL, df = Inf,
                 density.name = 'kde', value.name = 'Value',
                 discrete = FALSE, discrete.warn = TRUE,
                 to.environment = FALSE, envir = .GlobalEnv) {

  #Get preliminary call/name information
  CALL         <- sys.call()
  DATA.NAME    <- deparse(substitute(data))
  WEIGHTS.NAME <- deparse(substitute(weights))

  #Check input data and weights
  if (WEIGHTS.NAME == 'NULL') { WEIGHTS.NAME <- NA }
  if (!is.numeric(data))                  stop('Error: Input data should be numeric vector')
  if (length(data) == 0)                  stop('Error: Input data must not be empty')
  if (length(data) == 1)                  stop('Error: Input data must have at least two values to estimate the density')
  if (!is.null(weights)) {
    WEIGHTED <- TRUE
    if (!is.numeric(weights))             stop('Error: Input weights should be numeric')
    if (length(weights) != length(data))  stop('Error: Inputs data and weights should be the same length (unless the latter is NULL)')
    if (min(weights) < 0)                 stop('Error: Input weights cannot be negative')
    if (sum(weights) == 0)                stop('Error: Input weights cannot all be zero') } else {
    WEIGHTED <- FALSE }

  #Check bandwidth and df
  if (!is.null(bandwidth)) {
    BAND.EST <- FALSE
    if (!is.numeric(bandwidth))           stop('Error: Input bandwidth should be a numeric value')
    if (length(bandwidth) != 1)           stop('Error: Input bandwidth should be a single numeric value')
    if (min(bandwidth) <= 0)              stop('Error: Input bandwidth should be a positive value') } else {
    BAND.EST <- TRUE }
  if (!is.numeric(df))                    stop('Error: Input df should be a numeric value')
  if (length(df) != 1)                    stop('Error: Input df should be a single numeric value')
  if (min(df) <= 0)                       stop('Error: Input df should be a positive value')

  #Check discrete
  if (!is.logical(discrete))              stop('Error: Input discrete should be a logical value')
  if (length(discrete) != 1)              stop('Error: Input discrete should be a single logical value')
  if (!is.logical(discrete.warn))         stop('Error: Input discrete.warn should be a logical value')
  if (length(discrete.warn) != 1)         stop('Error: Input discrete.warn should be a single logical value')
  if ((discrete)&(discrete.warn)) {
    if (any(data != as.integer(data))) {
      warning('Non-integer data used to produce a discrete KDE over the integers') } }

  #Check other inputs
  if (!is.character(density.name))        stop('Error: Input density.name should be a character value')
  if (length(density.name) != 1)          stop('Error: Input density.name should be a single character value')
  if (!is.character(value.name))          stop('Error: Input value.name should be a character value')
  if (length(value.name) != 1)            stop('Error: Input value.name should be a single character value')
  if (!is.logical(to.environment))        stop('Error: Input to.environment should be a logical value')
  if (length(to.environment) != 1)        stop('Error: Input to.environment should be a single logical value')
  if (!is.environment(envir))             stop('Error: Input envir should be an environment')

  ###################################### Generate and return the output ######################################

  #Set the means and weights for KDE and estimate bandwidth
  DATA <- data
  k <- length(DATA)
  DF <- df
  if (is.null(weights))    { WEIGHTS <- NULL } else { WEIGHTS <- weights/sum(weights) }
  if (is.null(bandwidth))  { BAND <- stats::density(DATA, weights = WEIGHTS)$bw } else {
                             BAND <- bandwidth }

  #Create output list
  OUT <- vector(mode = 'list', length = 17)
  names(OUT)[1:4]  <- PROB.NAMES <- paste0(c('d', 'p', 'q', 'r'), density.name)
  names(OUT)[5:16] <- c('data.name', 'data', 'weighted', 'weights.name', 'weights',
                        'bandwidth', 'bandwidth.est', 'df', 'call', 'discrete',
                        'value.name', 'envir')
  OUT[[5]]   <- DATA.NAME
  OUT[[6]]   <- DATA
  OUT[[7]]   <- WEIGHTED
  OUT[[8]]   <- WEIGHTS.NAME
  if (WEIGHTED) { OUT[[9]] <- WEIGHTS }
  OUT[[10]]  <- BAND
  OUT[[11]]  <- BAND.EST
  OUT[[12]]  <- DF
  OUT[[13]]  <- deparse(CALL)
  OUT[[14]]  <- discrete
  OUT[[15]]  <- value.name
  OUT[[16]]  <- envir
  class(OUT) <- 'kde'

  #Generate probability functions for continuous KDE
  if (!discrete) {
    PROB.FUNCS <- .KDE.continuous(means = DATA, weights = WEIGHTS, bandwidth = BAND, df = DF)
    OUT[[1]] <- PROB.FUNCS[[1]]
    OUT[[2]] <- PROB.FUNCS[[2]]
    OUT[[3]] <- PROB.FUNCS[[3]]
    OUT[[4]] <- PROB.FUNCS[[4]] }

  #Generate probability functions for continuous KDE
  if (discrete) {
    PROB.FUNCS <- .KDE.discrete(means = DATA, weights = WEIGHTS, bandwidth = BAND, df = DF)
    OUT[[1]] <- PROB.FUNCS[[1]]
    OUT[[2]] <- PROB.FUNCS[[2]]
    OUT[[3]] <- PROB.FUNCS[[3]]
    OUT[[4]] <- PROB.FUNCS[[4]] }

  #Load functions to global environment (if required)
  if (to.environment) {

    #Check if functions already exist in the loading environment
    #Message user if functions are overwritten
    EXISTS <- rep(FALSE, 4)
    for (i in 1:4) { EXISTS[i] <- exists(PROB.NAMES[i], envir = envir) }
    if (any(EXISTS)) {
      if (identical(envir, .GlobalEnv)) {
        message('\n', '    Message from call: ', deparse(CALL), '\n',
                '    The following objects were overwritten in the global environment: \n \n',
                '    ', paste(PROB.NAMES[EXISTS], collapse = ', '), '\n')
      } else {
        message('\n', '    Message from call: ', deparse(CALL), '\n',
                '    The following objects were overwritten in the environment \'',
                deparse(substitute(envir)), '\' :', '\n \n',
                '    ', paste(PROB.NAMES[EXISTS], collapse = ', '), '\n') } }

    #Load functions to the global environment
    list2env(OUT[1:4], envir = envir) }

  #Return output
  OUT }

#' Utilities for KDE fits
#'
#' @name KDE_utils
NULL

#' \code{KDE.load} copies KDE distribution functions from a KDE object to a target environment.
#'
#' @rdname KDE_utils
#' @param object,x A KDE object
#' @param envir The target environment
#' @param overwrite If FALSE, aborts if the function names are already present in the target environment
#' @return KDE.load returns \code{envir}
KDE.load <- function (object, envir = NULL, overwrite = TRUE) {

  #Check object class
  if (!('kde' %in% class(object)))    stop('Error: This print method is only used for objects of class \'kde\'')

  #Check inputs envir and overwrite
  if (!is.null(envir)) {
    if (!is.environment(envir))             stop('Error: Input envir should be an environment') }
  if (!is.logical(overwrite))               stop('Error: Input overwrite should be a logical value')
  if (length(overwrite) != 1)               stop('Error: Input overwrite should be a single logical value')

  #Get information
  CALL       <- sys.call()
  PROB.FUNCS <- object[1:4]
  PROB.NAMES <- names(PROB.FUNCS)
  if (is.null(envir)) { envir <- object$envir }

  #Check if functions already exist in the loading environment

  EXISTS <- rep(FALSE, 4)
  for (i in 1:4) { EXISTS[i] <- exists(PROB.NAMES[i], envir = envir) }

  #Set functions to load to the loading environment
  if (overwrite) { LOAD <- PROB.FUNCS } else { LOAD <- PROB.FUNCS[!EXISTS] }

  #Message user if functions are/were already in the loading environment
  if (overwrite) {
    if (any(EXISTS)) {
      if (identical(envir, .GlobalEnv)) {
        message('\n', '    Message from call: ', deparse(CALL), '\n',
                '    The following objects were overwritten in the global environment: \n \n',
                '    ', paste(PROB.NAMES[EXISTS], collapse = ', '), '\n')
      } else {
        message('\n', '    Message from call: ', deparse(CALL), '\n',
                '    The following objects were overwritten in the environment \'',
                deparse(substitute(envir)), '\' :', '\n \n',
                '    ', paste(PROB.NAMES[EXISTS], collapse = ', '), '\n') } } }
  if (!overwrite) {
    if (any(EXISTS)) {
      if (identical(envir, .GlobalEnv)) {
        message('\n', '    Message from call: ', deparse(CALL), '\n',
                '    The following objects are already in the global environment (and were not overwritten): \n \n',
                '    ', paste(PROB.NAMES[EXISTS], collapse = ', '), '\n')
      } else {
        message('\n', '    Message from call: ', deparse(CALL), '\n',
                '    The following objects are already in the environment \'',
                deparse(substitute(envir)), '\' (and were not overwritten):', '\n \n',
                '    ', paste(PROB.NAMES[EXISTS], collapse = ', '), '\n') } } }

  #Load functions to the loading environment (if appropriate)
  invisible(list2env(LOAD, envir = envir)) }

#' \code{print.kde} prints the KDE object and returns it invisibly.
#'
#' @rdname KDE_utils
#' @param digits Number of digits to print
#' @return \code{print.kde} returns \code{x}, invisibly
print.kde <- function(x, digits = 6, ...) {

  #Check object class
  if (!inherits(x, 'kde'))            stop('Error: This print method is only used for objects of class \'kde\'')
  if (!is.numeric(digits))            stop('Error: Input digits should be a numeric value')
  if (length(digits) != 1)            stop('Error: Input digits should be a single numeric value')
  if (as.integer(digits) != digits)   stop('Error: Input digits should be an integer')
  if (min(digits) <= 1)               stop('Error: Input digits must be positive')

  #Extract information
  k  <- length(x$data)
  bw <- x$bandwidth
  df <- x$df
  nn <- names(x)[1:4]
  data.name    <- x$data.name
  weights.name <- x$weights.name
  weighted     <- x$weighted
  band.est     <- x$bandwidth.est
  discrete     <- x$discrete
  to.env       <- x$to.environment
  envir        <- x$envir

  #Print title
  if (!discrete) { cat('\n  Kernel Density Estimator (KDE) \n \n') }
  if (discrete)  { cat('\n  Discrete Kernel Density Estimator (KDE) \n \n') }

  #Print description
  cat('  Computed from', k, 'data points in the input', paste0('\'', data.name, '\'', '\n'))
  if (weighted) { cat('  KDE uses weights from input', paste0('\'', weights.name, '\''), '\n') }

  #Print bandwidth and df
  if (band.est) {
    cat('  Estimated bandwidth =', formatC(bw, digits = digits, format = 'f'), ' \n')
  } else {
    cat('  Input bandwidth =', formatC(bw, digits = digits, format = 'f'), ' \n') }
  cat('  Input degrees-of-freedom = ')
  if (df == Inf) {
    cat('Inf \n \n')
  } else {
    cat(formatC(df, digits = digits, format = 'f'), ' \n \n') }

  #Print information on probability functions
  STAR <- rep(FALSE, 4)
  if (exists(nn[1])) { if (identical(get(nn[1], envir = envir), x[[1]])) { STAR[1] <- TRUE } }
  if (exists(nn[2])) { if (identical(get(nn[2], envir = envir), x[[2]])) { STAR[2] <- TRUE } }
  if (exists(nn[3])) { if (identical(get(nn[3], envir = envir), x[[3]])) { STAR[3] <- TRUE } }
  if (exists(nn[4])) { if (identical(get(nn[4], envir = envir), x[[4]])) { STAR[4] <- TRUE } }
  if (!discrete) { cat('  Probability functions for the KDE are the following: \n \n') }
  if (discrete)  { cat('  Probability functions for the discrete KDE are the following: \n \n') }
  cat('      Density function:                  ', nn[1], ifelse(STAR[1], '*', ''), '\n',
       '     Distribution function:             ', nn[2], ifelse(STAR[2], '*', ''), '\n',
       '     Quantile function:                 ', nn[3], ifelse(STAR[3], '*', ''), '\n',
       '     Random generation function:        ', nn[4], ifelse(STAR[4], '*', ''), '\n', '\n')
  if (any(STAR)) {
    if (identical(envir, .GlobalEnv)) {
      cat('  * This function is presently loaded in the global environment \n \n')
    } else {
    cat('  * This function is presently loaded in the environment \'',
        deparse(substitute(envir)), '\' \n \n') } }

  invisible(x)
}

#' \code{plot.kde} draws a plot of the KDE.
#' @rdname KDE_utils
#' @param n number of bins
#' @param cut cutoffs for xaxis (in steps of bw)
#' @param fill.colour,fill.color  fill color of bars
#' @param ... unused
#' @returns \code{plot.kde} returns the plot as recorded by \code{recordPlot}
plot.kde <- function(x, digits = 6, n = 512, cut = 4,
                     fill.colour = 'dodgerblue', fill.color = fill.colour, ...) {

  #Check object class
  if (!('kde' %in% class(x)))    stop('Error: This plot method is only used for objects of class \'kde\'')

  #Check inputs
  if (!is.numeric(digits))            stop('Error: Input digits should be a numeric value')
  if (length(digits) != 1)            stop('Error: Input digits should be a single numeric value')
  if (as.integer(digits) != digits)   stop('Error: Input digits should be an integer')
  if (min(digits) <= 1)               stop('Error: Input digits must be positive')
  if (!is.numeric(n))                 stop('Error: Input n should be a numeric value')
  if (length(n) != 1)                 stop('Error: Input n should be a single numeric value')
  if (min(n) < 64)                    stop('Error: Input n is too small')
  if (!is.numeric(cut))               stop('Error: Input cut should be a numeric value')
  if (length(cut) != 1)               stop('Error: Input cut should be a single numeric value')
  if (min(cut) <= 0)                  stop('Error: Input cut must be positive')

  #Check fill.colour/color
  if (length(fill.colour) != 1)       stop('Error: Input fill.colour should be a single colour')
  if (length(fill.color) != 1)        stop('Error: Input fill.color should be a single colour')
  if (!fill.color %in% colours())     stop('Error: Input fill.colour/color must be in \'colours()\'')

  #Extract information
  df   <- x$df
  bw   <- x$bandwidth
  data <- x$data
  dens <- x[[1]]
  band.est   <- x$bandwidth.est
  discrete   <- x$discrete
  value.name <- x$value.name

  #Compute value range for the plot
  xmin <- min(data) - bw*cut
  xmax <- max(data) + bw*cut

  #Create the subtitle
  E1 <- 'Input degrees-of-freedom = '
  if (band.est) { E2 <- 'Estimated bandwidth = '          } else { E2 <- 'Input bandwidth = '          }
  SUBTITLE <- paste0( '(', E1, if (df == Inf) { '\U221E' } else { formatC(df, digits = digits, format = 'f') },
                     ', ', E2, formatC(bw, digits = digits, format = 'f'), ')')

  #Create the plot (continuous KDE)
  if (!discrete) {
    xx   <- xmin + (xmax-xmin)*((1:n)-1)/(n-1)
    yy   <- dens(xx)
    zz   <- rep(0, n)
    plot(xx, yy, type = 'l',
         xlab = value.name, ylab = 'Density')
    title(main = 'Kernel Density Estimator')
    mtext(side = 3, line = 0.25, cex = 0.8, SUBTITLE)
    polygon(c(xx, rev(xx)), c(zz, rev(yy)), col = fill.color) }

  #Create the plot (discrete KDE)
  if (discrete) {
    xx   <- (floor(xmin)):(ceiling(xmax))
    yy   <- dens(xx)
    zz   <- rep(0, n)
    barplot(names.arg = xx, height = yy,
            ylim = c(0, 1.1*max(yy)), col = fill.color,
            xlab = value.name, ylab = 'Probability')
    title(main = 'Discrete Kernel Density Estimator')
    mtext(side = 3, line = 0.25, cex = 0.8, SUBTITLE) }

  #Save the plot
  PLOT <- recordPlot()
  dev.off()

  #Return the plot
  PLOT }


.KDE.continuous <- function(means, weights = NULL, bandwidth, df) {

  #Set output list
  PROB.FUNCS <- vector(length = 4, mode = 'list')
  names(PROB.FUNCS) <- c('dkde', 'pkde', 'qkde', 'rkde')

  #Generate density function
  PROB.FUNCS[[1]] <- function(x, bandwidth, df, log = FALSE) {

    #Set mean-length and weights
    k <- length(means)
    weighted <- !is.null(weights)
    if (!weighted) { weights <- rep(1/k, k) }

    #Compute log-density values
    n <- length(x)
    LOGDENS   <- rep(0, n)
    LOGMATRIX <- matrix(0, nrow = n, ncol = k)
    for (i in 1:n) {
      LOGMATRIX[i,] <- dt((x[i]-means)/bandwidth, df = df, log = TRUE) - log(bandwidth) + log(weights)
      LOGDENS[i]    <- matrixStats::logSumExp(LOGMATRIX[i,]) }

    #Return output
    if (log) { LOGDENS } else { exp(LOGDENS) } }

  #Generate cumulative distribution function
  PROB.FUNCS[[2]] <- function(x, bandwidth, df, lower.tail = TRUE, log.p = FALSE) {

    #Set mean-length and weights
    k <- length(means)
    weighted <- !is.null(weights)
    if (!weighted) { weights <- rep(1/k, k) }

    #Compute log-density values
    n <- length(x)
    LOGCDF    <- rep(0, n)
    LOGMATRIX <- matrix(0, nrow = n, ncol = k)
    for (i in 1:n) {
      LOGMATRIX[i,] <- pt((x[i]-means)/bandwidth, df = df, lower.tail = lower.tail, log.p = TRUE) + log(weights)
      LOGCDF[i]     <- matrixStats::logSumExp(LOGMATRIX[i,]) }

    #Return output
    if (log.p) { LOGCDF } else { exp(LOGCDF) } }

  #Generate quantile function
  PROB.FUNCS[[3]] <- function(p, bandwidth, df, lower.tail = TRUE, log.p = FALSE) {

    #Set mean-length and weights
    k <- length(means)
    weighted <- !is.null(weights)
    if (!weighted) { weights <- rep(1/k, k) }

    #Compute approximate quantiles
    n <- length(p)
    if (log.p) { LOGP <- p } else { LOGP <- log(p) }
    if (!lower.tail) { LOGP <- VGAM::log1mexp(-LOGP) }
    QUANTILES.APPROX <- quantile(x = means, probs = exp(LOGP), names = FALSE)

    #Refine approximate quantiles via optimisation
    QUANTILES <- QUANTILES.APPROX
    for (i in 1:n) {

      #Check for special cases
      if (LOGP[i] == -Inf) { QUANTILES[i] <- -Inf }
      if (LOGP[i] == 0)    { QUANTILES[i] <- Inf }

      #Optimise for remaining cases
      if ((LOGP[i] != -Inf)&(LOGP[i] != 0)) {

        #Set objective function
        OBJECTIVE <- function(x) {
          LOGVALS <- pt((x-means)/bandwidth, df = df, lower.tail = lower.tail, log.p = TRUE) + log(weights)
          LOGCDF  <- matrixStats::logSumExp(LOGVALS)
          (LOGCDF - LOGP[i])^2 }

        #Refine quantile
        NLM <- nlm(f = OBJECTIVE, p = QUANTILES.APPROX[i])
        QUANTILES[i] <- NLM$estimate } }

    #Return output
    QUANTILES }

  #Generate random-generation function
  PROB.FUNCS[[4]] <- function(n, bandwidth, df) {

    #Set mean-length and weights
    k <- length(means)
    weighted <- !is.null(weights)
    if (!weighted) { weights <- rep(1/k, k) }

    #Generate random variables
    SAMPLE <- sample.int(k, size = n, replace = TRUE, prob = weights)
    MMM <- means[SAMPLE]
    TTT <- rt(n, df = df)
    OUT <- MMM + bandwidth*TTT

    #Return output
    OUT }

  #Link means and weights to environments of probability functions
  environment(PROB.FUNCS[[1]]) <- list2env(list(means = means, weights = weights),
                                           parent = environment(PROB.FUNCS[[1]]))
  environment(PROB.FUNCS[[2]]) <- list2env(list(means = means, weights = weights),
                                           parent = environment(PROB.FUNCS[[2]]))
  environment(PROB.FUNCS[[3]]) <- list2env(list(means = means, weights = weights),
                                           parent = environment(PROB.FUNCS[[3]]))
  environment(PROB.FUNCS[[4]]) <- list2env(list(means = means, weights = weights),
                                           parent = environment(PROB.FUNCS[[4]]))

  #Set arguments for probability functions
  formals(PROB.FUNCS[[1]])$x <- substitute()
  formals(PROB.FUNCS[[2]])$x <- substitute()
  formals(PROB.FUNCS[[3]])$p <- substitute()
  formals(PROB.FUNCS[[4]])$n <- substitute()
  formals(PROB.FUNCS[[1]])$df <- df
  formals(PROB.FUNCS[[2]])$df <- df
  formals(PROB.FUNCS[[3]])$df <- df
  formals(PROB.FUNCS[[4]])$df <- df
  formals(PROB.FUNCS[[1]])$bandwidth <- bandwidth
  formals(PROB.FUNCS[[2]])$bandwidth <- bandwidth
  formals(PROB.FUNCS[[3]])$bandwidth <- bandwidth
  formals(PROB.FUNCS[[4]])$bandwidth <- bandwidth

  #Return output
  PROB.FUNCS }


.KDE.discrete <- function (means, weights = NULL, bandwidth, df) {

  #Set output list
  PROB.FUNCS <- vector(length = 4, mode = 'list')
  names(PROB.FUNCS) <- c('dkde', 'pkde', 'qkde', 'rkde')

  #Generate density function
  PROB.FUNCS[[1]] <- function(x, bandwidth, df, log = FALSE) {

    #Set mean-length and weights
    k <- length(means)
    weighted <- !is.null(weights)
    if (!weighted) { weights <- rep(1/k, k) }

    #Compute log-density values
    n <- length(x)
    IMAG <- complex(imaginary = 1)
    LOGDENS   <- rep(-Inf, n)
    LOGMATRIX <- matrix(0, nrow = n, ncol = k)
    for (i in 1:n) {
      if (x[i] == as.integer(x[i])) {
      LOGMATRIX[i,] <- pt((x[i]-means+0.5)/bandwidth, df = df, log.p = TRUE) + log(weights) +
                        IMAG*(pt((x[i]-means-0.5)/bandwidth, df = df, log.p = TRUE) + log(weights))
                        L1 <- matrixStats::logSumExp(Re(LOGMATRIX[i,]))
                        L2 <- matrixStats::logSumExp(Im(LOGMATRIX[i,]))
                        LOGDENS[i] <- L1 + VGAM::log1mexp(L1-L2) } }

    #Return output
    if (log) { LOGDENS } else { exp(LOGDENS) } }

  #Generate cumulative distribution function
  PROB.FUNCS[[2]] <- function(x, bandwidth, df, lower.tail = TRUE, log.p = FALSE) {

    #Set mean-length and weights
    k <- length(means)
    weighted <- !is.null(weights)
    if (!weighted) { weights <- rep(1/k, k) }

    #Compute log-density values
    n <- length(x)
    LOGCDF    <- rep(0, n)
    LOGMATRIX <- matrix(0, nrow = n, ncol = k)
    for (i in 1:n) {
      LOGMATRIX[i,] <- pt((ceiling(x[i]-0.5)-means+0.5)/bandwidth, df = df, lower.tail = lower.tail, log.p = TRUE) + log(weights)
      LOGCDF[i]     <- matrixStats::logSumExp(LOGMATRIX[i,]) }

    #Return output
    if (log.p) { LOGCDF } else { exp(LOGCDF) } }

  #Generate quantile function
  PROB.FUNCS[[3]] <- function(p, bandwidth, df, lower.tail = TRUE, log.p = FALSE) {

    #Set mean-length and weights
    k <- length(means)
    weighted <- !is.null(weights)
    if (!weighted) { weights <- rep(1/k, k) }

    #Compute approximate quantiles
    n <- length(p)
    if (log.p) { LOGP <- p } else { LOGP <- log(p) }
    if (!lower.tail) { LOGP <- VGAM::log1mexp(-LOGP) }
    QUANTILES.APPROX <- quantile(x = means, probs = exp(LOGP), names = FALSE)

    #Refine approximate quantiles via optimisation
    QUANTILES <- QUANTILES.APPROX
    for (i in 1:n) {

      #Check for special cases
      if (LOGP[i] == -Inf) { QUANTILES[i] <- -Inf }
      if (LOGP[i] == 0)    { QUANTILES[i] <- Inf }

      #Optimise for remaining cases
      if ((LOGP[i] != -Inf)&(LOGP[i] != 0)) {

        #Set objective function
        OBJECTIVE <- function(x) {
          LOGVALS <- pt((x-means)/bandwidth, df = df, lower.tail = lower.tail, log.p = TRUE) + log(weights)
          LOGCDF  <- matrixStats::logSumExp(LOGVALS)
          (LOGCDF - LOGP[i])^2 }

        #Refine quantile
        NLM <- nlm(f = OBJECTIVE, p = QUANTILES.APPROX[i])
        QUANTILES[i] <- NLM$estimate } }

    #Return output
    QUANTILES <- ceiling(QUANTILES - 0.5)
    QUANTILES }

  #Generate random-generation function
  PROB.FUNCS[[4]] <- function(n, bandwidth, df) {

    #Set mean-length and weights
    k <- length(means)
    weighted <- !is.null(weights)
    if (!weighted) { weights <- rep(1/k, k) }

    #Generate random variables
    SAMPLE <- sample.int(k, size = n, replace = TRUE, prob = weights)
    MMM <- means[SAMPLE]
    TTT <- rt(n, df = df)
    OUT <- ceiling(MMM + bandwidth*TTT - 0.5)

    #Return output
    OUT }

  #Link means and weights to environments of probability functions
  environment(PROB.FUNCS[[1]]) <- list2env(list(means = means, weights = weights),
                                           parent = environment(PROB.FUNCS[[1]]))
  environment(PROB.FUNCS[[2]]) <- list2env(list(means = means, weights = weights),
                                           parent = environment(PROB.FUNCS[[2]]))
  environment(PROB.FUNCS[[3]]) <- list2env(list(means = means, weights = weights),
                                           parent = environment(PROB.FUNCS[[3]]))
  environment(PROB.FUNCS[[4]]) <- list2env(list(means = means, weights = weights),
                                           parent = environment(PROB.FUNCS[[4]]))

  #Set arguments for probability functions
  formals(PROB.FUNCS[[1]])$x <- substitute()
  formals(PROB.FUNCS[[2]])$x <- substitute()
  formals(PROB.FUNCS[[3]])$p <- substitute()
  formals(PROB.FUNCS[[4]])$n <- substitute()
  formals(PROB.FUNCS[[1]])$df <- df
  formals(PROB.FUNCS[[2]])$df <- df
  formals(PROB.FUNCS[[3]])$df <- df
  formals(PROB.FUNCS[[4]])$df <- df
  formals(PROB.FUNCS[[1]])$bandwidth <- bandwidth
  formals(PROB.FUNCS[[2]])$bandwidth <- bandwidth
  formals(PROB.FUNCS[[3]])$bandwidth <- bandwidth
  formals(PROB.FUNCS[[4]])$bandwidth <- bandwidth

  #Return output
  PROB.FUNCS }
