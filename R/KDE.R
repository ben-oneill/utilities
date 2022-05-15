#' Kernel Density Estimator
#'
#' \code{KDE} returns the probability function for the kernel density estimator
#'
#' The kernel density estimator for a set of input data is obtained by taking a mixture distribution
#' consisting of a (possibly weighted) combination of kernels.  In this function we compute the KDE
#' using the kernel of the T-distribution.  The degrees-of-freedom and the bandwidth for the KDE can
#' be specified in the inputs; if they are not specified then they are estimated using leave-one-out
#' maximum-likelihood estimation (LOO-MLE).  The output of the function is a list of class \code{kde}
#' that contains the probability functions for the KDE and associated information.  The output object
#' can be plotted to show the density function for the KDE.
#'
#' Note: The function has an option \code{to.environment} to allow the user to load the probability
#' functions to the global environment.  If this is set to \code{TRUE} then the probability functions
#' are loaded to the global environment in addition to appearing as elements of the output; there is
#' a message informing the user if existing objects in the global environment were overwritten.
#'
#' Note: For simplicity, the function uses a value \code{df.norm} which is the maximum value of the
#' degrees-of-freedom that is recognised as differing from normality.  If the inputted or estimated
#' degrees-of-freedom is higher than this value then the algorithm will round the value to infinity
#' (i.e., using the normal-distribution kernel instead of the T-distribution kernel).
#'
#' @usage \code{KDE}
#' @param data Input data for the kernel density estimator (a numeric vector)
#' @param na.rm Logical; if \code{TRUE} the function will remove \code{NA} values from the data
#' @param weights Weights for the kernel density estimator (a numeric vector with the same length as the data)
#' @param bandwidth Bandwidth for the KDE; if \code{NULL} it is estimated using LOO-MLE
#' @param df Degrees-of-freedom for the T-distribution; if \code{NULL} it is estimated using LOO-MLE
#' @param df.norm Maximum recognised value for the degrees-of-freedom; rounded up to Inf if above this value
#' @param density.name Name of the KDE distribution; used for naming of the probability functions (a character string)
#' @param value.name Name of the values in the data; used for naming the plot of the KDE
#' @param to.environment Logical; if \code{TRUE} the probability functions are attached to the global environment
#' @return A \code{kde} object containing the probability functions for the kernel density estimator

KDE <- function (data, na.rm = FALSE, weights = NULL,
                 bandwidth = NULL, df = Inf, df.norm = 1000,
                 density.name = 'kde', value.name = 'Value', to.environment = FALSE) {

  #Get preliminary call/name information
  CALL         <- sys.call()
  DATA.NAME    <- deparse(substitute(data))
  WEIGHTS.NAME <- deparse(substitute(weights))

  #Check input data and weights
  if (WEIGHTS.NAME == 'NULL') { WEIGHTS.NAME <- NA }
  if (!is.numeric(data))                  stop('Error: Input data should be numeric vector')
  if (length(data) == 0)                  stop('Error: Input data must not be empty')
  if (length(data) == 1)                  stop('Error: Input data must have at least two values to estimate the density')
  if (!is.logical(na.rm))                 stop('Error: Input na.rm should be a logical value')
  if (length(na.rm) != 1)                 stop('Error: Input na.rm should be a single logical value')
  if (!is.null(weights)) {
    WEIGHTED <- TRUE
    if (!is.numeric(weights))             stop('Error: Input weights should be numeric')
    if (length(weights) != length(data))  stop('Error: Inputs data and weights should be the same length (unless the latter is NULL)')
    if (min(weights) < 0)                 stop('Error: Input weights cannot be negative') } else {
    WEIGHTED <- FALSE }

  #Check bandwidth and df
  if (!is.null(bandwidth)) {
    BAND.EST <- FALSE
    if (!is.numeric(bandwidth))           stop('Error: Input bandwidth should be a numeric value')
    if (length(bandwidth) != 1)           stop('Error: Input bandwidth should be a single numeric value')
    if (min(bandwidth) <= 0)              stop('Error: Input bandwidth should be a positive value') } else {
    BAND.EST <- TRUE }
  if (!is.null(df)) {
    DF.EST <- FALSE
    if (!is.numeric(df))                  stop('Error: Input df should be a numeric value')
    if (length(df) != 1)                  stop('Error: Input df should be a single numeric value')
    if (min(df) <= 0)                     stop('Error: Input df should be a positive value') } else {
    DF.EST <- TRUE }
  if (!is.numeric(df.norm))               stop('Error: Input df.norm should be numeric value')
  if (length(df.norm) != 1)               stop('Error: Input df.norm should be a single numeric value')
  if (min(df.norm) < 0)                   stop('Error: Input df.norm cannot be negative')

  #Check other inputs
  if (!is.character(density.name))        stop('Error: Input density.name should be a character value')
  if (length(density.name) != 1)          stop('Error: Input density.name should be a single character value')
  if (!is.character(value.name))          stop('Error: Input value.name should be a character value')
  if (length(value.name) != 1)            stop('Error: Input value.name should be a single character value')
  if (!is.logical(to.environment))        stop('Error: Input to.environment should be a logical value')
  if (length(to.environment) != 1)        stop('Error: Input to.environment should be a single logical value')

  ############################################################################################################
  ########################################## Set the data and weights ########################################
  ############################################################################################################

  #Finalise the data and weights
  if (na.rm) {
    DATA <- data[!is.na(data)]
    WEIGHTS <- WEIGHTS[!is.na(data)] } else {
      DATA <- data }
  if (any(is.na(DATA))) {
    stop('Error: NA values in the data; set na.rm = TRUE if you want to remove these') }
  k <- length(DATA)
  if (!WEIGHTED) { WEIGHTS <- rep(1/k, k) } else { WEIGHTS <- weights/sum(weights) }

  ############################################################################################################
  ######################################## Estimate the KDE parameters #######################################
  ############################################################################################################

  #Estimate the parameters via LOO-MLE (if not given)
  #Case where df is given and bandwidth is to be estimated
  if ((BAND.EST)&(!DF.EST)) {
    NEGLOGLIKE <- function(p) {
      VEC <- rep(0, k)
      for (i in 1:k) {
        DDD <- dt((DATA[-i]-DATA[i])*exp(-p), df = df, log = TRUE) - p + log(WEIGHTS[-i])
        VEC[i] <- matrixStats::logSumExp(DDD) }
      -sum(VEC) }
    NLM  <- nlm(f = NEGLOGLIKE, p = 0)
    BAND <- exp(NLM$estimate)
    DF   <- df }

  #Case where df is to be estimated and bandwidth is given
  if ((!BAND.EST)&(DF.EST)) {
    NEGLOGLIKE <- function(p) {
      VEC <- rep(0, k)
      for (i in 1:k) {
        DDD <- dt((DATA[-i]-DATA[i])/bandwidth, df = exp(p), log = TRUE) + log(WEIGHTS[-i])
        VEC[i] <- matrixStats::logSumExp(DDD) }
      -sum(VEC) }
    NLM  <- nlm(f = NEGLOGLIKE, p = 0)
    BAND <- bandwidth
    DF   <- exp(NLM$estimate)
    if (DF > df.norm) { DF <- Inf } }

  #Case where df and bandwidth are both to be estimated
  if ((BAND.EST)&(DF.EST)) {
    NEGLOGLIKE <- function(p) {
      VEC <- rep(0, k)
      for (i in 1:k) {
        DDD <- dt((DATA[-i]-DATA[i])*exp(-p[1]), df = exp(p[2]), log = TRUE) - p[1] + log(WEIGHTS[-i])
        VEC[i] <- matrixStats::logSumExp(DDD) }
      -sum(VEC) }
    NLM  <- nlm(f = NEGLOGLIKE, p = c(0,0))
    BAND <- exp(NLM$estimate[1])
    DF   <- exp(NLM$estimate[2])
    if (DF > df.norm) { DF <- Inf } }

  ############################################################################################################
  ###################################### Generate and return the output ######################################
  ############################################################################################################

  #Create output list
  OUT <- vector(mode = 'list', length = 16)
  names(OUT)[1:4]  <- PROB.NAMES <- paste0(c('d', 'p', 'q', 'r'), density.name)
  names(OUT)[5:16] <- c('data.name', 'data', 'weighted', 'weights.name', 'weights',
                        'bandwidth', 'bandwidth.est', 'df', 'df.est', 'call', 'to.environment', 'value.name')
  OUT[[5]]   <- DATA.NAME
  OUT[[6]]   <- DATA
  OUT[[7]]   <- WEIGHTED
  OUT[[8]]   <- WEIGHTS.NAME
  OUT[[9]]   <- WEIGHTS
  OUT[[10]]  <- BAND
  OUT[[11]]  <- BAND.EST
  OUT[[12]]  <- DF
  OUT[[13]]  <- DF.EST
  OUT[[14]]  <- deparse(CALL)
  OUT[[15]]  <- to.environment
  OUT[[16]]  <- value.name
  class(OUT) <- 'kde'

  #Create computational components for KDE
  GEN.K       <- paste('k <-', length(DATA))
  GEN.DATA    <- paste('means <-', paste(deparse(DATA), collapse = ''))
  if (!WEIGHTED) {
    GEN.WEIGHTS <- NULL } else {
      GEN.WEIGHTS <- paste('weights <-', paste(deparse(WEIGHTS), collapse = '')) }

  #Create function-generator
  generate_function <- function(..., commands, envir = parent.frame(),
                                openbracket = FALSE, closebracket = FALSE) {
    body <- paste(commands, collapse = '\n')
    if (!openbracket)  { body <- paste(c('{', body), collapse = '\n') }
    if (!closebracket) { body <- paste(c(body, '}'), collapse = '\n') }
    as.function(c(..., str2lang(body)), envir = envir) }

  #Generate density function
  COMM <- c('',
            '#Set data and weights',
            GEN.K, GEN.DATA, GEN.WEIGHTS,
            '',
            '#Compute log-density values',
            'n <- length(x)',
            'LOGDENS <- rep(0, n)',
            'LOGMATRIX <- matrix(0, nrow = n, ncol = k)',
            if (!WEIGHTED) {
              c('  for (i in 1:n) {',
                '    LOGMATRIX[i,] <- dt((x[i]-means)/bandwidth, df = df, log = TRUE) - log(bandwidth) - log(k) ',
                '    LOGDENS[i] <- matrixStats::logSumExp(LOGMATRIX[i,]) }') } else {
                  c('  for (i in 1:n) {',
                    '    LOGMATRIX[i,] <- dt((x[i]-means)/bandwidth, df = df, log = TRUE) - log(bandwidth) + log(weights) ',
                    '    LOGDENS[i] <- matrixStats::logSumExp(LOGMATRIX[i,]) }') },
            '',
            '#Return output',
            'if (log) { LOGDENS } else { exp(LOGDENS) }')
  OUT[[1]] <- generate_function(x = NA, bandwidth = BAND, df = DF, log = FALSE, commands = COMM)

  #Generate cumulative distribution function
  COMM <- c('',
            '#Set data and weights',
            GEN.K, GEN.DATA, GEN.WEIGHTS,
            '',
            '#Compute log-CDF values',
            'n <- length(x)',
            'LOGCDF <- rep(0, n)',
            'LOGMATRIX <- matrix(0, nrow = n, ncol = k)',
            if (!WEIGHTED) {
              c('  for (i in 1:n) {',
                '    LOGMATRIX[i,] <- pt((x[i]-means)/bandwidth, df = df, lower.tail = lower.tail, log.p = TRUE) - log(k) ',
                '    LOGCDF[i] <- matrixStats::logSumExp(LOGMATRIX[i,]) }') } else {
                  c('  for (i in 1:n) {',
                    '    LOGMATRIX[i,] <- pt((x[i]-means)/bandwidth, df = df, lower.tail = lower.tail, log.p = TRUE) + log(weights) ',
                    '    LOGCDF[i] <- matrixStats::logSumExp(LOGMATRIX[i,]) }') },
            '',
            '#Return output',
            'if (log.p) { LOGCDF } else { exp(LOGCDF) }')
  OUT[[2]] <- generate_function(x = NA, bandwidth = BAND, df = DF, lower.tail = TRUE, log.p = FALSE, commands = COMM)

  #Generate quantile function
  COMM <- c('',
            '#Set data and weights',
            GEN.K, GEN.DATA, GEN.WEIGHTS,
            '',
            '#Get log-probabilities from inputs',
            'n <- length(p)',
            'if (log.p) { LOGP <- p } else { LOGP <- log(p) }',
            'if (!lower.tail) { LOGP <- VGAM::log1mexp(-LOGP) }',
            '',
            '#Compute approximate quantiles',
            'QUANTILES.APPROX <- quantile(x = means, probs = exp(LOGP), names = FALSE)',
            '',
            '#Refine approximate quantiles via optimisation',
            'QUANTILES <- QUANTILES.APPROX',
            '',
            'for (i in 1:n) {',
            '  ',
            '  #Check for special cases',
            '  if (LOGP[i] == -Inf) { QUANTILES[i] <- -Inf }',
            '  if (LOGP[i] == 0)    { QUANTILES[i] <- Inf }',
            '  ',
            '  #Optimise for remaining cases',
            '  if ((LOGP[i] != -Inf)&(LOGP[i] != 0)) {',
            '    ',
            '    #Set objective function',
            if (!WEIGHTED) { c(
              '    OBJECTIVE <- function(x) {',
              '      LOGVALS <- pt((x-means)/bandwidth, df = df, lower.tail = lower.tail, log.p = TRUE) - log(k) ',
              '      LOGCDF  <- matrixStats::logSumExp(LOGVALS) ',
              '      (LOGCDF - LOGP[i])^2 }') } else { c(
                '    OBJECTIVE <- function(x) {',
                '      LOGVALS <- pt((x-means)/bandwidth, df = df, lower.tail = lower.tail, log.p = TRUE) + log(weights) ',
                '      LOGCDF  <- matrixStats::logSumExp(LOGVALS) ',
                '      (LOGCDF - LOGP[i])^2 }') },
            '    ',
            '    #Refine quantile',
            '    NLM <- nlm(f = OBJECTIVE, p = QUANTILES.APPROX[i])',
            '    QUANTILES[i] <- NLM$estimate } }',
            '',
            '#Return output',
            'QUANTILES')
  OUT[[3]] <- generate_function(p = NA, bandwidth = BAND, df = DF, lower.tail = TRUE, log.p = FALSE, commands = COMM)

  #Generate random-generation function
  COMM <- c('',
            '#Set data and weights',
            GEN.K, GEN.DATA, GEN.WEIGHTS,
            '',
            '#Generate random variables',
            if (!WEIGHTED) {
              'SAMPLE <- sample.int(k, size = n, replace = TRUE)' } else {
                'SAMPLE <- sample.int(k, size = n, replace = TRUE, prob = weights)' },
            'MMM <- means[SAMPLE]',
            'TTT <- rt(n, df = df)',
            'OUT <- MMM + bandwidth*TTT',
            '',
            '#Return output',
            'OUT')
  OUT[[4]] <- generate_function(n = NA, bandwidth = BAND, df = DF, commands = COMM)

  #Load functions to global environment (if required)
  if (to.environment) {

    #Check if functions already exist in the global environment
    #Message user if functions are overwritten
    EXISTS <- rep(FALSE, 4)
    for (i in 1:4) { EXISTS[i] <- exists(PROB.NAMES[i], envir = .GlobalEnv) }
    if (any(EXISTS)) {
      message('\n', '    Message from call: ', deparse(CALL), '\n',
              '    The following objects were overwritten in the global environment:',
              '\n \n', '    ', paste(PROB.NAMES[EXISTS], collapse = ', '), '\n') }

    #Load functions to the global environment
    list2env(OUT[1:4], envir = .GlobalEnv) }

  #Return output
  OUT }


print.kde <- function(object, digits = 6) {

  #Check object class
  if (!('kde' %in% class(object)))    stop('Error: This print method is only used for objects of class \'kde\'')
  if (!is.numeric(digits))            stop('Error: Input digits should be a numeric value')
  if (length(digits) != 1)            stop('Error: Input digits should be a single numeric value')
  if (as.integer(digits) != digits)   stop('Error: Input digits should be an integer')
  if (min(digits) <= 1)               stop('Error: Input digits must be positive')

  #Extract information
  k  <- length(object$data)
  bw <- object$bandwidth
  df <- object$df
  nn <- names(object)[1:4]
  data.name    <- object$data.name
  weights.name <- object$weights.name
  weighted     <- object$weighted
  band.est     <- object$bandwidth.est
  df.est       <- object$df.est
  to.env       <- object$to.environment

  #Print title and description
  cat('\n  Kernel Density Estimator (KDE) \n \n')
  cat('  Computed from', k, 'data points in the input', paste0('\'', data.name, '\'', '\n'))
  if (weighted) { cat('  KDE uses weights from input', paste0('\'', weights.name, '\''), '\n') }

  #Print bandwidth and df
  if (band.est) {
    cat('  Estimated bandwidth =', formatC(bw, digits = digits, format = 'f'), ' \n')
  } else {
    cat('  Input bandwidth =', formatC(bw, digits = digits, format = 'f'), ' \n') }
  if (df.est) {
    cat('  Estimated degrees-of-freedom = ')
    if (df == Inf) { cat('Inf \n \n') } else { cat(formatC(df, digits = digits, format = 'f'), ' \n \n') }
  } else {
    cat('  Input degrees-of-freedom = ')
    if (df == Inf) { cat('Inf \n \n') } else { cat(formatC(df, digits = digits, format = 'f'), ' \n \n') } }

  #Print information on probability functions
  STAR <- rep(FALSE, 4)
  if (exists(nn[1])) { if (identical(get(nn[1], envir = .GlobalEnv), object[[1]])) { STAR[1] <- TRUE } }
  if (exists(nn[2])) { if (identical(get(nn[2], envir = .GlobalEnv), object[[2]])) { STAR[2] <- TRUE } }
  if (exists(nn[3])) { if (identical(get(nn[3], envir = .GlobalEnv), object[[3]])) { STAR[3] <- TRUE } }
  if (exists(nn[4])) { if (identical(get(nn[4], envir = .GlobalEnv), object[[4]])) { STAR[4] <- TRUE } }
  cat('  Probability functions for the KDE are the following: \n \n')
  cat('      Density function:                  ', nn[1], ifelse(STAR[1], '*', ''), '\n',
       '     Distribution function:             ', nn[2], ifelse(STAR[2], '*', ''), '\n',
       '     Quantile function:                 ', nn[3], ifelse(STAR[3], '*', ''), '\n',
       '     Random generation function:        ', nn[4], ifelse(STAR[4], '*', ''), '\n', '\n')
  if (any(STAR)) {
    cat('  * This function is presently loaded in the global environment \n \n') } }


plot.kde <- function(object, digits = 6, n = 512, cut = 4,
                     fill.colour = 'dodgerblue', fill.color = fill.colour) {

  #Check object class
  if (!('kde' %in% class(object)))    stop('Error: This plot method is only used for objects of class \'kde\'')

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
  df   <- object$df
  bw   <- object$bandwidth
  data <- object$data
  dens <- object[[1]]
  band.est   <- object$bandwidth.est
  df.est     <- object$df.est
  value.name <- object$value.name

  #Compute plot values
  xmin <- min(data) - bw*cut
  xmax <- max(data) + bw*cut
  xx   <- xmin + (xmax-xmin)*((1:n)-1)/(n-1)
  yy   <- dens(xx)
  zz   <- rep(0, n)

  #Create the subtitle
  if (df.est)   { E1 <- 'Estimated degrees-of-freedom = ' } else { E1 <- 'Input degrees-of-freedom = ' }
  if (band.est) { E2 <- 'Estimated bandwidth = '          } else { E2 <- 'Input bandwidth = '          }
  SUBTITLE <- paste0( '(', E1, if (df == Inf) { '\U221E' } else { formatC(df, digits = digits, format = 'f') },
                     ', ', E2, formatC(bw, digits = digits, format = 'f'), ')')

  #Create the plot
  plot(xx, yy, type = 'l',
       xlab = value.name, ylab = 'Density')
  title(main = 'Kernel Density Estimator')
  mtext(side = 3, line = 0.25, cex = 0.8, SUBTITLE)
  polygon(c(xx, rev(xx)), c(zz, rev(yy)), col = fill.color)

  #Save the plot
  PLOT <- recordPlot()
  dev.off()

  #Return the plot
  PLOT }
