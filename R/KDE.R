#' Kernel Density Estimator
#'
#' \code{KDE} returns the probability function for the kernel density estimator
#'
#' The kernel density estimator for a set of input data is obtained by taking a mixture distribution
#' consisting of a (possibly weighted) combination of kernels.  In this function we compute the KDE
#' using the kernel of the T-distribution with a specified value for the degrees-of-freedom (by default
#' the function uses infinite degrees-of-freedom giving the kernel of the normal distribution).  Given a
#' set of input data we optimise the bandwidth via leave-one-out maximum-likelihood estimation (LOO-MLE)
#' to generate the KDE.  The output of the function is a list of class \code{kde} that contains the
#' probability functions for the KDE and associated information.  The output objectcan be plotted to
#' show the density function for the KDE.
#'
#' Note: The function has an option \code{to.environment} to allow the user to load the probability
#' functions to the global environment.  If this is set to \code{TRUE} then the probability functions
#' are loaded to the global environment in addition to appearing as elements of the output; there is
#' a message informing the user if existing objects in the global environment were overwritten.
#'
#' @usage \code{KDE}
#' @param data Input data for the kernel density estimator (a numeric vector)
#' @param na.rm Logical; if \code{TRUE} the function will remove \code{NA} values from the data
#' @param weights Weights for the kernel density estimator (a numeric vector with the same length as the data)
#' @param bandwidth Input value for the bandwidth; if left unspecified it is estimated using LOO-MLE
#' @param df Degrees-of-freedom for the T-distribution (a positive numeric value)
#' @param density.name Name of the KDE distribution used for naming of the probability functions (a character string)
#' @param to.environment Logical; if \code{TRUE} the probability functions are attached to the global environment
#' @return A \code{kde} object containing the probability functions for the kernel density estimator

KDE <- function (data, na.rm = FALSE, weights = NULL,
                 bandwidth = NULL, df = Inf,
                 density.name = 'kde', to.environment = FALSE) {

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

  #Check other inputs
  if (!is.null(bandwidth)) {
    BAND.EST <- FALSE
    if (!is.numeric(bandwidth))           stop('Error: Input bandwidth should be a numeric value')
    if (length(bandwidth) != 1)           stop('Error: Input bandwidth should be a single numeric value')
    if (min(bandwidth) <= 0)              stop('Error: Input bandwidth should be a positive value') } else {
    BAND.EST <- TRUE }
  if (!is.numeric(df))                    stop('Error: Input df should be a numeric value')
  if (length(df) != 1)                    stop('Error: Input df should be a single numeric value')
  if (min(df) <= 0)                       stop('Error: Input df should be a positive value')
  if (!is.character(density.name))        stop('Error: Input density.name should be a character value')
  if (length(density.name) != 1)          stop('Error: Input density.name should be a single character value')
  if (!is.logical(to.environment))        stop('Error: Input to.environment should be a logical value')
  if (length(to.environment) != 1)        stop('Error: Input to.environment should be a single logical value')

  #Finalise the data and weights
  if (na.rm) {
    DATA <- data[!is.na(data)]
    WEIGHTS <- WEIGHTS[!is.na(data)] } else {
      DATA <- data }
  if (any(is.na(DATA))) {
    stop('Error: NA values in the data; set na.rm = TRUE if you want to remove these') }
  k <- length(DATA)
  if (!WEIGHTED) { WEIGHTS <- rep(1/k, k) } else { WEIGHTS <- weights/sum(weights) }

  #Estimate bandwidth via MLE (if not given)
  if (BAND.EST) {
    NEGLOGLIKE <- function(p) {
      VEC <- rep(0, k)
      for (i in 1:k) {
        DDD <- dt((DATA[-i]-DATA[i])*exp(-p), df = df, log = TRUE) - p + log(WEIGHTS[-i])
        VEC[i] <- matrixStats::logSumExp(DDD) }
      -sum(VEC) }
    NLM <- nlm(f = NEGLOGLIKE, p = 0)
    BAND <- exp(NLM$estimate)
  } else {
    BAND <- bandwidth }

  #Create computational components for KDE
  GEN.K       <- paste('k <-', length(DATA))
  GEN.DATA    <- paste('means <-', paste(deparse(DATA), collapse = ''))
  if (!WEIGHTED) {
    GEN.WEIGHTS <- NULL } else {
      GEN.WEIGHTS <- paste('weights <-', paste(deparse(WEIGHTS), collapse = '')) }

  #Create output list
  OUT <- vector(mode = 'list', length = 14)
  names(OUT)[1:4]  <- PROB.NAMES <- paste0(c('d', 'p', 'q', 'r'), density.name)
  names(OUT)[5:14] <- c('data.name', 'data', 'weighted', 'weights.name', 'weights',
                        'bandwidth', 'bandwidth.est', 'df', 'call', 'to.environment')
  OUT[[5]]   <- DATA.NAME
  OUT[[6]]   <- DATA
  OUT[[7]]   <- WEIGHTED
  OUT[[8]]   <- WEIGHTS.NAME
  OUT[[9]]   <- WEIGHTS
  OUT[[10]]  <- BAND
  OUT[[11]]  <- BAND.EST
  OUT[[12]]  <- df
  OUT[[13]]  <- deparse(CALL)
  OUT[[14]]  <- to.environment
  class(OUT) <- 'kde'

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
  OUT[[1]] <- generate_function(x = NA, bandwidth = BAND, df = df, log = FALSE, commands = COMM)

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
  OUT[[2]] <- generate_function(x = NA, bandwidth = BAND, df = df, lower.tail = TRUE, log.p = FALSE, commands = COMM)

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
  OUT[[3]] <- generate_function(p = NA, bandwidth = BAND, df = df, lower.tail = TRUE, log.p = FALSE, commands = COMM)

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
  OUT[[4]] <- generate_function(n = NA, bandwidth = BAND, df = df, commands = COMM)

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
  to.env       <- object$to.environment

  #Print title
  cat('\n  Kernel Density Estimator (KDE) \n \n')

  #Print description
  if (df == Inf) {
    cat('  Computed from', k, 'data points in the input', paste0('\'', data.name, '\''),
        'using the kernel of the normal distribution \n')
  } else {
    cat('  Computed from', k, 'data points in the input', paste0('\'', data.name, '\''),
        'using the kernel of the T-distribution with', df, 'degrees-of-freedom \n') }
  if (weighted) {
    cat('  KDE uses weights from input', paste0('\'', weights.name, '\''), '\n') }

  #Print bandwidth
  if (band.est) {
    cat('  Estimated bandwidth (LOO-MLE) =', formatC(bw, digits = digits, format = 'f'), ' \n \n')
  } else {
    cat('  Input bandwidth =', formatC(bw, digits = digits, format = 'f'), ' \n \n') }

  #Print information on probability functions
  if (to.env) {
    cat('  Probability functions for the KDE are the following (loaded to global environment): \n \n')
  } else {
    cat('  Probability functions for the KDE are the following: \n \n') }
  cat('      Density function:                  ', nn[1], '\n',
      '     Distribution function:             ', nn[2], '\n',
      '     Quantile function:                 ', nn[3], '\n',
      '     Random generation function:        ', nn[4], '\n', '\n') }


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
  band.est <- object$bandwidth.est

  #Compute plot values
  xmin <- min(data) - bw*cut
  xmax <- max(data) + bw*cut
  xx   <- xmin + (xmax-xmin)*((1:n)-1)/(n-1)
  yy   <- dens(xx)
  zz   <- rep(0, n)

  #Create the subtitle
  if (band.est) { EEE <- 'estimated bandwidth'} else {
    EEE <- 'input bandwidth' }
  if (df == Inf) {
    SUBTITLE <- paste0('(Kernel from the normal distribution; ', EEE, ' = ',
                       formatC(bw, digits = digits, format = 'f'), ')')
  } else {
    SUBTITLE <- paste0('(Kernel from the T-distribution with ', df, ' DF; ', EEE, ' = ',
                       formatC(bw, digits = digits, format = 'f'), ')') }

  #Create the plot
  plot(xx, yy, type = 'l',
       xlab = 'Value', ylab = 'Density')
  title(main = 'Kernel Density Estimator')
  mtext(side = 3, line = 0.25, cex = 0.8, SUBTITLE)
  polygon(c(xx, rev(xx)), c(zz, rev(yy)), col = fill.color)

  #Save the plot
  PLOT <- recordPlot()
  dev.off()

  #Return the plot
  PLOT }
